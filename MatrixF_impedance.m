%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fourier Matrix of the interior scattering operator 
%  for an Impedance boundary problem.
%       nf      : number of Fourier coefficients considered.
%       meshD   : mesh of the domain considered
%       meshC   : mesh of the interior domain
%       k       : value of k (should be such that Im(k) >0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MatrixF = MatrixF_impedance(nf, meshD, meshC, k)
typ = 'P1';
gss = 3;
tol = 0;
RC = norm(meshC.vtx(1,:));
%norm function
Rx = @(X) sqrt( X(:,1).^2 + X(:,2).^2 + X(:,3).^2 );
%Atan function
Atan = @(X) angle(X(:,1)+1i*X(:,2));
%impedance function
Eta = @(k) k/10;

FourierN = zeros(nf+1,nf+1);
% Green kernels
Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);
dyGxy{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
dyGxy{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
dyGxy{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);

const = @(X,n)  0.5*1i*pi*besselj(n,conj(k)*RC).*exp(1i*n*Atan(X));
dyGxyn{1} = @(X,n) const(X,n).*(-conj(k) * besselh(n+1,conj(k)*Rx(X))+ n./Rx(X).*besselh(n,conj(k)*Rx(X))) .* X(:,1) ./ Rx(X) ;
dyGxyn{2} = @(X,n) const(X,n).*(-conj(k) * besselh(n+1,conj(k)*Rx(X))+ n./Rx(X).*besselh(n,conj(k)*Rx(X))) .* X(:,2) ./ Rx(X) ;
dyGxyn{3} = @(X,n) const(X,n).*(-conj(k) * besselh(n+1,conj(k)*Rx(X))+ n./Rx(X).*besselh(n,conj(k)*Rx(X))) .* X(:,3) ./ Rx(X) ;
% Domain
sigma = dom(meshD,gss);

% Finite elements
u = fem(meshD,typ);
v = fem(meshD,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
S = (1i/4) .* integral(sigma,sigma,u,Gxy,v,tol);
% Regularization
Sr  = -1/(2*pi) .* regularize(sigma,sigma,u,'[log(r)]',v);
% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' dnyG(x,y) psi(y) dx dy 
D = (1i/4) .* integral(sigma,sigma,u,dyGxy,ntimes(v),tol);
%Regularization
Dr = -1/(2*pi) .* regularize(sigma,sigma,u,'grady[log(r)]',ntimes(v));
%id
Id = 0.5*integral(sigma,u,v);
%final operator
LHS = (D + Dr).'+ Id - 1i*Eta(k)*(S + Sr);

% LU factorization
[Lh,Uh] = lu(LHS);
% Finite element radiative operator --> \int_Sy G(x,y) psi(y) dy 
Sdom = 1i/4 .* integral(meshC.vtx,sigma,Gxy,v,tol);
% Regularization
Sreg = -1/(2*pi) .* regularize(meshC.vtx,sigma,'[log(r)]',v);
Sdom = Sdom + Sreg;
%Fourier tools
Nint = max(size(meshC.vtx(:,1)));
theta =  2*pi/Nint .* (1:Nint);
deltaTheta = 2*pi/Nint;
tnf = (0:nf)';
% Incident wave
N = max(size(meshD.vtx(:,1)));
RHS = zeros(N, nf);
for n=0:nf
    PW1{1} = @(X) dyGxyn{1}(X,n);
    PW1{2} = @(X) dyGxyn{2}(X,n);
    PW1{3} = @(X) dyGxyn{3}(X,n);
    PW2 = @(X) const(X,n).*1i.*Eta(conj(k)).*besselh(n,conj(k)*Rx(X));
    RHS(:,n+1) = - integral(sigma,ntimes(u),PW1) -integral(sigma,u,PW2);
end
RHS = conj(RHS);
% Solve linear system [S] * lambda = P0
lambda  = Uh \ (Lh \ RHS); % LHS \ RHS;
% Domain solution
Psca = Sdom * lambda;
Psca = conj(Psca);
FourierN(:,:) = exp(-1i*tnf*theta)*Psca*deltaTheta;
MatrixF = FourierN;
end




%LHS  = beta.*(-0.5*Id + (D+Dr).') - (H+Hr);

% Finite element incident wave trace --> \int_Sx psi(x) dnx(pw(x)) dx
%RHS = - integral(sigma,ntimes(u),gradxPW);