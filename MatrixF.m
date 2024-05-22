%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fourier Matrix of the interior scattering operator 
%  for a Dirichlet problem.
%       nf      : number of Fourier coefficients considered.
%       meshD   : mesh of the domain considered
%       meshC   : mesh of the interior domain
%       k       : value of k (should be such that Im(k) >0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MatrixF = MatrixF(nf, meshD, meshC, k)
 
typ = 'P1';
gss = 3;
tol = 0;
RC = norm(meshC.vtx(1,:));
%norm function
Rx = @(X) sqrt( X(:,1).^2 + X(:,2).^2 + X(:,3).^2 );
%Atan function
Atan = @(X) angle(X(:,1)+1i*X(:,2));
%Initialiaztion of the Matrix
FourierN = zeros(nf+1,nf+1);
% Green kernels
Gxy = @(X,Y) femGreenKernel(X,Y,'[H0(kr)]',k);

% Domain
sigma = dom(meshD,gss);

% Finite elements
u = fem(meshD,typ);
v = fem(meshD,typ);

% Finite element boundary operator --> \int_Sx \int_Sy psi(x)' G(x,y) psi(y) dx dy 
S = (1i/4) .* integral(sigma,sigma,u,Gxy,v,tol);
% Regularization
Sr  = -1/(2*pi) .* regularize(sigma,sigma,u,'[log(r)]',v);
LHS = S + Sr;

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
    PW = @(X) 0.5*1i*pi*besselj(n,conj(k)*RC)*besselh(n,conj(k)*Rx(X)).*exp(1i*n*Atan(X));
    RHS(:,n+1) = - integral(sigma,u,PW);
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
