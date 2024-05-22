%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function LSM_F computes the norm of A^-1*phi_k
%with phi_k the Fourier coefficients of the 
%fundamental solution in R^2
%       A       : matrix
%       meshD   : mesh of the considered domain
%       Nz      : number of exterior points (taken
%                 at delta_z of the boundary dD)
%       k       : complex value of k
%       RC      : radius of the interior circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function solution = LSM_F(A,meshD,Nz,k, RC)

HEST = 1e-2;

[Ueta,Seta,Veta]=svd(A);
SIGMAN=diag(Seta);
SQ=SIGMAN.^2;
Us=Ueta';
%Unused 
% gammamin=HEST*min(SIGMAN);
% gammamax=HEST*max(SIGMAN);
% morozov = @(x,GUP) sum((GUP.*(x^2-(HEST^2*SQ)))./(SQ+x).^2 );
% options=optimset('Display','off');

np = max(size(A,1)); %number of fourier coef
rhs = zeros(np,1);

%norm function
deltaz = 2*pi/abs(k)*0.25;
Rx = @(X) sqrt( X(:,1).^2 + X(:,2).^2 + X(:,3).^2 ) + deltaz;
%Atan function
Atan = @(X) angle(X(:,1)+1i*X(:,2));
ignorm = 0;

s = size(meshD.vtx(:,1));
N = s(1);
for ix=1:N/Nz:N  
      zpoint=meshD.vtx(ix,:);
      for n=0:np-1
        rhs(n+1) = -0.5*1i*pi.*besselj(n,conj(k)*RC).*besselh(n,conj(k)*Rx(zpoint)).*exp(1i*n*Atan(zpoint));
      end
        sc=Us*rhs;
        GUP=real(conj(sc).*sc); 
        %
        % Compute the Tikhonov parameter by Morozov's principle
        % a est le alpha de Morozov
        %[a,fval,exitflag]=fzero(@(x) morozov(x,GUP),[gammamin, gammamax]);
        a = 0.0;
        % disp(['a=',num2str(a),'    fval =',num2str(fval)])
%         if a<0|exitflag<0
%           disp(['Problem: a= ',num2str(a),' Resetting a=0'])
%           disp(['Problem: exitflag= ',num2str(exitflag)])
%           a=0;
%         end
        % Save 1/||g|| and the regularization parameter
        % Valxi=(SIGMAN.*sc)./(SQ+a);
        Valg= ((SIGMAN.*sc)./(SQ+a)); % Veta*Valxi;
        Valg= (conj(Valg).*Valg);
        %ValFact= GUP./(SIGMAN + sqrt(a));
        ignorm = ignorm + sqrt(sum(Valg));
end  
solution = ignorm;


          






