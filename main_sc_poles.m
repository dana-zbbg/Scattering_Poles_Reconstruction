%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For a range of complex values of k
%   Computes the Fourier Matrix of the interior scattering operator
%   And the norm of its inverse taken on the Fourier coefficients of 
%   the fundamental solutions in R^2. 
%   Peaks at the values of the scattering poles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleaning
clear all
close all
clc

% Gypsilab path
%run('../../addpathGypsilab.m')
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')

% ------------------Parameters
%values of k (!!!true values of k, conjugate taken later)
Nkre = 10; 
Nkim = 5;
Rek = linspace(0.31,0.51,Nkre);
Imk = linspace(-0.45,-0.37, Nkim);
%boundary condition : 1 = Dirichlet, 2 = Robin
boundary = 2;
%domain : 1 = Circle , 2 = Kite, 3 = ellipse
domain = 1;
%Array containing the norm of || MatrixF^-1 H^(1)(x-z) e^in theta(z) ||
ResNorm = zeros(Nkre, Nkim);
% Exterior 
Rz = 1; %Radius of the circle of exterior points z
Nz = 20; %Number of points z
%mshDisk(Nz, Rz);%+[ 0 0 3]

% ----------------------Domain
if domain == 1 %Circle
    N   = 1e3;
    RD = 1.3; %Radius of main Domain
    meshD = mshCircle(N,RD);
elseif domain == 2 %Kite
    N   = 1e3;
    h     = (2*pi)/N;
    theta = (0:h:(2*pi-h))';
    vtx  = [cos(theta) + 0.65*cos(2*theta) - 0.65, 1.5*sin(theta), zeros(N,1)];
    elt  = [[(2:N)';1] (1:N)'];
    meshD = msh(vtx,elt);
elseif domain == 3 %Ellipse
    N   = 1e3;
    h     = (2*pi)/N;
    theta = (0:h:(2*pi-h))';
    a = 1;
    b = 0.9;
    vtx  = [a*cos(theta), b*sin(theta), zeros(N,1)];
    elt  = [[(2:N)';1] (1:N)'];
    meshD = msh(vtx,elt);
end

% -------------------Interior circle
RC = 3/4;
modk_max = sqrt(max(Rek)^2+max(Imk)^2);     %module max of k
nf_max = ceil(modk_max*RC+4);               %greatest number of Fourier coefficents
Nint = 10*nf_max;                           %number of points on the interior circle
meshC = mshCircle(Nint,RC);

%Time 
T = zeros(2, Nkre*Nkim);
%------------------------------Loop over k
tstart = tic;
for i=1:Nkre
    for j=1:Nkim
        disp(['boucle ', num2str((i-1)*Nkim+j), ' out of ', num2str(Nkre*Nkim)])
        tic
        k = Rek(i)-1i*Imk(j);
        %Number of Fourier coefficients
        nf = ceil(abs(k)*RC+4);
        % Frequency adjusted to maximum esge size
        stp  = meshD.stp;
        kmax = 1/stp(2);
        if (k > kmax)
            warning('Wave number is too high for mesh resolution')
        end
        %Matrix of interior scattering operator
        tic
        if boundary == 1
            FourierN = MatrixF(nf, meshD, meshC, k);
        elseif boundary == 2
             FourierN = MatrixF_impedance(nf, meshD, meshC, k);
        end
        T(1, (i-1)*Nkim+j) = toc;
        tic
        ResNorm(i,j) = LSM_F(FourierN,meshD,Nz,k, RC);
        T(2, (i-1)*Nkim+j) = toc;
    end
end
tend = toc(tstart);
disp(['total time : ', num2str(tend)]);
disp(['average time for matrix N : ', num2str(sum(T(1,:))/(Nkre*Nkim))]);
disp(['average time for Inversion : ', num2str(sum(T(2,:))/(Nkre*Nkim))]);

%------------------Display
figure
surf(Rek,Imk,ResNorm');
xlabel('Re(k)');
ylabel('Im(k)');
zlabel('||R(k)||');
disp(ResNorm);
%Determine the coordiantes of the maximum
[M, I] = max(ResNorm);
[minim, jj] = max(M);
kres = Rek(I(jj))+1i*Imk(jj);

%%%%%%%% Color scale
%colors
col = zeros(256,3);
%Blue to Green
col(1:6,3) = 1;
col(7:20,2) = linspace(0,1,14);
col(7:20,3) = 1;
%Green toward red
col(21:34,1) = linspace(0,1,14);
col(21:34,2) = 1;
col(21:34,3) = linspace(1,0,14);
%
col(35:47,1) = 1;
col(35:47,2) = linspace(1,0,13);
%Red (peaks)
col(48:256,1) = linspace(1,0.5,209); 
% 
% figure
% surf(Rekcompl,Imkcompl,ResNorm_compl');
% title('Scattering poles for a kite with impedance boundary condition')
% xlabel('Re(k)');
% ylabel('Im(k)');
% zlabel('R(k)');
% colormap(col);
% colorbar

