%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Influence of the variation of eta on a pole of the kite
%  (computation of a pole for different values of the impedance parameter eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')

eta_nb = 20;%50;
eta = linspace(6,50,eta_nb);
kstart = 0.4022 - 0.9522i; %1.18 -0.63i;
Nkre = 10;%10;
Nkim = 10;%10;
Rek = zeros(Nkre,1);
Imk = zeros(Nkim,1);
space = 0.05;
SavingResNorm = zeros(Nkre,Nkim,eta_nb);
SavingPoles = zeros(eta_nb,1);

N   = 1e3;
h     = (2*pi)/N;
theta = (0:h:(2*pi-h))';
vtx  = [cos(theta) + 0.65*cos(2*theta) - 0.65, 1.5*sin(theta), zeros(N,1)];
elt  = [[(2:N)';1] (1:N)'];
meshD = msh(vtx,elt);
% Exterior
Rz = 1; %Radius of the circle of exterior points z
Nz = 20; %Number of points z
% -------------------Interior circle
RC = 3/4;
modk_max = sqrt(max(Rek)^2+max(Imk)^2);     %module max of k
nf_max = ceil(modk_max*RC+4);               %greatest number of Fourier coefficents
Nint = 10*nf_max;                           %number of points on the interior circle
meshC = mshCircle(Nint,RC);

for e=1:eta_nb
    Rek = linspace(real(kstart)-space,real(kstart)+space,Nkre);
    Imk = linspace(imag(kstart)-space,imag(kstart)+space, Nkim);
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
            FourierN = MatrixF_impedance_var(nf, meshD, meshC, k, eta(e));
            T(1, (i-1)*Nkim+j) = toc;
            tic
            SavingResNorm(i,j,e) = LSM_F(FourierN,meshD,Nz,k, RC);
            T(2, (i-1)*Nkim+j) = toc;
        end
    end
    tend = toc(tstart);
    disp(['total time : ', num2str(tend)]);
    disp(['average time for matrix N : ', num2str(sum(T(1,:))/(Nkre*Nkim))]);
    disp(['average time for Inversion : ', num2str(sum(T(2,:))/(Nkre*Nkim))]);
    %Determine the coordiantes of the maximum
    [M, I] = max(SavingResNorm(:,:,e));
    [minim, jj] = max(M);
    SavingPoles(e) = Rek(I(jj))+1i*Imk(jj);
    kstart = SavingPoles(e);
end
save('impedance_variation_6_50.mat','SavingPoles','SavingResNorm');
