%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Influence of the variation of the shape of the kite on a pole
% (computation of a pole of the kite for differente shapes)
% the shape of the kite can be parametrize and we use these parametrization
% to change the shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')
%run ('/home/dz289/Documents/gypsilab/gypsilab-master/addpathGypsilab.m')

eta = 0.1;
shape_nb = 15;
shape1 = linspace(0.7,1.2,shape_nb);
data = load('shape_variation3.mat');
%kstarts = data.SavingPoles;
%kstart = kstats(1);
kstart = 1.18 -0.63i;
Nkre = 8;
Nkim = 8;
Rek = zeros(Nkre,1);
Imk = zeros(Nkim,1);
spaceim = 0.012;
spacere = 0.012;
SavingResNorm = zeros(Nkre,Nkim,shape_nb);
SavingPoles = zeros(shape_nb,1);

%Constants of the domain
N   = 1e3;
h     = (2*pi)/N;
theta = (0:h:(2*pi-h))';
elt  = [[(2:N)';1] (1:N)'];

% Exterior
Rz = 1; %Radius of the circle of exterior points z
Nz = 20; %Number of points z
% -------------------Interior circle
RC = 3/4;
modk_max = sqrt(max(Rek)^2+max(Imk)^2);     %module max of k
nf_max = ceil(modk_max*RC+4);               %greatest number of Fourier coefficents
Nint = 10*nf_max;                           %number of points on the interior circle
meshC = mshCircle(Nint,RC);

for s=1:shape_nb
    %kstart = kstarts(s);
    %-----------------------Change domain
    vtx  = [cos(theta) + shape1(s)*cos(2*theta) - shape1(s), 1.5*sin(theta), zeros(N,1)];
    meshD = msh(vtx,elt);

    Rek = linspace(real(kstart)-spacere,real(kstart)+spacere,Nkre);
    Imk = linspace(imag(kstart)-spaceim,imag(kstart)+spaceim, Nkim);
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
            FourierN = MatrixF_impedance_var(nf, meshD, meshC, k, eta);
            T(1, (i-1)*Nkim+j) = toc;
            tic
            SavingResNorm(i,j,s) = LSM_F(FourierN,meshD,Nz,k, RC);
            T(2, (i-1)*Nkim+j) = toc;
        end
    end
    tend = toc(tstart);
    disp(['total time : ', num2str(tend)]);
    disp(['average time for matrix N : ', num2str(sum(T(1,:))/(Nkre*Nkim))]);
    disp(['average time for Inversion : ', num2str(sum(T(2,:))/(Nkre*Nkim))]);
    %Determine the coordiantes of the maximum
    [M, I] = max(SavingResNorm(:,:,s));
    [minim, jj] = max(M);
    SavingPoles(s) = Rek(I(jj))+1i*Imk(jj);
    kstart = SavingPoles(s);
end
save('shape_variation4.mat','SavingPoles','SavingResNorm');
