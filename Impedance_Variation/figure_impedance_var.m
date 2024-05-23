%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   figure impedance var kite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')
data = load('impedance_variation_savings.mat');
matrix = data.SavingResNorm;
poles = data.SavingPoles;
eta_nb = 70;
eta = log(data.eta);
cc1 = parula(eta_nb);
% cc = zeros(eta_nb,3);
% step1 = 5;
% cc(1:step1,1:2) = zeros(step1,2);
% cc(1:step1,3) = linspace(0.556,1,step1);
% 
% step2 = 12;
% cc(step1+1:step1+step2,1) = zeros(step2,1);
% cc(step1+1:step1+step2,2) = linspace(0.0556,1,step2);
% cc(step1+1:step1+step2,3) = ones(step2,1);
% 
% step3 = 12;
% cc(step2+1:step2+step3,1) = linspace(0.0556,1,step3);
% cc(step2+1:step2+step3,2) = linspace(1,1,step3);
% cc(step2+1:step2+step3,3) = linspace(1,0,step3);
% 
% step4 = 31;
% cc(step3+1:step3+step4,1) = linspace(1,1,step4);
% cc(step3+1:step3+step4,2) = linspace(1,0,step4);
% cc(step3+1:step3+step4,3) = linspace(0,0,step4);
% 
% step5 = 10;
% cc(step4+1:step4+step5,1) = linspace(1,0.611,step5);
% cc(step4+1:step4+step5,2) = linspace(0,0,step5);
% cc(step4+1:step4+step5,3) = linspace(0,0,step5);

figure
scatter(real(poles),imag(poles), [],eta,'filled');
hold on 
Dir1 = [0.39];
Dir2 =  [-1.164];
scatter(Dir1,Dir2,'red');text(Dir1,Dir2,{'  Dirichlet Pole'},'FontSize',11);
hold on
Neu1 = [1.22];
Neu2 = [-0.64];
scatter(Neu1,Neu2,'b');text(Neu1,Neu2,{' Neumann Pole'},'FontSize',11)
axis([0.2 1.6 -1.2 -0.6])
hcb=colorbar;
hcb.Label.String = "Impedance parameter eta";
hcb.FontSize = 11;
hcb.TickLabels = num2cell([0.1    0.4    1.0000    2.7    7.4   20]) ; 
%hcb.Limits = [0.1 50];
colormap turbo%jet
xlabel('Re(k)')
ylabel('Im(k)')
title('Influence of the variation of the impedance on a scattering pole', 'FontSize',11);

% 
% figure
% scatter(exp(eta),real(poles),50,'b');
% hold on 
% scatter(exp(eta),imag(poles),50,'red');