%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   figure shape var kite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')

data4 = load('shape_variation4.mat');
poles4 = data4.SavingPoles;
matrix4 = data4.SavingResNorm;
shape4 = data4.shape1;
shape_nb4 = 15;
cc = parula(shape_nb4);%jet(shape_nb4);

figure
scatter(real(poles4), imag(poles4), [], shape4);
hold on
theta = linspace(0,2*pi);
radius = 0.0034;
xc = 0;
yc = 0;
for point=1:shape_nb4
    xc = real(poles4(point));
    yc = imag(poles4(point));
    plot(xc+radius*cos(theta), yc+radius*sin(theta), 'color',cc(point,:));
    hold on
end

hcb=colorbar;
hcb.Label.String = "Shape parameter";
colormap parula%jet
xlabel('Re(k)')
ylabel('Im(k)')
title('Influence of the variation of the shape on a scattering pole');

% for i=2:15
%     figure
%     surf(matrix(:,:,i))
% end
% for i=2:15
%     figure
%     surf(matrix3(:,:,i))
% end

