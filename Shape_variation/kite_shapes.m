%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure kite 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run ('C:/Users/danaz/OneDrive/Documents/MATLAB/gypsilab-master/addpathGypsilab.m')

%%%%%%%Test aire en fonction du nb de points 
% a = 0.65;
% c = 1.5;
% nb_points = linspace(10,1000,100);
% airetot = zeros(100,1);
% distance = @(x1,y1,x2,y2) sqrt((x1-x2)^2 + (y1-y2)^2);
% aire = @(AB,AC,BC) 0.5*AB*sqrt(AC^2-((BC^2 - AC^2 - AB^2)/(2*AB))^2);
% for n=1:100
% theta2 = linspace(0,2*pi, nb_points(n));
% x2 = cos(theta2) + a*cos(2*theta2) - a;
% y2 = c*sin(theta2);
% for i=1:nb_points(n)-1
%     xA = x2(i);
%     yA = y2(i);
%     xB = x2(i+1);
%     yB = y2(i+1);
%     AB = distance(xA,yA,xB,yB);
%     AC = distance(xA,yA,0,0);
%     BC = distance(xB,yB,0,0);
%     airetot(n) = airetot(n) + aire(AB,AC,BC);
% end
% end

distance = @(x1,y1,x2,y2) sqrt((x1-x2)^2 + (y1-y2)^2);
aire = @(AB,AC,BC) 0.5*AB*sqrt(AC^2-((BC^2 - AC^2 - AB^2)/(2*AB))^2);

%%%% aire de reference
airetot = 0;
theta2 = linspace(0,2*pi, 200);
x2 = cos(theta2) + a*cos(2*theta2) - a;
y2 = c*sin(theta2);
for i=1:200-1
    xA = x2(i);
    yA = y2(i);
    xB = x2(i+1);
    yB = y2(i+1);
    AB = distance(xA,yA,xB,yB);
    AC = distance(xA,yA,0,0);
    BC = distance(xB,yB,0,0);
    airetot = airetot + aire(AB,AC,BC);
end

%%%% variation des parametres 
Na = 10;
Nc = 10;
a = linspace(0.65,1.2,Na);
c = linspace(0,2,Nc);

Z = zeros(Na,Nc);
theta2 = linspace(0,2*pi, 200);
for ia=1:Na
    for ic=1:Nc
    x2 = cos(theta2) + a(ia)*cos(2*theta2) - a(ia);
    y2 = c(ic)*sin(theta2);
    for i=1:200-1
        xA = x2(i);
        yA = y2(i);
        xB = x2(i+1);
        yB = y2(i+1);
        AB = distance(xA,yA,xB,yB);
        AC = distance(xA,yA,0,0);
        BC = distance(xB,yB,0,0);
        Z(ia,ic) = Z(ia,ic) + aire(AB,AC,BC);
    end
    end
end
errZ = abs(Z-airetot);

%%%%%%%%%%variation a
N   = 1e3;
h     = (2*pi)/N;
theta = (0:h:(2*pi-h))';
vtx  = [cos(theta) + 0.65*cos(2*theta) - 0.65, 1.5*sin(theta), zeros(N,1)];

circle = 0.5*[cos(theta), sin(theta), zeros(N,1)];

figure
N   = 1e3;
n = 15;
cc = parula(n);%jet(n);
a = zeros(16,1);
a(1) = 0.65;
a(2:16) = linspace(0.7,1.2,15);%linspace(0.65,1.2,n);
zmap = linspace(min(a),max(a),length(cc));
for ia = 1:n
    color = interp1(zmap,cc,a(ia));
    vtx  = [cos(theta) + a(ia)*cos(2*theta) - a(ia) , 1.5*sin(theta), zeros(N,1)];
    plot(vtx(:,1),vtx(:,2),'color',color);%cc(ia,:));
    hold on
end
plot(circle(:,1), circle(:,2), 'color', [0,0,0])
caxis([min(a), max(a)])
hcb=colorbar;
hcb.Label.String = "Shape parameter";
hcb.Limits = [0.65 1.2];
title('Variation of the shape of the kite (constant surface)')


% figure
% N   = 1e3;
% n = 10;
% cc = jet(n);
% a = linspace(0.65,0.85,n);
% b = a;%linspace(0.65,0.85,n);
% c = @(x) sqrt(R^2 - (x-xc).^2);
% for ia = 1:n
%     vtx  = [cos(theta) + a(ia)*cos(2*theta) - b(ia) , c(a(ia))*sin(theta), zeros(N,1)];
%     plot(vtx(:,1),vtx(:,2),'color',cc(ia,:));
%     hold on
% end
% plot(circle(:,1), circle(:,2), 'color', [0,0,0])
% colorbar
% 
% if false
% figure
% aref = 0.65;
% b = linspace(0.1,1,100);
% for ib = 1:100
%     vtx  = [cos(theta) + aref*cos(2*theta) - b(ib), 1.5*sin(theta), zeros(N,1)];
%     plot(vtx(:,1),vtx(:,2),'color',cc(ib,:));
%     hold on
% end
% colorbar
% 
% figure
% aref = 0.65;
% bref = 0.65;
% c = linspace(0.1,0.3,100);
% for ic = 1:100
%     vtx  = [cos(theta) + aref*cos(2*theta) - bref, c(ic)*sin(theta), zeros(N,1)];
%     plot(vtx(:,1),vtx(:,2),'color',cc(ic,:));
%     hold on
% end
% colorbar
% end