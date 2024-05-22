
###############################################
# Computes a Cauchy integral in a rectangle
################################################

% n : order of besselH
% l : order of the power in the Cauchy integral
% Re1, Re2, Im1, Im2 : bounds of the interval of Re(k), Im(k)
function integrale = Cauchy(l, n, Re1, Re2, Im1, Im2)

%discretisation
NpointsR = 1e2;
deltaR = abs(Re1-Re2)/NpointsR;
NpointsI = 1e1;
deltaI = abs(Im1-Im2)/NpointsI;
%z^l * H_n' / H_n
%%%%%%%%%%%%%%%%%%%%% DIrichlet Boundary Condition
%f = @(n,z) (z.^l).*(besselh(n-1, z)./besselh(n,z) - n./z); %Bessel H
%%%%%%%%%%%%%%%%%%%%% Impedance Boundary Condition
R = 1.3;
Eta = @(k) k./10;
DEta = 1/10;
b = @(n,k) -k.*besselh(n+1, k*R) + (1i*Eta(k)+ n/R).*besselh(n, k*R);
Db = @(n,k) R*k.*besselh(n+2, k*R) -(2*(n+1)+R*1i.*Eta(k)).*besselh(n+1, k*R) ...
+ (1i*DEta+ (1i*Eta(k)+n/R).*(n./k)).*besselh(n, k*R);
f = @(n,z) z.^l .* Db(n,z)./b(n,z);

%integrale sur 4 segments
segment1 = 1i*linspace(Im1, Im2, NpointsI) + Re1;
integrale = 1i*deltaI*sum(f(n,segment1));

segment2 = linspace(Re1, Re2, NpointsR) + 1i*Im2;
integrale = integrale + deltaR*sum(f(n,segment2));

segment3 = 1i*linspace(Im2, Im1, NpointsI) + Re2;
integrale = integrale - 1i*deltaI*sum(f(n,segment3));

segment4 = linspace(Re2, Re1, NpointsR) + 1i*Im1;
integrale = integrale - deltaR*sum(f(n,segment4));

end
