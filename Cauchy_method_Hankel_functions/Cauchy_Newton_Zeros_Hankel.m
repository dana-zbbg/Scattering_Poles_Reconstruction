
################################
# Find zeros of Bessel Function
# using Cauchy integral
################################


%Borne sup ordre fonction de bessel
nbessel = 20;
%Domaine de recherche : Rectangle
nb_box_re = 10; %nombre de boite selon Re
nb_box_im = 10; %nombre de boite selon Im
RealParts = linspace(0.31,3.05,nb_box_re);
ImagParts = linspace(-0.96,-0.32,nb_box_im);
%poles
Ntotpoles = 0;
zerosHankel = zeros(nbessel+1,50);
nbzeros_per_n = zeros(nbessel, 1); %nombre de zeros pour chaque ordre

%mesh
N = 1e3;
R = 1;
meshD = mshCircle(N,R);

npoles_besseln_max = 0;
npoles_outbox = 0;
for n=0:nbessel
    nbpoles_besseln = 0;
    for re=1:nb_box_re-1
         Re1 = RealParts(re);
         Re2 = RealParts(re+1);
        for im=1:nb_box_im-1
            Im1 = ImagParts(im);
            Im2 = ImagParts(im+1);
            out_box = @(z) real(z)<Re1 || real(z)>Re2 || imag(z)<Im1 || imag(z)>Im2 ;
            %Nombre de poles dans la boite Re1,Re2,Im1,Im2
            Npoles = Cauchy(0,n, Re1, Re2, Im1, Im2)/(2*1i*pi);
            %zeros of bessel functions
            %Npoles = Cauchy_matrix(meshD,0,Re1,Re2,Im1,Im2);
            if abs(Npoles)> 0.4
                Nz = ceil(abs(Npoles));
                nbpoles_besseln = nbpoles_besseln + Nz;
                Ntotpoles = Ntotpoles+Nz;
                cauchy_sums = zeros(Nz,1);
                %Calcul de l'integrale de Cauchy pour autant d'ordre que de
                %zeros dans la boite
                for l=1:Nz
                    cauchy_sums(l) = -Cauchy(l,n,Re1, Re2, Im1, Im2)/(2*1i*pi);%Cauchy_matrix(meshD,l,Re1, Re2,Im1,Im2);%
                end
                %Calcul des valeurs des poles
                if Nz == 1
                    zerosHankel(n+1,nbpoles_besseln) = cauchy_sums(1);
                    disp([Nz, cauchy_sums(1)]);
                    if out_box(cauchy_sums(1))
                        disp('pole en dehors de la boite');
                        npoles_outbox = npoles_outbox + 1;
                    end
                elseif Nz == 2
                    discr = 4*cauchy_sums(1)^2- 8*(cauchy_sums(1)^2-cauchy_sums(2));
                    r1 = cauchy_sums(1)/2-sqrt(discr)/4;
                    r2 = cauchy_sums(1)/2+sqrt(discr)/4;
                    if out_box(r1)&& out_box(r2)
                        disp(['racines hors du domaine ', num2str(Re1), num2str(Re2), num2str(Im1), num2str(Im2)]);
                        npoles_outbox = npoles_outbox + 1;
                    elseif out_box(r1)&& not(out_box(r2))
                        zerosHankel(n+1,nbpoles_besseln-1) = r2;
                        zerosHankel(n+1,nbpoles_besseln) = cauchy_sums(1) - r2;
                    elseif out_box(r2)&& not(out_box(r1))
                        zerosHankel(n+1,nbpoles_besseln-1) = r1;
                        zerosHankel(n+1,nbpoles_besseln) = cauchy_sums(1) - r1;
                    else
                        disp('les deux racines sont dans le domaine...');
                        disp([r1,r2]);
                    end
                    if out_box(zerosHankel(n+1,nbpoles_besseln))
                        npoles_outbox = npoles_outbox + 1;
                    end
                    disp([Nz, zerosHankel(n+1,nbpoles_besseln), zerosHankel(n+1,nbpoles_besseln-1)])
                else
                    disp('Attention Nz > 2');
                end
            end
        end
    end
    nbzeros_per_n(n+1) = nbpoles_besseln;
    npoles_besseln_max = max(npoles_besseln_max, nbpoles_besseln);
end
disp(['Nombre total de poles : ', num2str(Ntotpoles)]);
zerosHankel2 = zerosHankel(:,1:npoles_besseln_max);

poles = zeros(Ntotpoles,1);
ind = 1;
for i=1:nbessel
    for j=1:npoles_besseln_max
        if abs(zerosHankel2(i,j))>1e-6
            poles(ind) = zerosHankel2(i,j);
            ind = ind + 1;
        end
    end
end
scatter(real(poles), imag(poles));


%%%%%%%%%%%%%%%%%%%%%%%% Tracé Dirichlet
image = false;
if image
    N = 1000;
    X = linspace(1,1.5,N);
    Y = linspace(-1.8, -1.3, N);
    Z = zeros(N, N);
    for i=1:N
        for j=1:N
            Z(i,j) = abs(besselh(3,X(i)+1i*Y(j)));
        end
    end

    %surf(X,Y,Z);
    [M, I] = min(Z);
    [minim, jj] = min(M);
    kt = X(I(jj))+1i*Y(jj);
end
%%%%%%%%%%%%%%%%% Tracé Impedance
image2 = false;
if image2
    R = 1.3;
    Eta = @(k) k./10;
    b =  @(n,k) -k.*besselh(n+1, k*R) + (1i*Eta(k)+ n/R).*besselh(n, k*R);
    N = 1000;
    X = linspace(0.1,2,N);
    Y = linspace(-2, -0.1, N);
    Z = zeros(N, N);
    for i=1:N
        for j=1:N
            k = (X(i)+1i*Y(j));
            Z(i,j) = abs(b(2,k));
        end
    end
end
