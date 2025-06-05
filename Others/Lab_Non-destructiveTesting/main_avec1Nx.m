clear all;
close all;

load('dataSAFT_exp.mat');  

% 2. Ouverture du fichier
[Nt, Nel] = size(A);
d=0.5; %mm

t = linspace(0,Nt-1,Nt);
u = linspace(0,(Nel-1)*d,Nel);

figure(1)
imagesc(u,t,A);
title('Représentation des données');
xlabel('axe spatial u');
ylabel('axe temporel t');

figure(2)
ascan = A(:,1);
plot(t,ascan);
title('Ascan (colonne 1 de A)');
xlabel('temps');


% 3. Définition de la grille de reconstruction

Nx = 512; % nb de points
xmin = 0;
xmax = u(Nel) ; 
x = linspace(xmin,xmax,Nx); %Construction du vecteur de position x

dz = c /(2*Fs);
Nz = Nt;
zmin = 0;
zmax = (Nz-1)*dz;
z = linspace(zmin,zmax,Nz);

% On récupère l indice qui correspond à z=30
ind_z30=round(30/dz);


% 4. Reconstruction de l'image par une méthode temporelle

tic
O = zeros(Nz,Nx); % Initialisation de l'image O

for ix=1:1:Nx
    for iz=1:1:Nz
        somme = 0;
        for i=1:1:Nel
            r = sqrt(z(iz)^2 +(x(ix)-u(i))^2);
            tau = 2*r / c ; 
            tau_round = round(tau*Fs);
            ind = tau_round + 1;
            if (ind <= Nt) 
                a = A(ind,i);
                somme = somme + a;
            end
        end
        O(iz,ix) = somme;
    end
end
toc

figure(3)
imagesc(x,z,O);
colorbar;
axis equal; 
title('Image reconstruite avec Nx=512');
xlabel('x (mm)');
ylabel('z (mm)');


% 5. Post-traitement de l'image

O2 = abs(hilbert(O));
O2 = O2./max(O2);

figure(4)
imagesc(x,z,O2);
colorbar;
axis equal; 
title('Image reconstruite post-traitement (Nx=512)');
xlabel('x (mm)');
ylabel('z (mm)');

z30_MAT=O2(ind_z30,:);

% MEX
tic
O3 = MEX_SAFT(c,Fs,Nt,Nel,u,Nx,Nz,x,z,A);
O3 = reshape(O3,Nz,Nx);

figure(5)
imagesc(x,z,O3);
colorbar;
axis equal; 
title('Image reconstruite via MEX  SAFT');
xlabel('x (mm)');
ylabel('z (mm)');


O3 = abs(hilbert(O3));
O3 = O3./max(abs(O3));
toc

z30_MEX=O3(ind_z30,:);

figure(6)
imagesc(x,z,O3);
colorbar;
axis equal; 
title('Image reconstruite via MEX  SAFT');
xlabel('x (mm)');
ylabel('z (mm)');


figure(7)
subplot(1,2,1);
plot(x,z30_MAT);
title('Ligne pour z=30mm pour Matlab');
subplot(1,2,2);
plot(x,z30_MEX);
title('Ligne pour z=30mm pour MEX');


er=z30_MEX-z30_MAT;
figure(8)
plot(x,er);
title('Différence MEX-MAT pour la ligne z=30mm');