clear all;
close all;

load('dataSAFT_exp.mat');  

Nx_tab = [64 , 128 , 256 , 512 , 1024 , 2048];

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

% Nx = 512; % nb de points
% xmin = 0;
% xmax = u(Nel) ; 
% x = linspace(xmin,xmax,Nx); %Construction du vecteur de position x
% 
% dz = c /(2*Fs);
% Nz = Nt;
% zmin = 0;
% zmax = (Nz-1)*dz;
% z = linspace(zmin,zmax,Nz);


% 4. Reconstruction de l'image par une méthode temporelle
temps_mat = zeros([1,6]);

for i_nx=1:1:6
    Nx=Nx_tab(i_nx);
    %défininition de la grille de reconstruction
    xmin = 0;
    xmax = u(Nel) ; 
    x = linspace(xmin,xmax,Nx); %Construction du vecteur de position x
    
    dz = c /(2*Fs);
    Nz = Nt;
    zmin = 0;
    zmax = (Nz-1)*dz;
    z = linspace(zmin,zmax,Nz);
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
    %Post-traitement de l'image
    O2 = abs(hilbert(O));
    O2 = O2./max(O2);
    temps_mat(i_nx)=toc;
end

figure(3)
plot(Nx_tab,temps_mat);
title('Temps de calcul en fonction de Nx via Matlab');
ylabel('Temps (s)');
xlabel('Nx');

% figure(3)
% imagesc(x,z,O);
% colorbar;
% axis equal; 
% title('Image reconstruite avec Nx=512');
% xlabel('x (mm)');
% ylabel('z (mm)');


% 5. Post-traitement de l'image

% O2 = abs(hilbert(O));
% O2 = O2./max(O2);

% figure(4)
% imagesc(x,z,O2);
% colorbar;
% axis equal; 
% title('Image reconstruite post-traitement (Nx=512)');
% xlabel('x (mm)');
% ylabel('z (mm)');


% MEX
temps_mex = zeros([1,6]);

for i_nx=1:1:6
    Nx=Nx_tab(i_nx);
    tic
    O3 = MEX_SAFT(c,Fs,Nt,Nel,u,Nx,Nz,x,z,A);
    O3 = reshape(O3,Nz,Nx);
    
    O3 = abs(hilbert(O3));
    O3 = O3./max(abs(O3));
    temps_mex(i_nx)=toc;
end

figure(4)
plot(Nx_tab,temps_mex);
title('Temps de calcul en fonction de Nx via Mex');
ylabel('Temps (s)');
xlabel('Nx');

disp(temps_mat./temps_mex);

% figure(5)
% imagesc(x,z,O3);
% colorbar;
% axis equal; 
% title('Image reconstruite via MEX  SAFT');
% xlabel('x (mm)');
% ylabel('z (mm)');