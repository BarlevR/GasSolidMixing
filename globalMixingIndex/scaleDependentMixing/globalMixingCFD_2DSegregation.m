%% -----------------------------------------------------------------------
% Global mixing index
% Accessed region solid index for solid-solid mixing
% Segregation of particles

% Filter field with Gaussian filter and then compute global mixing index
% on the filtered solid fraction fields

% Author: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear 
close all
clc

%% Load and plot the volume fraction fields

load ../../CFDdata/particleSegregation/alpha30.dat

% Sort data 
D = sortrows(alpha30,[2 1]);

% Coordinates (x,y)
x = D(1:180,1);
y = D(1:180:end,2);

% Volume fractions
p1 = D(:,4); % phase 1
p2 = D(:,5); % phase 2
p = D(:,6); % total volume fraction ()

% No. of cells in each diretion
Nx = length(x);
Ny = length(y);

% Organize data for plotting in physical (x,y)
P1 = reshape(p1,180,420)'; % phase 1
P2 = reshape(p2,180,420)'; % phase 2
P = reshape(p,180,420)'; % total volume fraction

% plot of particle phase 1
figure(1)
surf(x,y,P1,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
xlabel('X, m')
ylabel('Y, m')
title('Volume fraction Phase 1')
colorbar
colormap(jet(1024))
clim([0 0.45])
set(gca,'FontSize',14)

%% Calculate and plot solid fraciton field

% Solid fraction
xi = P1./(P1+P2);

for j = 1:Nx
    for i = 1:Ny
        if isnan(xi(i,j))
            xi(i,j) = 0;
        end
    end
end

%% Define filter function

% Filter width
filt = 0.04; 

% Coordinates for gauss filter
delr = x(2) - x(1);
dels = y(2) - y(1);
N = filt*5;
r = -1*N:delr:N;
s = -1*N:dels:N;

nx = length(r);
ny = length(s);

% Gauss filter
G = zeros(ny,nx);
Gsum = 0;
for j = 1:nx
    for k = 1:ny
        G(k,j) = (6/(pi*(filt^2)))*exp(-6*((r(j))^2+(s(k))^2)/(filt^2));
        Gsum = Gsum + G(k,j); 
    end
end
G = G/Gsum;

%% Filter solid fraction and plot

% Filter solid fraction field with convolution 
xi_f = imfilter(xi,G,'conv','replicate');

%% Compute solid mixing index

xi_sum = 0;
v = ones(Ny,Nx);
vtot = 0;
ncells = 0;
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            xi_sum = xi_sum + xi_f(k,j)*v(k,j);
            vtot = vtot + v(k,j);
            ncells = ncells + 1;
        end
    end
end
xi_f_mean = xi_sum/vtot;

M = zeros(Ny,Nx); % initialize
for j = 1:Nx
    for k = 1:Ny
        if (xi_f(k,j) < xi_f_mean)
            M(k,j) =  xi_f(k,j)/(2*xi_f_mean);
        else
            M(k,j) = 1 - 0.5*(1-xi_f(k,j))/(1-xi_f_mean);
        end
    end
end

%% Compute local accessed region mixing index

A = zeros(Ny,Nx); % initialize

% Select accessed region range [Xi_f_min, Xi_f_max]
xi_f_min = 0.3;
xi_f_max = 0.7;
 
for j = 1:Nx
    for k = 1:Ny
        if xi_f(k,j) < xi_f_min
           A(k,j) = 0; 
        elseif (xi_f(k,j)>= xi_f_min && xi_f(k,j) < xi_f_mean)
            A(k,j) =  (xi_f(k,j) - xi_f_min)/(2*(xi_f_mean - xi_f_min));
        elseif (xi_f(k,j) >= xi_f_mean && xi_f(k,j) < xi_f_max)
            A(k,j) = 1 - (xi_f_max - xi_f(k,j))/(2*(xi_f_max - xi_f_mean));
        else
            A(k,j) = 1;
        end
    end
end

%% Compute Rose index on filtered solid fraction

xi_f_d = zeros(Ny,Nx);
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 0.0001
            xi_f_d(k,j) = xi_f(k,j) - xi_f_mean;
        else
            xi_f_d(k,j) = 0;
        end
    end
end
sigma_m = sqrt(sum(xi_f_d.^2,'all')/ncells);
sigma_um = sqrt(xi_f_mean*(1-xi_f_mean));
GRI = 1 - sigma_m/sigma_um

%% Compute Rose index on local mixing index

M_sum = 0;
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            M_sum = M_sum + M(k,j)*v(k,j);
        end
    end
end
M_mean = M_sum/vtot;

M_d = zeros(Ny,Nx);
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 0.0001
            M_d(k,j) = M(k,j) - M_mean;
        else
            M_d(k,j) = 0;
        end
    end
end

Msigma_m = sqrt(sum(M_d.^2,'all')/ncells);
Msigma_um = sqrt(M_mean*(1-M_mean));
MGM = 1 - Msigma_m/Msigma_um


%% Compute Rose index from local accessed mixing
A_sum = 0;
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            A_sum = A_sum + A(k,j)*v(k,j);
        end
    end
end
A_mean = A_sum/vtot;

A_d = zeros(Ny,Nx);
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 0.0001
            A_d(k,j) = A(k,j) - A_mean;
        else
            A_d(k,j) = 0;
        end
    end
end

Asigma_m = sqrt(sum(A_d.^2,'all')/ncells);
Asigma_um = sqrt(A_mean*(1-A_mean));
AGM = 1 - Asigma_m/Asigma_um