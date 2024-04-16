%% -----------------------------------------------------------------------
% Local mixing index
% Solid index for solid-solid mixing
% Scale-dependent mixing
% Segregation of particles

% Filter field with Gaussian filter and then compute local mixing index
% Both solid mixing and accessed region indices used in this script

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

% Clip region with no particles for plotting
Xi = zeros(Ny,Nx);
for j = 1:Nx
    for i = 1:Ny
        if P(i,j) > 0.0001
            Xi(i,j) = xi(i,j);
        else
            Xi(i,j) = NaN;
        end
    end
end

% Plot solid fraction field
figure(2)
surf(x,y,Xi,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Solid fraction')
set(gca,'FontSize',14)

%% Define filter function

% Filter width
filt = 0.03; 

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

%Plot filter fucntion
figure(3)
surf(r,s,G,'EdgeColor', 'none')

%% Filter solid fraction and plot

% Filter solid fraction field with convolution 
xi_f = imfilter(xi,G,'conv','replicate');

% Clip region with no particles for plotting
Xi_f = zeros(Ny,Nx);
for j = 1:Nx
    for i = 1:Ny
        if P(i,j) > 0.0001
            Xi_f(i,j) = xi_f(i,j);
        else
            Xi_f(i,j) = NaN;
        end
    end
end

%Plot filtered solid fraction field
figure(3)
surf(x,y,Xi_f,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Filtered solid fraction function')
set(gca,'FontSize',14)

%% Compute solid mixing index

Xi_sum = 0;
v = ones(Ny,Nx);
vtot = 0; 
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            Xi_sum = Xi_sum + Xi_f(k,j)*v(k,j);
            vtot = vtot + v(k,j);
        end
    end
end
Xi_f_mean = Xi_sum/vtot;

M = zeros(Ny,Nx); % initialize

for j = 1:Nx
    for k = 1:Ny
        if (Xi_f(k,j) < Xi_f_mean)
            M(k,j) =  Xi_f(k,j)/(2*Xi_f_mean);
        else
            M(k,j) = 1 - 0.5*(1-Xi_f(k,j))/(1-Xi_f_mean);
        end
    end
end

% Plot solid mixing index
figure(4)
surf(x,y,M,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Mixing index')
set(gca,'FontSize',14)

%% Accessed region mixing index

A = zeros(Ny,Nx); % initialize

% Select accessed region range [Xi_f_min, Xi_f_max]
Xi_f_min = 0.3;
Xi_f_max = 0.7;
 
for j = 1:Nx
    for k = 1:Ny
        if Xi_f(k,j) < Xi_f_min
           A(k,j) = 0; 
        elseif (Xi_f(k,j)>= Xi_f_min && Xi_f(k,j) < Xi_f_mean)
            A(k,j) =  (Xi_f(k,j) - Xi_f_min)/(2*(Xi_f_mean - Xi_f_min));
        elseif (Xi_f(k,j) >= Xi_f_mean && Xi_f(k,j) < Xi_f_max)
            A(k,j) = 1 - (Xi_f_max - Xi_f(k,j))/(2*(Xi_f_max - Xi_f_mean));
        else
            A(k,j) = 1;
        end
    end
end

% Clip region with no particles
for j = 1:Nx
    for i = 1:Ny
        if P(i,j) > 0.0001
            A(i,j) = A(i,j);
        else
            A(i,j) = NaN;
        end
    end
end

figure(5)
surf(x,y,A,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Accessed region mixing index')
set(gca,'FontSize',14)