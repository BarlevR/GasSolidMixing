%% -----------------------------------------------------------------------
% Global mixing index

% Solid mixing based no solid fraction fields

% Indices:
%   Rose
%   Lacey

% CFD segregation of particles on a 2D uniform mesh

% Authors: Barlev R. Nagawkar, Iowa State University

% Reference: 
%   Nagawkar BR, Kotrike VP, Passalacqua A, Subramaniam S. 
%   An index to characterize gas-solid and solid-solid mixing from average 
%   volume fraction fields. AIChE J. 2022;68(6):e17639. 
%   doi:10.1002/aic.17639
% ------------------------------------------------------------------------

clear
close all
clc

%% Load and compute solid fraction

load ../../CFDdata/particleSegregation/alpha20.dat

% Sort data 
D = sortrows(alpha20,[2 1]);

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

% Construct solid fraction
xi = zeros(Ny,Nx);
for j = 1:Nx
    for i = 1:Ny
        if P(i,j) > 0.0001
            xi(i,j) = P1(i,j)/(P1(i,j)+P2(i,j));
        else
            xi(i,j) = NaN;
        end
    end
end

% Plot solid fraction field
figure(2)
surf(x,y,xi,'EdgeColor', 'none')
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

%% Compute global index

% Initialize 
xi_sum = 0;
v = ones(Ny,Nx);
vtot = 0; 
ncells = 0;

% Calculate volume sum
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            xi_sum = xi_sum + xi(k,j)*v(k,j);
            vtot = vtot + v(k,j);
            ncells = ncells + 1;
        end
    end
end

% Average
xi_mean = xi_sum/vtot;

% Difference from mean for variance and standard deviation calc
xi_d = zeros(Ny,Nx);
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 0.0001
            xi_d(k,j) = xi(k,j) - xi_mean;
        else
            xi_d(k,j) = 0;
        end
    end
end

% Rose index
sigma_m = sqrt(sum(xi_d.^2,'all')/ncells);
sigma_um = sqrt(xi_mean*(1-xi_mean));
RI = 1 - sigma_m/sigma_um

% Lacey index (using variance)
sig = sum(xi_d.^2,'all')/ncells;
sigo = xi_mean*(1-xi_mean);
sigr = xi_mean*(1-xi_mean)/ncells;
Li = (sigo-sig)/(sigo-sigr)