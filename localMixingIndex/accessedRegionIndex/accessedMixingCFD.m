%% -----------------------------------------------------------------------
% Local mixing index
% Accessed region solid index for solid-solid mixing
% Segregation of particles

% Author: Barlev R. Nagawkar, Iowa State University

% Reference: 
%   Nagawkar BR, Kotrike VP, Passalacqua A, Subramaniam S. 
%   An index to characterize gas-solid and solid-solid mixing from average 
%   volume fraction fields. AIChE J. 2022;68(6):e17639. 
%   doi:10.1002/aic.17639
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

%% Calculate the accessed region mixing index

% Compute volume average of solid fraciton phase 
% for the region where particles exist

% Initialize 
xi_sum = 0;
v = ones(Ny,Nx);
vtot = 0; 

% Calculate volume sum
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            xi_sum = xi_sum + xi(k,j)*v(k,j);
            vtot = vtot + v(k,j);
        end
    end
end

% Average
xi_mean = xi_sum/vtot;

% Calculate index

% Pick accessed region
xi_min = 0.3;
xi_max = 0.7;

A = zeros(Ny,Nx); % initialize
for j = 1:Nx
    for k = 1:Ny
        if xi(k,j) < xi_min
           A(k,j) = 0; 
        elseif (xi(k,j)>= xi_min && xi(k,j) < xi_mean)
            A(k,j) =  (xi(k,j) - xi_min)/(2*(xi_mean - xi_min));
        elseif (xi(k,j) >= xi_mean && xi(k,j) < xi_max)
            A(k,j) = 1 - (xi_max - xi(k,j))/(2*(xi_max - xi_mean));
        else
            A(k,j) = 1;
        end
    end
end

% Clip region where particles don't exist
for j = 1:Nx
    for i = 1:Ny
        if P(i,j) > 0.0001
            A(i,j) = A(i,j);
        else
            A(i,j) = NaN;
        end
    end
end

% Plot accessed mixing index
figure(3)
surf(x,y,A,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Solid mixing index')
set(gca,'FontSize',14)