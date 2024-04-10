%% -----------------------------------------------------------------------
% Local mixing index
% Richness index for gas-solid mixing
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

% Plot of particle phase 1
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

% Plot of particle phase 2
figure(2)
surf(x,y,P2,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
xlabel('X, m')
ylabel('Y, m')
title('Volume fraction Phase 2')
colorbar
colormap(jet(1024))
clim([0 0.45])
set(gca,'FontSize',14)

% Plot of total particle phase
figure(3)
surf(x,y,P,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
xlabel('X, m')
ylabel('Y, m')
title('Total particle volume fraction')
colorbar
colormap(jet(1024))
clim([0 0.65])
set(gca,'FontSize',14)


%% Calculate richness index

% Compute volume average of particle volume fraciton phase 
% for the region where particles exist

% Initialize 
p_sum = 0;  
v = ones(Ny,Nx); 
vtot = 0; 

% Calculate volume sum
for j = 1:Nx
    for k = 1:Ny
        if P(k,j) > 1e-4
            p_sum = p_sum + P(k,j)*v(k,j);
            vtot = vtot + v(k,j);
        end
    end
end

% Average
P_mean = p_sum/vtot;

% Initialize Richness field
R = zeros(Ny,Nx);
P1_max = max(max(P));

% Calculate richness index
for j = 1:Nx
    for k = 1:Ny
        if (P(k,j) < P_mean)
            R(k,j) =  P(k,j)/(2*P_mean);
        else
            R(k,j) = 0.75 + 0.25*((2*P(k,j)-(P1_max + P_mean))/(P1_max - P_mean));
        end
    end
end

% Plot richness index field
figure(4)
surf(x,y,R,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Richness Index')
set(gca,'FontSize',14)
