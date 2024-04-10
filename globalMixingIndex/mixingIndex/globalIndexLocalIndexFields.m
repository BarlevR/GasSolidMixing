%% -----------------------------------------------------------------------
% Global mixing index

% Mixing based on local mixing index fields

% Index:
%   Rose

% CFD mixing of biomass and sand in fuidized bed on a non-uniform 3D mesh
% Volume of cells are needed

% The volume fraction of two fluidized bed reactor designs are given for
% each second, and is located in: 
%   ../../CFDdata/biomassSandMixingFluidizedBed
% The files have the same names so change design name in script as needed
%   (design1 and design2)

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

%% load and set data

load ../../CFDdata/biomassSandMixingFluidizedBed/design2/alpha10.txt

alpha = alpha10;
p1 = alpha(:,1); % Sand
p2 = alpha(:,2); % Biomass
p = p1+p2;       % Total particle
V = alpha(:,4);  % Cell volume
N = length(p);

%% Solid fraction and mean

xi = zeros(size(p));
p_res = 0.0001; 
for i = 1:N
    if p(i) > p_res
        xi(i) = p2(i)/p(i);
    else
        xi(i) = 0;
    end
end

xi_sum = 0;
Vsum = 0;
Ncells = 0;
p2sum = 0;

for i = 1:N
    if p(i) > p_res
        xi_sum = xi_sum + xi(i)*V(i);
        p2sum = p2sum + p2(i)*V(i);
        Vsum = Vsum + V(i);
        Ncells = Ncells+1;
    end
end

%Vsum = sum(V);
xi_mean = xi_sum/Vsum;
p2_mean = p2sum/Vsum;

%% Global Rose index on solid fraction

xi_d = zeros(N,1);
% Difference from mean for standard deviation calc
for i = 1:N
    if p(i) > p_res
        xi_d(i) = xi(i)-xi_mean;
    else
        xi_d(i) = 0;
    end
end

sigma_m = sqrt(sum(xi_d.^2,'all')/Ncells);
sigma_um = sqrt(xi_mean*(1-xi_mean));
GRI = 1 - sigma_m/sigma_um

%% Local mixing index

M = zeros(N,1); % initialize

for i = 1:N
    if (xi(i) < xi_mean)
        M(i) =  xi(i)/(2*xi_mean);
    else
        M(i) = 1 - 0.5*(1-xi(i))/(1-xi_mean);
    end
end

%% Global Rose on local mixing index

% Calculate average of local index
M_sum = 0;
Vsum = 0;
for i = 1:N
    if p(i) > p_res
        M_sum = M_sum + M(i)*V(i);
        Vsum = Vsum + V(i);
    end
end
M_mean = M_sum/Vsum;

% Difference from mean for standard deviation calc
M_d = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        M_d(i) = M(i)-M_mean;
    else
        M_d(i) = 0;
    end
end

Msigma_m = sqrt(sum(M_d.^2,'all')/Ncells);
Msigma_um = sqrt(M_mean*(1-M_mean));
MGRI = 1 - Msigma_m/Msigma_um

%% Accessed region local index

A = zeros(N,1); % initialize

% Select accessed region range [a_min, a_max]
a_min = 0;
a_max = 0.5;%max(max(p1));
a_mean = xi_mean;
a = xi;
% Compute accessed region mixing index field
for i = 1:N
    if a(i) < a_min
       A(i) = 0; 
    elseif (a(i)>= a_min && a(i) < a_mean)
       A(i) =  (a(i) - a_min)/(2*(a_mean - a_min));
    elseif (a(i) >= a_mean && a(i) < a_max)
       A(i) = 1 - (a_max - a(i))/(2*(a_max - a_mean));
    else
       A(i) = 1;
    end
end

%% Global Rose index based on local accessed region index

% Calculate average of accessed region field
A_sum = 0;
Vsum = 0;
for i = 1:N
    if p(i) > p_res
        A_sum = A_sum + A(i)*V(i);
        Vsum = Vsum + V(i);
    end
end
A_mean = A_sum/Vsum;

% Difference from mean for standard deviation calc
A_d = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        A_d(i) = A(i)-A_mean;
    else
        A_d(i) = 0;
    end
end

Asigma_m = sqrt(sum(A_d.^2,'all')/Ncells);
Asigma_um = sqrt(A_mean*(1-A_mean));
AGRI = 1 - Asigma_m/Asigma_um