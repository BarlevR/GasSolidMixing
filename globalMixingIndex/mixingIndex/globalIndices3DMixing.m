%% -----------------------------------------------------------------------
% Global mixing index

% Mixing based on volume fraction fields

% Indices:
%   Rose
%   Lacey
%   Relative standard deviation (RSD)
%   Variance amongst bimodal bin counts (VBBC)
%   Mixing segragation index (MSI)

% CFD mixing of biomass and sand in fuidized bed on a non-uniform 3D mesh
% Volume of cells are needed

% The volume fraction of two fluidized bed reactor designs are given for
% each second, and is located in: 
%   ../../CFDData/biomassSandMixingFluidizedBed
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

load ../../CFDData/biomassSandMixingFluidizedBed/design1/alpha10.txt

alpha = alpha10;
p1 = alpha(:,1); % Sand
p2 = alpha(:,2); % Biomass
p = p1+p2;       % Total particle
V = alpha(:,4);  % Cell volume
N = length(p);

%% Compute averages

p_res = 0.0001; % for computations of regions that only have particles

Vsum = 0;
Ncells = 0;
p1sum = 0;
p2sum = 0;

for i = 1:N
    if p(i) > p_res
        p1sum = p1sum + p1(i)*V(i);
        p2sum = p2sum + p2(i)*V(i);
        Vsum = Vsum + V(i);
        Ncells = Ncells+1;
    end
end

%Vsum = sum(V);
p1_mean = p1sum/Vsum;
p2_mean = p2sum/Vsum;

%% Rose Index

% Difference from mean for variance and standard deviation calc
p2_d = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        p2_d(i) = p2(i)-p2_mean;
    else
        p2_d(i) = 0;
    end
end

sigma_m = sqrt(sum(p2_d.^2,'all')/Ncells);
sigma_um = sqrt(p2_mean*(1-p2_mean));
RI = 1 - sigma_m/sigma_um

%% Lacey Index

sig = sum(p2_d.^2,'all')/Ncells;
sigo = p2_mean*(1-p2_mean);
sigr = p2_mean*(1-p2_mean)/Ncells;
LI = (sigo-sig)/(sigo-sigr)


%% Relative standard deviation RSD

RSD = sigma_m/p2_mean

%% Variance amongst bimodal bin counts VBBC
Nt=1;
p_d = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        p_d(i) = (p1(i) - (p1_mean/p2_mean)*p2(i))^2;
    else
        p_d(i) = 0;
    end
end

sigma_wei = (sum(p_d(:))/Ncells)

% Note: This index requires that values are rescaled with respect to
% the maximum value of sigma_wei in time

%% Mixing segragation index

xi = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        xi(i) = p2(i)/(p1(i)+p2(i));
    else
        xi(i) = 0;
    end
end

xi_sum = 0;
Vsum = 0;
Ncells = 0;
for i = 1:N
    if p(i) > p_res
        xi_sum = xi_sum+ xi(i)*V(i);
        Vsum = Vsum + V(i);
        Ncells = Ncells+1;
    end
end
%Vsum = sum(V);
xi_mean = xi_sum/Vsum;

Xi = zeros(N,1);

for i = 1:N
    Xi(i) = xi(i)/xi_mean;
end

Xi_sum = 0;
Vsum = 0;
Ncells = 0;

for i = 1:N
    if p(i) > p_res
        Xi_sum = Xi_sum+ Xi(i)*V(i);
        Vsum = Vsum + V(i);
        Ncells = Ncells+1;
    end
end
Xi_mean = Xi_sum/Vsum;

Xi_d = zeros(N,1);
for i = 1:N
    if p(i) > p_res
        Xi_d(i) = (Xi(i) - Xi_mean)^2;
    else
        Xi_d(i) = 0;
    end
end

MSI = sqrt(sum(Xi_d(:))/(Ncells - 1))