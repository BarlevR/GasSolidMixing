%% -----------------------------------------------------------------------
% Global mixing index analysis
% Segregation of particles

% Filter solid fraction field with Gaussian filter and check histogram 

% Author: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear 
close all
clc

%% Load and plot the volume fraction fields

load ../../../CFDdata/particleSegregation/alpha30.dat

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

%% Fully mixed and segregated states

% Load fully mixed state (initial condition)
load ../../../CFDdata/particleSegregation/alpha0.dat

% Sort data
Dm = sortrows(alpha0,[2 1]);

% Coordinates (x,y)
xm = Dm(1:180,1);
ym = Dm(1:180:end,2);

% Volume fractions
p1m = Dm(:,4);
p2m = Dm(:,5);
pm = Dm(:,6);

% Organize data for plotting in physical (x,y)
P1m = reshape(p1m,180,420)';
P2m = reshape(p2m,180,420)';
Pm = reshape(pm,180,420)';

Xim = P1m./Pm;
figure(2)
surf(xm,ym,Xim,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Solid fraction, perfectly mixed')
set(gca,'FontSize',14)

% Create the fully segregated stated with equal portions of phases
PUm1 = zeros(420,180);
PUm2 = zeros(420,180);
PUm1(1:168,:) = 0.61;
PUm2(169:336,:) = 0.61;
PUm = PUm1+PUm2;

Xium = PUm1./PUm;
figure(3)
surf(xm,ym,Xium,'EdgeColor', 'none')
view(0,90)
axis([-0.06 0.06 0 0.5])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Solid fraction, fully segregated')
set(gca,'FontSize',14)

%% Histogram

% Unfiltered 
figure(4)
h1 = histogram(xi,100);
ycoords1 = h1.Values;
edgecoords = h1.BinEdges;
xcoords1 = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Unfiltered solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)
axis([0 1 0 11000])
ax=gca; 
ax.YAxis.Exponent = 3;

% Filtered
figure(5)
h2 = histogram(xi_f,100);
ycoords2 = h2.Values;
edgecoords = h2.BinEdges;
xcoords2 = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Filtered solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)
axis([0 1 0 11000])
ax=gca; 
ax.YAxis.Exponent = 3;

% Mixed state
figure(6)
hm = histogram(Xim,100);
ycoords_m = hm.Values;
edgecoords = hm.BinEdges;
xcoords_m = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)


% Segregated state
figure(7)
hum = histogram(Xium,100);
ycoords_um = hum.Values;
edgecoords = hum.BinEdges;
xcoords_um = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)



