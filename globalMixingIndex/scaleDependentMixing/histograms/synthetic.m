%% -----------------------------------------------------------------------
% Global mixing index analysis
% Synthetic fields

% Filter solid fraction field with Gaussian filter and check histogram 

% Author: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear 
close all
clc

%% Create field and plot

% Resolution
res = 0.05;

% Coordinates fortest function
x = 0:res:50;
y = 0:res:100;

% Size of domain
Nx = length(x);
Ny = length(y);

% Initialize 
f = zeros(Ny,Nx);

%Amplitude
A1 = 1;
A2 = 1;

% Test function
for j = 1:Nx
    for k = 1:Ny
        f(k,j) = (sin(2*pi*x(j)/100).*sin(2*pi*y(k)/100) + A1*sin(2*pi*x(j)/10).*sin(2*pi*y(k)/10) + A2*sin(2*pi*x(j)/1).*sin(2*pi*y(k)/1));
    end
end

% Scale test function between a and b to mimic  particle phase 
% max[a,b]: [0,1]
a = 0;
b = 1;
minf = min(min(f));
maxf = max(max(f));
delf = maxf-minf;
for j = 1:length(x)
    for k = 1:length(y)
        f(k,j) = (b-a)*(f(k,j)-minf)/delf + a;
    end
end

%% Define filter function

% Filter width
filt = 20;
N = filt*2; %filt/2 + 2*filt;
r = -1*N:res:N;
s = -1*N:res:N;

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

%% Filter solid fraction

% Filter solid fraction field with convolution 
Xi_f = imfilter(f,G,'conv','circular');

%% Histogram

% Unfiltered 
figure(1)
h1 = histogram(f,100);
ycoords1 = h1.Values;
edgecoords = h1.BinEdges;
xcoords1 = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Filtered solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)
%axis([0 1 0 11000])
ax=gca; 
ax.YAxis.Exponent = 3;

% Filtered
figure(2)
h2 = histogram(Xi_f,100);
ycoords2 = h2.Values;
edgecoords = h2.BinEdges;
xcoords2 = (edgecoords(1:end-1)+edgecoords(2:end))./2;
hold on
xlabel('Filtered solid fraction')
ylabel('Number of cells')
set(gca,'FontSize',14)
%axis([0 1 0 11000])
ax=gca; 
ax.YAxis.Exponent = 3;