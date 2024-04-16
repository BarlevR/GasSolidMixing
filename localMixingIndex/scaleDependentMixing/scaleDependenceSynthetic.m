%% -----------------------------------------------------------------------
% Local mixing index
% Solid index for solid-solid mixing
% Scale-dependent mixing
% Synthetic fields

% Filter field with Gaussian filter and then compute local mixing index
% Both solid mixing and accessed region indices used in this script

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
a = 0.45;
b = 0.95;
minf = min(min(f));
maxf = max(max(f));
delf = maxf-minf;
for j = 1:length(x)
    for k = 1:length(y)
        f(k,j) = (b-a)*(f(k,j)-minf)/delf + a;
    end
end

% Plot scaled test function
figure(1)
surf(x,y,f,'EdgeColor', 'none')
hold on
xlabel('X')
ylabel('Y')
title('Test function')
colorbar
colormap(jet(1024))
clim([0 1])
view(0,90)
daspect([1 1 1])
axis ([0 50 0 100])
set(gca,'FontSize',14)

%% Define filter function

% Filter width
filt = 1;
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

%Plot filter fucntion
figure(2)
surf(r,s,G,'EdgeColor', 'none')

%% Filter solid fraction and plot

% Filter solid fraction field with convolution 
Xi_f = imfilter(f,G,'conv','circular');

%Plot filtered solid fraction field
figure(3)
surf(x,y,Xi_f,'EdgeColor', 'none')
view(0,90)
axis([0 50 0 100])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Filtered solid fraction function')
set(gca,'FontSize',14)

%% Mixing index

Xi_f_mean = mean(mean(Xi_f)); 
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

% Plot mixing index
figure(4)
surf(x,y,M,'EdgeColor', 'none')
hold on
xlabel('X')
ylabel('Y')
title('Mixing index')
colorbar
clim ([0 1])
colormap(jet(1024))
view(0,90)
daspect([1 1 1])
axis ([0 50 0 100])
set(gca,'FontSize',14)

%% Accessed region mixing index

A = zeros(Ny,Nx); % initialize

% Select accessed region range [Xi_f_min, Xi_f_max]
Xi_f_min = 0.45;
Xi_f_max = 0.95;
 
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

figure(5)
surf(x,y,A,'EdgeColor', 'none')
view(0,90)
axis([0 50 0 100])
daspect([1 1 1])
colorbar
colormap(jet(1024))
clim([0 1])
xlabel('X, m')
ylabel('Y, m')
title('Accessed region mixing index')
set(gca,'FontSize',14)