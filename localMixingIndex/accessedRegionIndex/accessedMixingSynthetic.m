%% -----------------------------------------------------------------------
% Local mixing index
% Accessed region solid mixing index for solid-solid mixing
% Synthetic field

% Authors: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear
close all
clc

%% Create field and plot

% Coordinates fortest function
x = 0:0.05:50;
y = 0:0.05:100;

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
a = 0.1;
b = 0.75;
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

%% Calculate the richness Index


f_mean = mean(mean(f));
% Pick accessed region
f_min = a;
f_max = b;

A = zeros(Ny,Nx); % initialize
for j = 1:Nx
    for k = 1:Ny
        if f(k,j) < f_min
           A(k,j) = 0; 
        elseif (f(k,j)>= f_min && f(k,j) < f_mean)
            A(k,j) =  (f(k,j) - f_min)/(2*(f_mean - f_min));
        elseif (f(k,j) >= f_mean && f(k,j) < f_max)
            A(k,j) = 1 - (f_max - f(k,j))/(2*(f_max - f_mean));
        else
            A(k,j) = 1;
        end
    end
end

% Plot mixing index
figure(2)
surf(x,y,A,'EdgeColor', 'none')
hold on
xlabel('X')
ylabel('Y')
title('Mixing index')
colorbar
colormap(jet(1024))
clim([0 1])
view(0,90)
daspect([1 1 1])
axis ([0 50 0 100])
set(gca,'FontSize',14)

