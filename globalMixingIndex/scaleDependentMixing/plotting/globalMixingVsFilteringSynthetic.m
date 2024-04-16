%% -----------------------------------------------------------------------
% Global mixing index
% Scale-dependent mixing
% Synthetic fields

% Plotting global mixing index vs filter width
% Values are obtained from scale dependent mixing global index scritps

% Author: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear 
close all
clc

%% Global index based on solid fraction Lacey vs Rose

filt = [0.01 0.1 0.5 1 2 3 5 8 12 20];
GI = [0.7069 0.7100 0.7503 0.7644 0.7751 0.7895 0.8161 0.8328 0.8378 0.8482];
LI = [0.9141 0.9159 0.9377 0.9445 0.9494 0.9557 0.9663 0.9720 0.9737 0.9770];
x1 = [1 1];
x10 = [10 10];
x100 = [100 100];
x0 = [0 0];
y = [0 1];

figure(1)
plot(filt,LI, '-*k', 'Linewidth',2)
hold on
plot(filt,GI, '--*k', 'Linewidth',2)
plot(x0,y,'r', 'Linewidth',2)
plot(x1,y,'r', 'Linewidth',2)
plot(x10,y,'r', 'Linewidth',2)
plot(x100,y,'r', 'Linewidth',2)
xlabel('Filter width')
ylabel('Global index')
axis([0 21 0.7 1])
legend('Lacey', 'Rose')
set(gca,'FontSize',18)
%% Global index based on local mixing Rose

filt = [0.01 0.1 2 5 20];
gm = [0.8402 0.8419 0.8773 0.8988 0.9172];
mgm = [0.8206 0.8225 0.8624 0.8864 0.9076];
amgm = [0.7071 0.7102 0.7751 0.8165 0.8482];

figure(2)
plot(filt,gm, '-*k', 'Linewidth',2)
hold on
plot(filt,mgm, '--*k', 'Linewidth',2)
plot(filt,amgm, '-ok', 'Linewidth',2)
plot(x0,y,'r', 'Linewidth',2)
plot(x1,y,'r', 'Linewidth',2)
plot(x10,y,'r', 'Linewidth',2)
plot(x100,y,'r', 'Linewidth',2)
xlabel('Filter width')
ylabel('Global index')
legend('Solid fraction', 'Solid mixing index, [0,1]', 'Accessed mixing index, [0.45. 0.95]')
set(gca,'FontSize',18)
axis([0 21 0.7 1])