%% -----------------------------------------------------------------------
% Global mixing index
% Scale-dependent mixing
% 2D CFD particle segregation

% Plotting global mixing index with fitlering and in time
% Values are obtained from scale dependent mixing global index scritps

% Author: Barlev R. Nagawkar, Iowa State University
% ------------------------------------------------------------------------

clear 
close all
clc

%% Global index based on filtered solid fraciton

time = 0:5:30;

GI_30 = [0.7914 0.7916 0.7956 0.8038 0.8176];
GI_25 = [0.7589 0.7591 0.7620 0.7699 0.7903]; 
GI_20 = [0.8044 0.8045 0.8080 0.8129 0.8208];
GI_15 = [0.8012 0.8013 0.8053 0.8122 0.8244];
GI_10 = [0.7771 0.7771 0.7798 0.7834 0.7942];
GI_05 = [0.7869 0.7870 0.7908 0.7955 0.8068]; 
GI_00 = [1 1 0.9777 0.9695 0.9471];

GI = [GI_00; GI_05; GI_10; GI_15; GI_20; GI_25; GI_30];

figure(1)
plot(time,GI(:,1), '--*', 'Color','#0072BD',  'LineWidth',3)
hold on
plot(time,GI(:,2),'-*', 'Color', '#D95319', 'LineWidth',2)
plot(time,GI(:,3),'-*', 'Color', '#EDB120', 'LineWidth',2)
plot(time,GI(:,4),'-*', 'Color', '#7E2F8E', 'LineWidth',2)
plot(time,GI(:,5),'-*', 'Color', '#77AC30', 'LineWidth',2)
xlabel('Time, s')
ylabel('Global index')
set(gca,'FontSize',18)
hlegend=legend('Unfiltered','\Delta = 0.001','\Delta = 0.004','\Delta = 0.01','\Delta = 0.03');

%% Global index based on accessed region mixing

filt = [0.006 0.012 0.02 0.04];
gm = [0.7974 0.8048 0.8131 0.8171]; % Solid fraction
amgm = [0.5785 0.5853 0.5920 0.5893]; % Accessed mixing 

figure(2)
plot(filt,gm, '-*k', 'Linewidth',2)
hold on
plot(filt,amgm, '--*k', 'Linewidth',2)
xlabel('Filter width, m')
ylabel('Global index')
legend('Solid fraciton', 'Accessed mixing index, [0.3. 0.7]')
axis([0 0.04 0.55 1])
set(gca,'FontSize',18)