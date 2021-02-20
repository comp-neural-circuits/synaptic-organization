%% Figure 4 Panel D
% Generates the illustration the lack of orientation clustering for the
% mouse simulations
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;
addpath(genpath('../tools'));
close all
clear all
% load the simulation used in the paper
load('../data/NOBAP_MOUSE_BIOLOGICAL_737530.85233633569441735744.mat')

% plotting
figure; 
subplot(1,2,1)
[~ , sID] = sort(pos);
plot([0 , 0] , [0 , L] , 'LineWidth',5 , 'Color' , rgb('grey'));
hold on;
xJitt = (rand(N,1)-0.5)/30;
s = scatter(xJitt , pos , 50 , storeConfig.thetas(:,1) , 'filled');
s.MarkerEdgeColor = rgb('black');

xlim([-0.5 , 0.5])
ylim([0 , L])
caxis([0,2*pi])
colorbar;
colormap(getWilsonMap(100));
ylabel('Distance \mu m')
title('Synapses on Tree')
box off
subplot(1,2,2)
[~ , sID] = sort(pos);
plot([0 , 0] , [0 , L] , 'LineWidth',5 , 'Color' , rgb('grey'));
hold on;
xJitt = (rand(N,1)-0.5)/30;
s = scatter(xJitt , pos , 50 , thetas , 'filled');
s.MarkerEdgeColor = rgb('black');
xlim([-0.5 , 0.5])
ylim([0 , L])
caxis([0,2*pi])
colorbar;
colormap(getWilsonMap(100));
ylabel('Distance \mu m')
title('Synapses on Tree')
box off