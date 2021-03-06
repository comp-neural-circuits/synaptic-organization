%% Figure 6 Panel A
% Generates the illustration of circular dispersion and receptive field
% center offset for one simulation with a morphologically reconstructed
% dendritic tree and a large  receptive field center offset
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;

addpath(genpath('../tools'));
close all
clear all
%%
% load the simulation use for the figure
load('../../master/data/BAP_MOUSE_TREE_BIOLOGICAL_737645.26905233610887080431.mat')
% plotting
pIDS = [1 , cell2mat(arrayfun(@(x) find(full(tr.dA(x ,:)) == 1 , 1 , 'first') , 1:size(tr.dA , 1) , 'UniformOutput' , 0))];
dTree = [(tr.X - tr.X(pIDS)) , (tr.Y - tr.Y(pIDS)) , (tr.Z - tr.Z(pIDS))]; 
dTree = normalize(dTree , 2 , 'norm');
aPos = [tr.X(pos) , tr.Y(pos) , tr.Z(pos)] + (subpos).*dTree(pos , :);
[~ , apt] = sub_tree(tr , find(tr.R == 3 | tr.R == 4 | tr.R == 1));
atr = apt;
tMAP = flipud(cbrewer('seq' , 'Greens' , length(tr.X) + 100)); tMAP = tMAP(1:end-100 , :);
figure; 
s1 = subplot(1,2,1);
plot_tree(atr, rgb('grey') , [] , [] , [] , '-2l -thick'); 
hold on; 
s = scatter(aPos(: , 1) , aPos(: , 2) , 25 , thetas , 'filled');
s.MarkerEdgeColor = rgb('black');
s.LineWidth = 0.25;
colormap(s1 , getWilsonMap(100));
axis off
shading interp

s2 = subplot(1,2,2);
PRGNmap = cbrewer('div' , 'PRGn' , 200); PRGNmap = PRGNmap(50:150 , :);
plot_tree(atr, rgb('grey') , [] , [] , [] , '-2l -thick'); 
hold on; 
s = scatter(aPos(: , 1) , aPos(: , 2) , 40 , 62.5*sqrt(sum(MUs.^2 , 2))/pi , 'filled');
caxis([ 10 , 30])
s.MarkerEdgeColor = rgb('black');
s.LineWidth = 0.25;
colormap(s2 , PRGNmap);%cbrewer('div' , 'PRGn' , 100));%'seq' , 'Purples' , 100))
colorbar
axis off
shading interp
