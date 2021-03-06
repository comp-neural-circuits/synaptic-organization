%% Figure 5 Panel E
% Generates the illustration of the clustered morphological neuron with
% homogeneous orientation preference at proximal synapses
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;


addpath(genpath('../tools'));
close all
clear all
% load('../data/BAP_FERRET_TREE_BIOLOGICAL_737572.00262520276010036469.mat')
% load('../data/TREE_FERRET_LINEAR_737832.10825115267653018236.mat')
load('/gpfs/gjor/personal/kirchnerj/SIM_OUT/REDO_FERRET/Fig6/TREE_FERRET_LINEAR_737831.27746278641279786825.mat')
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
plot_tree(atr, rgb('grey') , [] , [] , [] , '-2l -thick'); 
hold on; 
s = scatter(aPos(: , 1) , aPos(: , 2) , 25 , 62.5*sqrt(sum(MUs.^2 , 2))/pi , 'filled');
caxis([ 0 , 12.5])
s.MarkerEdgeColor = rgb('black');
s.LineWidth = 0.25;
colormap(s2 , cbrewer('seq' , 'Purples' , 100))
axis off
shading interp
