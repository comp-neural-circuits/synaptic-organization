%% Figure 5 Panel E insets
% Generates the illustration of the averaged circular dispersion and
% receptive field offset across many simulations
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% October 2020;


addpath(genpath('../tools'));
close all
%%
close all;
if exist('../data/Figure5PanelDEinset.mat')
    load('../data/Figure5PanelDEinset.mat')
else
    % get all the morphological simulations with small receptive field
    % center spread
    fList = rdir('/gpfs/gjor/personal/kirchnerj/SIM_OUT/REDO_FERRET/Fig6/TREE*');%../sims/Fig5/BAP*.mat');
    if isempty( fList )
        fprintf('run batch_process_ferret_morphological.m first\n')
    end
    N = 504;

    CIRCDISP = zeros(length(fList) , N); 
    MUDIST = zeros(length(fList) , N); 
    POS = zeros(length(fList) , N); 
    somCONST = zeros(length(fList) , N);
    for xx = 1:length(fList)
        cFile = fList(xx).name
        dat = load(cFile , '-regexp' , '(compSomDist)|(pos)|(subpos)|(thetas)|(MUs)|(somConst)|(muVar)');
        thetas =  2*mod(dat.thetas,pi) - pi;
        % compute circular dispersion
        circ_dist_to_mean = 90*(abs(circ_dist2(thetas , circ_mean(thetas))))/pi ;
        CIRCDISP(xx , :) = circ_dist_to_mean;
        % compute receptive field center distance
        MUDIST(xx , :) = 62.5*sqrt(sum(dat.MUs.^2 , 2))/pi;
        % accumulate synapse positions and attenuation constants
        POS(xx , :) = dat.pos;
        somCONST(xx , :) = dat.somConst;
    end
end
%%
% compute averages on different parts of the tree
tab = table(); tab.POS = POS(:); tab.MUDIST = MUDIST(:); tab.CIRCDISP = CIRCDISP(:); 
tab.somCONST = somCONST(:);
outTab = varfun(@nanmean , tab , 'GroupingVariables' , {'somCONST','POS'} , 'InputVariables' , {'MUDIST'  , 'CIRCDISP'});
%%
% plot the averages on the tree
t = load_tree('../../master/data/L23.swc');
resConst = 10; 
tr = resample_tree(t , resConst);
%%
uSOMCON = unique(somCONST); 
% get positions
pIDS = [1 , cell2mat(arrayfun(@(x) find(full(tr.dA(x ,:)) == 1 , 1 , 'first') , 1:size(tr.dA , 1) , 'UniformOutput' , 0))];
dTree = [(tr.X - tr.X(pIDS)) , (tr.Y - tr.Y(pIDS)) , (tr.Z - tr.Z(pIDS))]; 
dTree = normalize(dTree , 2 , 'norm');
[~ , atr] = sub_tree(tr , find(tr.R == 3 | tr.R == 4 | tr.R == 1));
%%
close all;
figure; 
for xx2 = 1:length(uSOMCON)
    sTab = outTab(outTab.somCONST == uSOMCON(xx2) , :);
    aPos = [tr.X(sTab.POS) , tr.Y(sTab.POS) , tr.Z(sTab.POS)];
    tMAP = flipud(cbrewer('seq' , 'Greens' , length(tr.X) + 100)); tMAP = tMAP(1:end-100 , :);
    s1 = subplot(2,length(uSOMCON),xx2);
    plot_tree(atr, rgb('grey') , [] , [] , [] , '-2l -thick'); 
    hold on; 
    s = scatter(aPos(: , 1) , aPos(: , 2) , 15 , sTab.nanmean_CIRCDISP , 'filled');
    s.MarkerEdgeColor = 'none';%
    s.LineWidth = 0.1;
    colormap(s1 , cbrewer('seq' , 'Greens' , 100));
    axis off
    shading interp
    caxis([20 , 45])

    s2 = subplot(2,length(uSOMCON),xx2+length(uSOMCON));
    plot_tree(atr, rgb('grey') , [] , [] , [] , '-2l -thick'); 
    hold on; 
    s = scatter(aPos(: , 1) , aPos(: , 2) , 15 , sTab.nanmean_MUDIST , 'filled');
    caxis([ 0 , 12.5])

    s.MarkerEdgeColor = 'none';%
    s.LineWidth = 0.1;
    colormap(s2 , cbrewer('seq' , 'Purples' , 100))
    axis off
    shading interp
end

