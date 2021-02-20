%% Figure 4 Panel F
% Computes and displays the correlations between pairs of synapses at
% different distances
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;
addpath(genpath('../tools'));
close all
%%
N = 30; L = 150;
close all;
if exist('../data/Figure4PanelF.mat')
    load('../data/Figure4PanelF.mat')
else
    % load all the mouse simulations
    fList = rdir('../sims/Fig4/NO*.mat');
    if isempty( fList )
        fprintf('run batch_process_mouse.m first\n')
    end

    accTOTCORRMATS = zeros(length(fList) , N*N); 
    accDMATS = zeros(length(fList) , N*N);

    for xx = 1:length(fList)
        xx
        dat = load(fList(xx).name , '-regexp' , '(fxs)|(pos)|(thetas)|(Sexcerpt)');
        S = dat.Sexcerpt(: , end-10000:end);
        % compute correlation
        totCor = corrcoef(smoothdata(S' , 'movmean' ,  60));
        accTOTCORRMATS(xx , :) = totCor(:);
        % compute distances
        sPos = mod((dat.pos + L/2),L);
        dMat = min(pdist2(dat.pos,dat.pos) , pdist2(sPos,sPos)) + diag(nan(N,1)); 
        accDMATS(xx , :) = dMat(:);
    end
    %%
    accDMATS = accDMATS';
    accTOTCORRMATS = accTOTCORRMATS';
        save('../data/Figure4PanelF.mat')
end
%%
%load experimental data
winnubstFig1e = csvread('../data/winnubstFig1E.csv'); winnubstFig1e = sortrows(winnubstFig1e , 1);

%%
% plotting
clear h
figure; 

h = gramm('x' , [accDMATS(:) ; winnubstFig1e(: , 1)] , 'y' , [100*accTOTCORRMATS(:) ; winnubstFig1e(: , 2)],...
    'color' , categorical([repmat("sim" , length(accDMATS(:)) , 1) ; repmat("data" , length(winnubstFig1e(:,1)) , 1)]));
h.stat_summary('geom' , 'point' , 'type' , 'std' , 'bin_in' , 40)
h.stat_summary('geom' , 'line' , 'type' , 'std' , 'bin_in' , 40)
h.set_color_options('map' , 'brewer_dark' )
h.set_order_options('color' , categorical(["sim" , "data"]))
h.set_names('x' , 'Pairwise distance' , 'y' , 'Correlation' , 'color' , 'data type')
h.axe_property('XLim' , [0 , 30] , 'YLim' , [0. , 40] , 'YTick' , [0.1 , 0.2 , 0.3 , 0.4]*100 , 'XTick' , [0:5:40])

h.draw
