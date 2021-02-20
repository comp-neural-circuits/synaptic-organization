%% Figure 4 Panel E and G
% Computes and plots the delta orientation and the receptive field overlap between
% pairs of synapses at increasing distances
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;

addpath(genpath('../tools'));
close all
%%
N = 30; L = 150;
close all;
if 0;% exist('../data/Figure5PanelCE.mat')
    load('../data/Figure4PanelEG.mat')
else
    % load all the mouse simulations
    fList = rdir('../sims/Fig4/NO*.mat');
    if isempty( fList )
        fprintf('run batch_process_mouse.m first\n')
    end

    accSPCORRMATS = zeros(length(fList) , N*N); 
    accDMATS = zeros(length(fList) , N*N);
    accTMATS = zeros(length(fList) , N*N);

    for xx = 1:length(fList)
        xx
        dat = load(fList(xx).name , '-regexp' , '(fxs)|(pos)|(thetas)');
        % compute receptive field overlap
        tmpFX = reshape(dat.fxs  , [73*73 , N]); 
        spCORRMAT = corrcoef(tmpFX ) + diag(nan(N,1));
        accSPCORRMATS(xx , :) = spCORRMAT(:);
        % computes distances
        sPos = mod((dat.pos + L/2),L);
        dMat = min(pdist2(dat.pos,dat.pos) , pdist2(sPos,sPos)) + diag(nan(N,1)); 
        accDMATS(xx , :) = dMat(:);
        % compute the delta orientation modulo 180 degree
        thetas = 2*mod(dat.thetas , pi);
        tMat = abs(circ_dist2(thetas - pi , thetas - pi))/2 + diag(nan(N,1));
        accTMATS(xx , :) = tMat(:);
    end
        save('../data/Figure4PanelEG.mat')
end


%%
% load experimental data
useID = ~isnan(accDMATS(:)) & ~isnan(accTMATS(:)) & ~isnan(accSPCORRMATS(:)) &...
    ~isinf(accDMATS(:)) & ~isinf(accTMATS(:)) & ~isinf(accSPCORRMATS(:));
clear g h f
iacaTab1g = csvread('../data/IacarusoFig1g.csv'); iacaTab1g = sortrows(iacaTab1g , 1);
iacaTab1h = csvread('../data/IacarusoFig1h.csv'); iacaTab1h = sortrows(iacaTab1h , 1);


%%
% plotting
figure; 

g = gramm('x' , [accDMATS(useID) ; iacaTab1g(:,1)] , 'y' , [accSPCORRMATS(useID) ; iacaTab1g(:,2)] ,...
    'color' , categorical([repmat("sim" , length(accDMATS(useID)) , 1) ; repmat("data" , length(iacaTab1g(:,1)) , 1) ; ]));
g.stat_summary('geom' , 'line' , 'type' , 'std' , 'bin_in' , 20)
g.stat_summary('geom' , 'point' , 'type' , 'std' , 'bin_in' , 20)

g.set_names('x' , 'Pairwise distance' , 'y' , 'Spatial correlation')
g.axe_property('XLim' , [0 , 30] , 'YLim' , [0 , 0.3])
g.set_order_options('color' , categorical(["sim" , "data"]))

g.draw

figure; 

h = gramm('x' , [accDMATS(useID) ; iacaTab1h(:,1)] , 'y' , [180*accTMATS(useID)/pi ; iacaTab1h(:,2) ],...
    'color' , categorical([repmat("sim" , length(accDMATS(useID)) , 1) ; repmat("data" , length(iacaTab1h(:,1)) , 1) ; ]));
h.stat_summary('geom' , 'line' , 'type' , 'std' , 'bin_in' , 20)
h.stat_summary('geom' , 'point' , 'type' , 'std' , 'bin_in' , 20)
h.set_names('x' , 'Pairwise distance' , 'y' , 'Delta orientation')
h.axe_property('XLim' , [0 , 30] , 'YLim' , [0 , 50] , 'YTicks' , [0 , 20 , 40])
h.set_order_options('color' , categorical(["sim" , "data"]))

h.draw

%%
autoArrangeFigures
