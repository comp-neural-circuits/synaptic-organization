%% Figure 3 Panel F
% Computes and displays the difference in orientation preference for pairs
% of synapses at different distances
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
N = 30; L = 150; normSTD = 6;
close all;
if 0; %exist('../data/Figure4PanelF.mat')
    load('../data/Figure3PanelF.mat')
else
    % load all retinal wave simulations
    fListFerr = rdir('../sims/Fig3/NO*.mat');
    if isempty( fListFerr )
        fprintf('run batch_process_ferret.m first\n')
    end
    % load all white noise simulations
    fListWN = rdir('../sims/Fig3WN/NO*.mat');
    if isempty( fListFerr )
        fprintf('run batch_process_wn.m first\n')
    end
    fList = [fListFerr ; fListWN];
    accSIMTYPE = [zeros(length(fListFerr) , N*N) ; ones(length(fListWN) , N*N)];
    accSPCORRMATS = zeros(length(fList) , N*N); 
    accDMATS = zeros(length(fList) , N*N);
    accTMATS = zeros(length(fList) , N*N);
    for xx = 1:length(fList)
        xx
        dat = load(fList(xx).name , '-regexp' , '(pos)|(thetas)');
        % compute distances
        sPos = mod((dat.pos + L/2),L);
        dMat = min(pdist2(dat.pos,dat.pos) , pdist2(sPos,sPos)) + diag(nan(N,1)); 
        accDMATS(xx , :) = dMat(:);
        % compute the differences in orientation preference modulo 180
        % degree
        thetas = 2*mod(dat.thetas , pi);
        tMat = abs(circ_dist2(thetas - pi , thetas - pi))/2 + diag(nan(N,1));
        accTMATS(xx , :) = tMat(:);
    end
    save('../data/Figure3PanelF.mat')
end


%%
% load experimental data
useID = ~isnan(accDMATS(:)) & ~isnan(accTMATS(:)) & ...
    ~isinf(accDMATS(:)) & ~isinf(accTMATS(:)) & accDMATS(:) < 30;
WilsonSupFig7a = csvread('../data/WilsonSupFig7a.csv'); WilsonSupFig7a = sortrows(WilsonSupFig7a , 1);
%%
% plotting
figure; 
h = gramm('x' , accDMATS(useID)  , 'y' , 180*accTMATS(useID)/pi , 'color' , accSIMTYPE(useID));
h.stat_summary('geom' , 'line' , 'type' , 'std' , 'bin_in' , 15);
h.set_names('x' , 'Pairwise distance' , 'y' , 'Delta orientation');
h.geom_hline('yintercept' , 45 );
h.axe_property('XLim' , [0 , 30] , 'YLim' , [20 , 50] , 'XTick' , 0:10:30 , 'YTick' , 20:10:60,'PlotBoxAspectRatio' , [1 , 1 , 1]);
h.draw;
%%
%%
figure;
g = gramm('x' , WilsonSupFig7a(:,1) , 'y' , cumsum(WilsonSupFig7a(:,2)));
g.geom_bar
g.set_color_options('map' , rgb('black'))
snapnow;
g.update('x' , 180*accTMATS(accDMATS(:) < 3)/pi , 'subset' , accSIMTYPE(accDMATS(:) < 3) == 0);
g.stat_bin('nbins' , 4 , 'normalization' , 'cdf' , 'geom' , 'line')
g.axe_property('XLim' , [0 , 90] , 'YLim' , [0 , 1] , ...
    'XTick' , [0 , 45 , 90] , 'YTick' , 0:0.25:1, ...
    'PlotBoxAspectRatio' , [1 , 1 , 1])
g.set_color_options('map' , rgb('red'))
g.set_names('x' , 'Delta orientation')
g.draw;