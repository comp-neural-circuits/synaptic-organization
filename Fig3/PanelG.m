%% Figure 3 Panel G
% Computes and displays the correlation between the activity of pairs of synapses at
% different distance
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
% close all
%%
N = 30; L = 150;
% close all;
if 0;%exist('../data/Figure4PanelG.mat')
    load('../data/Figure3PanelG.mat')
else
    % load all retinal wave simulation
    fListFerr = rdir('../sims/Fig3/NO*.mat');
    if isempty( fListFerr )
        fprintf('run batch_process_ferret.m first\n')
    end
    % load all white noise simulation
    fListWN = rdir('../sims/Fig3WN/NO*.mat');
    if isempty( fListFerr )
        fprintf('run batch_process_wn.m first\n')
    end
    fList = [fListFerr ; fListWN];
    accSIMTYPE = [zeros(length(fListFerr) , N*N) ; ones(length(fListWN) , N*N)]';
    accTOTCORRMATS = zeros(length(fList) , N*N); 
    accDMATS = zeros(length(fList) , N*N);
    %%
    % iterate over all simulations
    for xx = 1:length(fList)
        xx
        dat = load(fList(xx).name , '-regexp' , '(pos)|(thetas)|(Sexcerpt)');
        % compute correlation matrix
        S = dat.Sexcerpt;
        totCor = corrcoef(smoothdata(S' , 'movmean' ,  60)) + diag(nan(N,1)); %   corrcoef(R'); %   
        accTOTCORRMATS(xx , :) = totCor(:);
        % compute distances
        sPos = mod((dat.pos + L/2),L);
        dMat = min(pdist2(dat.pos,dat.pos) , pdist2(sPos,sPos)) + diag(nan(N,1)); 
        accDMATS(xx , :) = dMat(:);

    end
    %%
    accDMATS = accDMATS';
    accTOTCORRMATS = accTOTCORRMATS';
    save('../data/Figure3PanelG.mat')
end
%%
% load experimental data
schollFig5C = csvread('../data/SchollFig5c.csv'); schollFig5C = sortrows(schollFig5C , 1);
useID = ~isnan(accDMATS(:))  & ~isnan(accTOTCORRMATS(:)) &...
    ~isinf(accDMATS(:))  & ~isinf(accTOTCORRMATS(:)) & accDMATS(:) < 20;
%%
% plotting
figure; 

h = gramm('x' , [accDMATS(useID) ; schollFig5C(: , 1)] , 'y' , [accTOTCORRMATS(useID) ; schollFig5C(: , 2)],...
    'color' , [categorical(accSIMTYPE(useID)) ; categorical(repmat("data" , length(schollFig5C(:,1)) , 1))]); % 'color' , categorical([repmat("sim" , length(accDMATS(:)) , 1) ; repmat("data" , length(schollFig5C(:,1)) , 1) ]),...
h.stat_summary('geom' , 'point' , 'type' , 'std' , 'bin_in' , 10)
h.stat_summary('geom' , 'line' , 'type' , 'std' , 'bin_in' , 10)
h.set_color_options('map' , 'brewer_dark' )
h.set_names('x' , 'Pairwise distance' , 'y' , 'Correlation' , 'color' , 'data type')
h.axe_property('XLim' , [0 , 20] , 'YLim' , [-0.025 , 0.3] , 'YTick' , [-0.1 , 0 , 0.1 , 0.2 , 0.3 , 0.4] , 'XTick' , [0:5:40] , 'PlotBoxAspectRatio' , [1,1,1] )

h.draw
