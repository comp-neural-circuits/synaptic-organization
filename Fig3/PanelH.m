%% Figure 3 Panel H
% Computes and displays the decay of the correlation vs. distance curve for
% simulations with different calcium diffusion constants
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
N = 30; L = 150; normSTD = 6;
close all;
if 0;% exist('../data/Figure4PanelH.mat')
    load('../data/Figure3PanelH.mat')
else
    % load all retinal wave simulations that have different values for the
    % calcium constant
    fList = rdir('../sims/Fig3VaryCal/NO*.mat');
    if isempty( fList )
        fprintf('run batch_process_varyCal.m first\n')
    end
    accNORMSTD = zeros(length(fList) , N*N); 
    accTOTCORRMATS = zeros(length(fList) , N*N); 
    accDMATS = zeros(length(fList) , N*N);
    %%
    for xx = 1:length(fList)
        xx
        dat = load(fList(xx).name , '-regexp' , '(pos)|(thetas)|(Sexcerpt)|(normSTD)');
        % compute correlations
        S = dat.Sexcerpt;
        totCor = corrcoef(smoothdata(S' , 'movmean' ,  60));   
        accTOTCORRMATS(xx , :) = totCor(:);
        % store all the diffusion constants
        accNORMSTD(xx , :) = dat.normSTD;
        % compute distances
        sPos = mod((dat.pos + L/2),L);
        dMat = min(pdist2(dat.pos,dat.pos) , pdist2(sPos,sPos)) + diag(nan(N,1)); 
        accDMATS(xx , :) = dMat(:);
    end
    %%
    accDMATS = accDMATS'; accDMATS = accDMATS(:);
    accTOTCORRMATS = accTOTCORRMATS'; accTOTCORRMATS = accTOTCORRMATS(:);
    accNORMSTD = accNORMSTD'; accNORMSTD = accNORMSTD(:);
    save('../data/Figure3PanelH.mat')
end

%%
% plotting
uNORMSTD = unique(accNORMSTD(:));
cMAP = cbrewer('seq' , 'Greys' , 30); cMAP = cMAP(11:end , :);
% fit a Gaussian
gaussEqn = 'a*exp(-(x.^2/(2*l)^2))';
dX = linspace(0 , 75);
figure; 
hold on;
ls_sim = []; as_sim = []; 
for ii = 1:length(uNORMSTD)
    useID = ~isnan(accDMATS) & ~isnan(accTOTCORRMATS);
    distances = accDMATS(accNORMSTD == uNORMSTD(ii) & useID);
    correlations = accTOTCORRMATS(accNORMSTD == uNORMSTD(ii) & useID);
    % subtract the baseline
    correlationsMBL = correlations - nanmean(correlations(distances > prctile(distances , 50) & distances < prctile(distances , 75)));
    f = fit( distances , correlationsMBL,  gaussEqn , 'Lower' , [0 , 0 ] ,'Upper' , [inf, inf] , 'StartPoint' , [0.11 , uNORMSTD(ii) ] );
    ls_sim = [ls_sim , f.l];     as_sim = [as_sim , f.a];    
end
% scatter the decay constants
plot(0:15 , 0:15 , 'Color' , rgb('black'))
for ii = 1:length(uNORMSTD)
    scatter(uNORMSTD(ii) , ls_sim(ii) , 'MarkerFaceColor' , cMAP(ii , :) , 'MarkerEdgeColor' , rgb('black') )
end
axis square; xlabel('sigma_c'); ylabel('cluster sigma'); 
