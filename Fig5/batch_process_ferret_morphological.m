%% batch process ferret morphological
% runs simulations required for the morphological dendrite with receptive field
% spread as for the ferret
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

KK = 50; % number of simulations
addpath(genpath('../tools/'))
rng('shuffle');

%%
load('../data/wavmat.mat');

%%
ferretSpreadInDegree = 5.3; fullVisualFieldDegree = 100;
muSTD = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree); % receptive field center spread converted to custom coordinate system

RFsize = 13.4;
sCs = [25 , 75 , 125];
for xx = 1:KK
    for xx2 = 1:length(sCs)
        dynamic_clustering_morpho(SWM , muSTD, RFsize , sCs(xx2) , 'Fig5/BAP_FERRET_MORPHO');
    end
end
