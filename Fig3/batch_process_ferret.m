%% batch process ferret
% runs simulations required for the linear dendrite with receptive field
% spread as for the ferret
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;


KK = 250; % number of simulations
addpath('../tools/')
rng('shuffle');

%%
load('../data/wavmat.mat'); 

%%
ferretSpreadInDegree = 5.3; fullVisualFieldDegree = 100;
muSTD = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree); % receptive field center spread converted to custom coordinate system

RFsize = 13.4;

for xx = 1:KK
    dynamic_clustering(SWM , muSTD , RFsize , 'Fig3/NOBAP_FERRET_LINEAR');
end
%%