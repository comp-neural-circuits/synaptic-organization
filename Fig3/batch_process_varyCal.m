%% batch process ferret
% runs simulations required for the linear dendrite with receptive field
% spread as for the ferret
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;


KK = 50; % number of simulations
addpath('../tools/')
rng('shuffle');

%%
load('../data/wavmat.mat'); WM = SWM; clear SWM

%%
ferretSpreadInDegree = 5.3; fullVisualFieldDegree = 100;
muSTD = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree); % receptive field center spread converted to custom coordinate system

RFsize = 13.4;
cspread = 1:15;
for xx = 1:KK
    for xx2 = 1:length(cspread)
        dynamic_clustering_clu_size(WM , muSTD , RFsize , cspread(xx2) , 'Fig3VaryCal/NOBAP_FERRET_LINEAR');
    end
end
