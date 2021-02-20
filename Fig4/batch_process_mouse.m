%% batch process mouse
% runs simulations required for the linear dendrite with receptive field
% spread as for the mouse
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;


KK = 250; % number of simulations
addpath('../tools/')

%%
load('../data/wavmat.mat'); WM = SWM; clear SWM

%%
mouseSpreadInDegree = 26; fullVisualFieldDegree = 100;
muSTD = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree); % receptive field center spread converted to custom coordinate system

RFsize = 20;

for xx = 1:KK
    dynamic_clustering(WM , muSTD , RFsize , 'Fig5/NOBAP_MOUSE_LINEAR');
end
