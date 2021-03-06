%% Figure 6 Panel E, F
% Computes and displays the circular dispersion in co-axial and orthogonal
% space
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;

% close all
addpath(genpath('../tools'));
close all
%%
if exist('../data/Figure6PanelEFG.mat')
    load('../data/Figure6PanelEFG.mat')
else
    % get all simulations with morphology and large receptive field center
    % spread
    fList = rdir('../sims/Fig6/BAP*.mat');
    if isempty( fList )
        fprintf('run batch_process_mouse.m first\n')
    end
    N = 504;

    MUDIST = zeros(length(fList) , N); 
    MUROT = zeros(length(fList) , N , 2); 
    MUANGLE = zeros(length(fList) , N); 
    THETAS = zeros(length(fList) , N); 

    POS = zeros(length(fList) , N); 
    somCONST = zeros(length(fList) , N);
    for xx = 1:length(fList)
        cFile = fList(xx).name
        dat = load(cFile , '-regexp' , '(compSomDist)|(pos)|(subpos)|(thetas)|(MUs)|(somConst)|(somMU)');
        % collect centered orientation preferences
        thetas =  dat.thetas - pi;%
        centeredThetas = circ_dist2(thetas , circ_mean(thetas));
        mTheta = circ_mean(thetas);
        % collect centered receptive field centers
        MUs = dat.MUs;
        MUDIST(xx , :) = 62.5*sqrt(sum(MUs.^2 , 2))/pi;
        MUANGLE(xx , :) = 180*angle(MUs(: , 1) + i*MUs(:,2))/pi;
        R = [cos(-mTheta) -sin(-mTheta) ;sin(-mTheta) cos(-mTheta)] ;
        MUROT(xx , : , :) =  (R*dat.MUs')';
        THETAS(xx , :) = centeredThetas;
        % collect positions and attenuation constants
        POS(xx , :) = dat.compSomDist(dat.pos) + dat.subpos;
        somCONST(xx , :) = dat.somConst;
    end
end
%%
% compute circular dispersion
uCONST = unique(somCONST(:));
modTHETAS  = mod(THETAS(:),pi);
modCIRCDISP = min(modTHETAS , abs(modTHETAS - pi));
MUROTX = MUROT(: , : , 1); MUROTY = MUROT(: , : , 2);
VANGLE = 180*angle(MUROTX(:) + MUROTY(:)*i)/pi;
%%
% plotting
f = figure;
k = gramm('x' , 62.5*MUROTX(:)/pi , 'y' , 62.5*MUROTY(:)/pi , 'color' , 180*THETAS(:)/pi , ... 
    'subset' , logical((somCONST(:) == uCONST(2)).*(1:N*length(fList) < 10*N)'));
k.geom_point;
k.set_point_options('base_size' , 2 );
k.axe_property('PlotBoxAspectRatio' , [1 , 1 , 1] , 'YLim' , [-50 , 50] , ...
    'XTick' , [-50 , 0 , 50] , 'XLim' , [-50 , 50] , 'YTick' , [-50 , 0 , 50]);
k.set_continuous_color('CLim' , [-180 , 180] , 'colormap' , 'custom' , 'customColormap' , getWilsonMap(256))
k.set_names('x' , 'Azimuth' , 'y' , 'Elevation')
k.draw;
k.results.geom_point_handle(1).MarkerEdgeColor = 'flat';
k.results.geom_point_handle(1).MarkerFaceColor = 'none';

%%
f = figure;
k = gramm('x' , 62.5*MUROTX(:)/pi , 'y' , 62.5*MUROTY(:)/pi , 'color' , 180*modCIRCDISP(:)/pi , ... 
    'subset' , logical((somCONST(:) == uCONST(2))));
k.facet_grid([] , (VANGLE > -45 & VANGLE < 45) | (VANGLE > 135 )  | VANGLE < -135)
k.geom_point;
k.set_point_options('base_size' , 2 );
k.axe_property('PlotBoxAspectRatio' , [1 , 1 , 1] , 'YLim' , [-50 , 50] , ...
    'XTick' , [-50 , 0 , 50] , 'XLim' , [-50 , 50] , 'YTick' , [-50 , 0 , 50]);
k.set_continuous_color('CLim' , [0 , 90] , 'colormap' , 'custom' , 'customColormap' , flipud(jet(256)))
k.set_names('x' , 'Azimuth' , 'y' , 'Elevation')
k.geom_abline('intercept' , 0 , 'slope' , 1)
k.geom_abline('intercept' , 0 , 'slope' , -1)
k.draw;
k.results.geom_point_handle(1).MarkerEdgeColor = 'flat';
k.results.geom_point_handle(1).MarkerFaceColor = 'none';
k.results.geom_point_handle(2).MarkerEdgeColor = 'flat';
k.results.geom_point_handle(2).MarkerFaceColor = 'none';
%%
% get experimental data
iacaTab3d = csvread('../data/IacarusoFig3d.csv');
selectID = logical(somCONST(:) == uCONST(2));
COAX = ((VANGLE > -45 & VANGLE < 45) | (VANGLE > 135 )  | VANGLE < -135);
NUMCOAX = histcounts(180*modCIRCDISP(selectID & COAX)/pi , [0 , 30 , 60 , 90]);
NUMORTH = histcounts(180*modCIRCDISP(selectID & ~COAX)/pi , [0 , 30 , 60 , 90]);
totNUM = sum(selectID);
%%

f = figure;
k2 = gramm('x' , [15 , 45 , 75 , 15 , 45 , 75] , 'y' , [NUMCOAX/totNUM , NUMORTH/totNUM]);
k2.facet_grid([] , [0 , 0 , 0 , 1 , 1 , 1]);
k2.geom_bar
k2.set_color_options('map' , rgb('black'))
k2.axe_property('PlotBoxAspectRatio' , [1 , 1 , 1],...
    'XTick' , [15 , 45 , 75] , 'XLim' , [0 , 90] , ...
    'XTickLabel' , [0 , 45 , 90] , 'YTick' , 0:0.1:0.4 , ...
    'YLim' , [0 , 0.4] );
k2.set_names('x' , 'circular dispersion' , 'y' , 'Fraction' , 'column' , 'region of space')
snapnow;
k2.update('x' , iacaTab3d(:,1) , 'y' , iacaTab3d(: , 2) , 'subset' , [] )
k2.geom_point;
k2.geom_line;
k2.set_color_options('map' , rgb('red'))
k2.draw;
