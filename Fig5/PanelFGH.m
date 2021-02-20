%% Figure 5 Panel F, G and H
% Computes and displays the circular dispersion and receptive field offset
% as a function of the distance from the soma
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;


addpath(genpath('../tools'));
close all
%%
if exist('../data/Figure5PanelFGH.mat')
    load('../data/Figure5PanelFGH.mat')
else
    % get all simulations with morphology and small receptive field center
    % spread
    fList = rdir('../sims/Fig5/BAP*.mat');
    if isempty( fList )
        fprintf('run batch_process_ferret_morphological.m first\n')
    end
    N = 504;

    CIRCDISP = zeros(length(fList) , N); 
    MUDIST = zeros(length(fList) , N); 
    POS = zeros(length(fList) , N); 
    somCONST = zeros(length(fList) , N);
    muVAR = zeros(length(fList) , N);
    for xx = 1:length(fList)
        cFile = fList(xx).name
        dat = load(cFile , '-regexp' , '(compSomDist)|(pos)|(subpos)|(thetas)|(MUs)|(somConst)|(muVar)');
        % compute circular dispersion
        thetas =  dat.thetas - pi;
        circ_dist_to_mean = 90*(abs(circ_dist2(thetas , circ_mean(thetas))))/pi ;
        CIRCDISP(xx , :) = circ_dist_to_mean;
        % compute receptive field center offset
        MUDIST(xx , :) = 62.5*sqrt(sum(dat.MUs.^2 , 2))/pi;
        % accumulate distances and bAP attenuation constants
        POS(xx , :) = dat.compSomDist(dat.pos) + dat.subpos;
        somCONST(xx , :) = dat.somConst;
    end
end
%%
cMAP = cbrewer('seq' , 'Oranges' , 6); cMAP = cMAP(2:end-1 , :);

figure; 
g = gramm('x' , POS(:) , 'y' , CIRCDISP(:) , 'color' , somCONST(:) );
g.stat_summary('bin_in' , 10 , 'geom' , 'area');
g.set_color_options('map' , cMAP);
g.axe_property('PlotBoxAspectRatio' , [1 , 1 , 1] , 'YLim' , [0 , 50] ,...
    'XTick' , [0 , 175 , 350] , 'XLim' , [0 , 350] , 'YTick' , [0 , 15 , 30 , 45]);
g.set_names('x' , 'Distance to soma' , 'y' , 'Circular dispersion');
g.draw;
%%
figure; 
h = gramm('x' , POS(:) , 'y' , MUDIST(:) , 'color' , somCONST(:) );
h.stat_summary('bin_in' , 10 , 'geom' , 'area');
h.set_color_options('map' , cMAP );
h.axe_property('PlotBoxAspectRatio' , [1 , 1 , 1] , 'YLim' , [0 , 25] , ...
    'XTick' , [0 , 175 , 350] , 'XLim' , [0 , 350] , 'YTick' , [0 , 12.5 , 25]);
h.set_names('x' , 'Distance to soma' , 'y' , 'RF center offset');
h.draw;
%%

figure; hold on;
uCONST = unique(somCONST(:));
ps = [];
for ii = 1:length(uCONST)
    [FYESBAP,XYESBAP] = ecdf([CIRCDISP(somCONST(:) == uCONST(ii)  & POS(:) < 2*uCONST(ii)) ; 90]);
    p = plot(XYESBAP , FYESBAP , 'Color' , cMAP(ii , :)); ps = [ps , p];
end
for ii = 1:length(uCONST)
    [FYESBAP,XYESBAP] = ecdf([CIRCDISP(somCONST(:) == uCONST(ii)  & POS(:) > 2*uCONST(ii)) ; 90]);
    p = plot(XYESBAP , FYESBAP , 'Color' , cMAP(ii , :) , 'LineStyle' , '--'); ps = [ps , p];
end
legend(ps , {'Distance < 150 micron' , 'Distance > 150 micron'  } , 'Location' , 'southeast' ) ; legend boxoff
xlabel('circular dispersion (°)')
ylabel('cumulative fraction')
xticks([0 , 45 , 90])
yticks([0 , 0.5 , 1])
axis square
colormap(cMAP(1:4 , :))
colorbar