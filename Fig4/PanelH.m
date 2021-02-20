%% Figure 4 Panel H
% Computes average correlation as a function of the receptive field center
% spread for orientation & overlap clustering
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;
load('../data/WMSHORT.mat')
%%
WM = WM(: , : , 1:10000);
%%
addpath(genpath('../tools'));
% close all
%%
N = 30; L = 150;
scaleSpaceFerret = 0.125:0.125:6;%0.5:0.5:5;
scaleSpaceMouse = ((scaleSpaceFerret)*5.3)/26;

% close all;
if exist('../data/Figure4PanelG.mat')
    load('../data/Figure4PanelH.mat')
else
    % load all retinal wave simulation
    fListFerret = rdir('../Fig4/NO*.mat');

    accAVGCOR = []; 
    accSCALES = [];
    accSIMTYPE = [];

    %%
    % iterate over all simulations
    for xx = 1:length(fListFerret)
        xx
        dat = load(fListFerret(xx).name , '-regexp' , '(dMat)|(thetas)|(MUs)');
        
        [accRFCOR] = getScaleCorrelation(dat , 13.4 , scaleSpaceFerret , WM);
        dMat = dat.dMat - eye(N);
        dMat = dMat./sum(dMat);
        avgCor = accRFCOR'*dMat(:)/N;
        accAVGCOR = [accAVGCOR ; avgCor]; 
        accSCALES = [accSCALES ; 5.3*scaleSpaceFerret']; 
        accSIMTYPE = [accSIMTYPE ; zeros(length(scaleSpaceFerret) , 1)];
    end
    
    fListMouse = rdir('../sims/Fig5_94/NO*.mat');

    %%
    % iterate over all simulations
    for xx = 1:length(fListMouse)
        xx
        dat = load(fListMouse(xx).name , '-regexp' , '(dMat)|(thetas)|(MUs)');
        [accRFCOR] = getScaleCorrelation(dat , 20 , scaleSpaceMouse , WM);
        dMat = dat.dMat - eye(N);
        dMat = dMat./sum(dMat);
        avgCor = accRFCOR'*dMat(:)/N;
        accAVGCOR = [accAVGCOR ; avgCor]; 
        accSCALES = [accSCALES ; 26*scaleSpaceMouse']; 
        accSIMTYPE = [accSIMTYPE ; ones(length(scaleSpaceMouse) , 1)];
    end
    
    %%
    save('../data/Figure4PanelH.mat')
end
%%
ovCor = load('../data/overlapCorrelation.mat');
%%
figure;

g = gramm('x' , accSCALES , 'y' , accAVGCOR , 'color' , accSIMTYPE); % interp1(ovCor.RFOV,ovCor.CORR,accAVGCOR)
% g.geom_point;
g.stat_summary;
g.axe_property('XLim' , [0 , 30] ,...
               'YLim' , [0.0 , 0.15] ,...
               'XTick' , 0:10:30 ,...
               'YTick' , 0.00:0.05:0.15,...
               'PlotBoxAspectRatio' , [1 , 1 , 1]);
g.set_names('x' , 'RF spread' , 'y' , 'average correlation')
g.no_legend;
g.draw;