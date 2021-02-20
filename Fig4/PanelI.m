%% Figure 4 Panel I
% Computes average correlation as a function of the receptive field center
% spread and the diameter for orientation & overlap clustering
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;
load('/home/kirchnerj/synaptic_clustering/master/data/WMSHORT.mat')
%%
WM = WM(: , : , 1:10000);
%%
addpath(genpath('../tools'));
% close all
%%
N = 30; L = 150;
spreadSpaceFerret = 0.125:0.125:6;
spreadSpaceMouse = ((spreadSpaceFerret)*5.3)/26;

scaleSpace = 1:1.25:30;

% close all;
if exist('../data/Figure4PanelI.mat')
    load('../data/Figure4PanelI.mat')
else
    % load all retinal wave simulation
    fListFerret = rdir('/media/kirchnerj@MPIBR/Elements/SIM_OUT/SIM_OUT/REDO_FERRET/Fig4/NO*.mat');

    accAVGCOR = []; 
    accSCALES = [];
    accSIZE = [];
    accSIMTYPE = [];

    %%
    % iterate over all simulations
    for xx = 1:length(fListFerret)
        xx
        dat = load(fListFerret(xx).name , '-regexp' , '(dMat)|(thetas)|(MUs)');
        for xx2 = 1:length(scaleSpace)
            RFsize = scaleSpace(xx2);
            [accRFCOR] = getScaleCorrelation(dat , RFsize , spreadSpaceFerret , WM);
            dMat = dat.dMat - eye(N);
            dMat = dMat./sum(dMat);
            avgCor = accRFCOR'*dMat(:)/N;
            accAVGCOR = [accAVGCOR ; avgCor]; 
            accSCALES = [accSCALES ; 5.3*spreadSpaceFerret']; 
            accSIZE = [accSIZE ; RFsize*ones(length(spreadSpaceFerret) , 1)];
            accSIMTYPE = [accSIMTYPE ; zeros(length(spreadSpaceFerret) , 1)];
        end
    end
    
    fListMouse = rdir('/home/kirchnerj/synaptic_clustering/master/sims/Fig5_94/NO*.mat');

    %%
    % iterate over all simulations
    for xx = 1:length(fListMouse)
        xx
        dat = load(fListMouse(xx).name , '-regexp' , '(dMat)|(thetas)|(MUs)');
        for xx2 = 1:length(scaleSpace)
            RFsize = scaleSpace(xx2);
            [accRFCOR] = getScaleCorrelation(dat , RFsize , spreadSpaceMouse , WM);
            dMat = dat.dMat - eye(N);
            dMat = dMat./sum(dMat);
            avgCor = accRFCOR'*dMat(:)/N;
            accAVGCOR = [accAVGCOR ; avgCor]; 
            accSCALES = [accSCALES ; 26*spreadSpaceMouse']; 
            accSIZE = [accSIZE ; RFsize*ones(length(spreadSpaceFerret) , 1)];
            accSIMTYPE = [accSIMTYPE ; ones(length(spreadSpaceFerret) , 1)];
        end
    end
    
    %%
    save('../data/Figure4PanelI.mat')
end

%%
tab = table( accSCALES , accAVGCOR , accSIZE , accSIMTYPE , ...
    'VariableNames' , {'spread' , 'corr' , 'size' , 'type'} );
statarray = grpstats(tab,{'spread','size','type'},'mean','DataVars','corr');

figure;

s1 = subplot(1,3,1);
c0Mat = reshape(statarray.mean_corr(statarray.type == 0) , [length(scaleSpace) , length(spreadSpaceFerret)]);
p = pcolor((spreadSpaceFerret)*5.3 , scaleSpace , c0Mat);
p.LineStyle = 'None';
shading interp  
axis square
set(gca,'YDir','normal');
colormap(s1 , cbrewer('seq' , 'Greens' , 100))
caxis([0 , 0.3])

s2 = subplot(1,3,2);
c1Mat = reshape(statarray.mean_corr(statarray.type == 1) , [length(scaleSpace) , length(spreadSpaceMouse)]);
p = pcolor((spreadSpaceFerret)*5.3 , scaleSpace , c1Mat);
p.LineStyle = 'None';
shading interp  
axis square
set(gca,'YDir','normal');
colormap(s2 , cbrewer('seq' , 'Purples' , 100))
caxis([0 , 0.3])
s3 = subplot(1,3,3);
p = pcolor((spreadSpaceFerret)*5.3 , scaleSpace , 1.0*(c1Mat > c0Mat));
p.LineStyle = 'None';
ylim([0 , 30])
xlim([0 , 30])

axis square
set(gca,'YDir','normal');
colormap(s3 , flipud(cbrewer('div' , 'PRGn' , 3)))
hold on
RFsize = [20 , 13.4 , 4 , 5 , 2];
V1Area = [3 , 82.5 , 80 , 380 , 1189];
V1diam = sqrt(V1Area);
RFspread = 26*V1diam(1)./V1diam;

animal = {'mouse'  , 'ferret'  , 'rabbit' ,  'cat' , 'macaque'};
hold on
ss = [];
for ii = 1:length(RFsize)
    s = scatter(RFspread(ii) , RFsize(ii) , 'filled');
    ss = [ss , s];
end
