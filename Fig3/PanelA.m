%% Figure 3 Panel A
% Generates an illustration of the retinal wave movie and the LN model
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
% load excerpt from retinal wave movie
load('../data/WMSHORT.mat')
%%
% has been downsampled by factor 2
DOWNSAMPLE = 2;
%%
% set up receptive fields
M = size(WM , 1); N = 250;
ferretSpreadInDegree = 5.3; fullVisualFieldDegree = 120;
muVar = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree)^2;
% size of RF
S1 = 1/15; S2 = 1/8;
% LN parameters
a = 0.01; b = 9.4;
% initial receptive field centers
MUs = mvnrnd([0 , 0],eye(2)*muVar,N);
% intial receptive field orientations
thetas = sort(rand(N,1)*2*pi);
% generate 2D Gabors
RFsize = 13.4; 
[XSPACE , YSPACE] = meshgrid(linspace(-pi , pi , M) , linspace(-pi , pi , M));
fxs = zeros( M , M , N);
for ii = 1:N
    [fx] = mixGauss([XSPACE(:) , YSPACE(:)] , MUs(ii , :) , thetas(ii) , S1 , S2 , RFsize);
    fxs(: , : , ii) = reshape(fx , [M , M]);
end
%%
%%
figure; hold on;
% considered time period
T = 700/DOWNSAMPLE; cMAPTIME = 120/DOWNSAMPLE;
% retinal wave colormaps
cMAPfull = {cbrewer('seq' , 'Purples' , cMAPTIME),...
        cbrewer('seq' , 'Oranges' , cMAPTIME),...
        cbrewer('seq' , 'Greens' , cMAPTIME),...
        cbrewer('seq' , 'Blues' , cMAPTIME)};
% count waves
cc = 1; tt2 = 1;
% construct continuous colormap
cMAPtime = []; cCol = rgb('black');
% time points of waves
WAVEIDS = floor([10:70 , 200:240 , 430:480 , 530:580 ]/DOWNSAMPLE);
% iterate through video
accLINAC = zeros(N , T);
for tt = 1:T
    subplot(1,2,1)
    % if wave occurs
    if sum(tt == WAVEIDS) > 0 
        if pause
            % if new wave starts, switch colormap and reset counter
            cc = mod(cc,4) + 1;
            cMAP = cMAPfull{cc};
            tt2 = 1;
        end
        tt2 = mod(tt2 + 1 , size(cMAP , 1)) + 1;
        cCol = cMAP(tt2 , :);
        % layer frames of wave on top of each other
        alphamask(WM(: , : , tt) > 0.25 , cMAP(tt2 , :) , 0.15);
        cMAPtime = [cMAPtime ; cMAP(tt2 , :)];
        pause = 0;
    else
        pause = 1;
        cMAPtime = [cMAPtime ; 1,1,1];
        cCol = rgb('lightgray');
    end
    subplot(1,2,2); hold on; % get spike trains
    [Sin , sCWM] = getSin(WM(: , : , tt) , fxs , a , b); 
    accLINAC(: , tt) = sCWM;
    spEv = find(Sin);
    for ii = 1:length(spEv)
        if pause
            plot( [tt , tt]*50*2/1000 , [spEv(ii)-2.5 , spEv(ii)+2.5], 'Color' , cCol , 'LineWidth' , 1 )
        else
            plot( [tt , tt]*50*2/1000 , [spEv(ii)-2.5 , spEv(ii)+2.5], 'Color' , cCol , 'LineWidth' , 4 )
        end
    end
end
% improve formatting
subplot(1,2,1)
viscircles([36 , 36] , 36 , 'Color' , rgb('black') , 'EnhanceVisibility' , 0)
axis square
set(gca,'YDir','normal')
caxis(([1,T])*50*2/1000)
colormap(cMAPtime)
colorbar
subplot(1,2,2)
axis square
ylim([0 , 250])
%%
% once more for each individual wave
for tt = 1:T
    if sum(tt == WAVEIDS) > 0 
        if pause
            figure; hold on;
            axis square;
            set(gca,'YDir','normal');
            caxis(([1,T])*50*2/1000);
            colormap(cMAPtime);
            cc = mod(cc,4) + 1;
            cMAP = cMAPfull{cc};
            tt2 = 1;
        end
        tt2 = mod(tt2 + 1 , size(cMAP , 1)) + 1;
        cCol = cMAP(tt2 , :);
        alphamask(WM(: , : , tt) > 0.25 , cMAP(tt2 , :) , 0.15);
        set(gca,'YDir','normal');
        cMAPtime = [cMAPtime ; cMAP(tt2 , :)];
        pause = 0;
    else
        pause = 1;
        cMAPtime = [cMAPtime ; 1,1,1];
        cCol = rgb('lightgray');
    end
end
%%
% plot linear and nonlinear activation traces
cMAP = getWilsonMap(N);
figure; 
subplot(2,1,1)
hold on;
for ii = 1:1:N
    plot(DOWNSAMPLE*(1:4:T)*50/1000 , accLINAC(ii , 1:4:end)' , 'Color' , cMAP(ii , :))
end
axis square;
xlim([0 , 40]); ylim([-0.5 , 0.5]);
xticks([0 , 20 , 40]);yticks([-0.5 , 0 , 0.5])
subplot(2,1,2)
hold on;
for ii = 1:1:N
    plot(DOWNSAMPLE*(1:4:T)*50/1000 , min(a*exp(b*accLINAC(ii , 1:4:end)'),1) , 'Color' , cMAP(ii , :))
end
axis square;
xlim([0 , 40]); ylim([0 , 0.3])
xticks([0 , 20 , 40]); yticks([0 , 0.15 , 0.3])