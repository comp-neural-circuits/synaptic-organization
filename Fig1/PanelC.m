%% Figure 1 Panel C
% Generates traces of the neurotrophin model for one burst stimulation
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;

addpath(genpath('../tools'));
close all
%%
% Model parameters
Winit = 0.5; tauw = 6000; alpha = -1; beta = 1;
mmp9init = 0; taum = 600; pBDNFinit = 0; eta = 0.45; taup = 5;
bdnfinit = 0; taub = 5; tauy = 300; phi = 3;
L = 150; normSTD = 6;
T = 10000; 

% two position offsets, 0 and 8 micron apart
poffsets = [0 , 8];

%%
for xx = 1:length(poffsets)
    % position two synapses on branch of length 150 micron, offset delta p
    posOffset = poffsets(xx);
    pos = [50, 50 + posOffset]';
    % compute proximity matrix
    sPos = mod((pos + L/2),L);
    dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
    dMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
    % set up input bursts, duration 1 second
    X = zeros(T , 2); Iext = zeros(T , 1); 
    Tinit = 2000;
    for kk = 0:49
        X(Tinit+kk:100:Tinit+1010 , 1) = 1;
    end

    % set up empty accumulator arrays
    Y = zeros(T , 2); W = ones(T,2)*Winit;
    pBDNF = zeros(T,2); pBDNF(1,:) = pBDNFinit;
    bdnf = zeros(T,2); bdnf(1,:) = bdnfinit;
    mmp9 = zeros(T,2); mmp9(1,:) = mmp9init;

    % iterate neurotrophin model for duration of simulation
    for tt = 2:T
        W(tt,:) = max(W(tt-1,:) + (1./tauw)*(alpha*pBDNF(tt-1,:) + beta*bdnf(tt-1,:)) , 0);
        Waug = repmat(W(tt-1,:) , [2 , 1]).*dMat; 
        Y(tt,:) =  Y(tt-1,:) + (1./tauy)*(-Y(tt-1,:) +(Waug*X(tt,:)')' + 5*Iext(tt));
        mmp9(tt,:) = mmp9(tt-1,:) + (1./taum)*(-mmp9(tt-1,:) + phi*X(tt-1,:));
        pBDNF(tt,:) = max(pBDNF(tt-1,:) + (1./taup)*(-pBDNF(tt-1,:) + (1-eta)*Y(tt-1,:) - mmp9(tt-1,:).*pBDNF(tt-1,:)) , 0);
        bdnf(tt,:) = min(bdnf(tt-1,:) + (1./taub)*(-bdnf(tt-1,:) + eta*Y(tt-1,:) + mmp9(tt-1,:).*pBDNF(tt-1,:)) , 1);
    end
    %%
    % show accumulators
    figure(1);
    subplot(2 , 2 , 1); hold on
    plot((1:T)/1000 , mmp9(: , 1)/3 , 'Color' , rgb('orange'))
    plot((1:T)/1000 , Y(: , 1) , 'Color' , rgb('purple'))
    ylim([0 , 0.5]); xlim([1.5 , 4.5]); xticks([2 , 3 , 4]); axis square
    subplot(2 , 2 ,2); hold on
    plot((1:T)/1000 , mmp9(: , 2)  , 'Color' , rgb('orange'))
    plot((1:T)/1000 , Y(: , 2) , 'Color' , rgb('purple'))
    ylim([0 , 0.5]); xlim([1.5 , 4.5]); xticks([2 , 3 , 4]); axis square
    subplot(2 , 2 , 3); hold on
    plot((1:T)/1000 , bdnf(: , 1) -  pBDNF(: , 1) , 'Color' , rgb('green'))
%     plot((1:T)/1000 , pBDNF(: , 1) , 'Color' , rgb('red'))
    ylim([-0.15 , 0.15]); xlim([1.5 , 4.5]); xticks([2 , 3 , 4]); axis square
    subplot(2 , 2 , 4); hold on
    plot((1:T)/1000 , bdnf(: , 2) - pBDNF(: , 2) , 'Color' , rgb('green'))
%     plot((1:T)/1000 , pBDNF(: , 2) , 'Color' , rgb('red'))
    ylim([-0.15 , 0.15]); xlim([1.5 , 4.5]); xticks([2 , 3 , 4]); axis square
end