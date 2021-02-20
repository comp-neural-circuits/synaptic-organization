%% Figure 1 Panel D
% Computes expected change in synaptic efficacy for different input rates
% and different intersynaptic distances
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
T = 2400000; 
XDUR = 50;
% derived cosntant for the analytics with minimal model
rho = (2*eta - 1)/(2 - 2*eta);
tauW = tauw * (1/(2*beta*(1 - eta)));
tauY = tauy/XDUR; tauR = taum/XDUR; 
% four position offsets
poffsets = fliplr([0 , 5 , 10 , 15]);
% twenty different input rates
prate = fliplr(linspace(1 , 20 , 20));
% input rates in millisecond^-1
realRate = prate/1000/60;
% set up accumulator for mean change
mVals = zeros(length(poffsets) , length(prate) , 2);
% colormap for plotting
cMAP = flipud(cbrewer('seq' , 'Greys' , 5));

%%
if exist('../data/Figure1PanelD.mat')
    load('../data/Figure1PanelD.mat')
else
    for xx = 1:length(poffsets)
        xx
        for xx2 = 1:length(prate)
            xx2
            % position synapses and compute proximity matrix
            posOffset = poffsets(xx);
            pos = [50, 50 + posOffset]';
            sPos = mod((pos + L/2),L);
            dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
            dMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
            % get firing rate and convert to millisecond
            FR = prate(xx2)/1000/60;
            % set up input array and create smoothed Poisson input for
            % synapse one
            X = zeros(T , 2); 
            X(: , 1) = (movmean(sampleCovPoisson(FR,FR,T) , XDUR)> 0);
            % set up accumulator arrays
            Y = zeros(T , 2); 
            W = ones(T,2)*Winit; W(:) = Winit;
            pBDNF = zeros(T,2); pBDNF(1,:) = pBDNFinit;
            bdnf = zeros(T,2); bdnf(1,:) = bdnfinit;
            mmp9 = zeros(T,2); mmp9(1,:) = mmp9init;
            %%
            % iterate system for duration of simulation. Do not update W.
            for tt = 2:T
                Waug = repmat(W(tt-1,:) , [2 , 1]).*dMat; 
                Y(tt,:) =  Y(tt-1,:) + (1./tauy)*(-Y(tt-1,:) +(Waug*X(tt,:)')' );
                mmp9(tt,:) = mmp9(tt-1,:) + (1./taum)*(-mmp9(tt-1,:) + phi*X(tt-1,:));
                pBDNF(tt,:) = max(pBDNF(tt-1,:) + (1./taup)*(-pBDNF(tt-1,:) + (1-eta)*Y(tt-1,:) - mmp9(tt-1,:).*pBDNF(tt-1,:)) , 0);
                bdnf(tt,:) = min(bdnf(tt-1,:) + (1./taub)*(-bdnf(tt-1,:) + eta*Y(tt-1,:) + mmp9(tt-1,:).*pBDNF(tt-1,:)) , 1);
            end
            % compute expected change given the rate and fixed weights
            mVals(xx , xx2 , :) = mean(beta*bdnf + alpha*pBDNF);
        end
    end
    save('../data/Figure1PanelD.mat')
end
%%
% change in efficacy for the unstimulated synapse
figure; hold on;
for xx = 1:length(poffsets)
    % compute proximity variable
    posOffset = poffsets(xx);
    pos = [50, 50 + posOffset]';
    sPos = mod((pos + L/2),L);
    sMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
    sMat = normpdf(sMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
    % scatter simulation results and analytic prediction
    scatter(prate , (XDUR/tauw)*mVals(xx ,: , 2), 15 , 'MarkerFaceColor' , cMAP(xx , :) , 'MarkerEdgeColor' , 'none' );
    plot(prate , (XDUR/tauW)*XDUR*rho*realRate*Winit*sMat(1,2) , 'Color' , cMAP(xx , :))
end
axis square
xlabel('Input rate (1/min)'); ylabel('bdnf - proBDNF (a.u.)'); 
%%
% change in efficacy for the stimulated synapse
scaleFactor = 0.6;
figure; hold on;
for xx = 1:length(poffsets)
    scatter(prate , (XDUR/tauw)*mVals(xx ,: , 1), 15 ,'MarkerFaceColor' , rgb('white') , 'MarkerEdgeColor' , rgb('black') );
    plot(prate , scaleFactor*(XDUR/tauW)*(phi*Winit*(XDUR*realRate/(tauR+tauY) + (XDUR*realRate).^2) + XDUR*rho*Winit*realRate) , 'Color' , cMAP(xx , :)) 
end 
axis square
xlabel('Input rate (1/min)'); ylabel('bdnf - proBDNF (a.u.)'); 
ylim([0 , 0.6*10^-5])
