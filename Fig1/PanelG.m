%% Figure 1 Panel G
% Computes and displays change in synaptic efficacy for different values of
% Delta T. Superimposes experimental data from the developing
% retinogeniculate nucleus.
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
L = 150; normSTD = 6; bAMP = 5;
T = 400000; 
%%
% 100 temporal offsets
doffset = -2000:100:2000;
% set up accumulator for total change in weight
dWs = zeros(length(doffset),1);
for ii = 1:length(doffset)
    % set up BTDP protocol
    offset = doffset(ii);
    X = zeros(T + 10000 , 1); 
    Iext = zeros(T + 10000 , 1); 
    for tt = 12001:40*1000:T
        for kk = 0:11
            X(tt+kk:100:tt+kk+900 , 1) = 1;
            Iext(tt+offset+kk:100:tt+offset+kk+900 ) = 5;
        end
    end
    % set up accumulator
    Y = zeros(T , 1); W = ones(T,1)*Winit; W(1,:) = Winit;
    pBDNF = zeros(T,1); pBDNF(1,:) = pBDNFinit;
    bdnf = zeros(T,1); bdnf(1,:) = bdnfinit;
    mmp9 = zeros(T,1); mmp9(1,:) = mmp9init;
    % iterate neurotrophin model
    for tt = 2:T
        W(tt,1) = max(W(tt-1,1) + (1./tauw)*(alpha*pBDNF(tt-1,1) + beta*bdnf(tt-1,1)) , 0);
        Y(tt) =  Y(tt-1) + (1./tauy)*(-Y(tt-1) +W(tt,:)*X(tt,:)' + Iext(tt));
        mmp9(tt,:) = mmp9(tt-1,:) + (1./taum)*(-mmp9(tt-1,:) + phi*X(tt-1,:));
        pBDNF(tt,:) = max(pBDNF(tt-1,:) + (1./taup)*(-pBDNF(tt-1,:) + (1-eta)*Y(tt-1) - mmp9(tt-1,:).*pBDNF(tt-1,:)) , 0);
        bdnf(tt,:) = min(bdnf(tt-1,:) + (1./taub)*(-bdnf(tt-1,:) + eta*Y(tt-1) + mmp9(tt-1,:).*pBDNF(tt-1,:)) , 1);
    end
    % compute percentage change in weight
    dWs(ii) = 100*(W(end , 1) - W(1 , 1))/Winit;
end

%% 
% plot with gramm
% ButtsFig3 = csvread('../data/ButtsFigure3.csv'); ButtsFig3 = sortrows(ButtsFig3 , 1);
figure; 

g = gramm('x' , doffset'/1000  , 'y' , dWs  , 'color' , categorical(repmat("sim" , length(doffset) , 1) ));
g.geom_line
g.set_color_options('map' , rgb('black') )
% snapnow;
% g.update('x' , ButtsFig3(: , 1)/1000 , 'y' , ButtsFig3(: , 2) , 'color' , categorical(repmat("data" , length(ButtsFig3(:,1)) , 1)))
% g.geom_point;
g.set_color_options('map' , 'brewer_dark');
g.set_names('x' , 'Temporal offset Delta T' , 'color' , 'data type' ,...
    'y' , '  change in EPSC size or efficacy Wi (%)');
g.draw
