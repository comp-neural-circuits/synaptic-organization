%% Figure 1 Panel F
% Generates traces of the neurotrophin model for BTDP stimulation protocol
% for two values of Delta T
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
% Model parameters
Winit = 0.5; tauw = 6000; alpha = -1; beta = 1;
mmp9init = 0; taum = 600; pBDNFinit = 0; eta = 0.45; taup = 5;
bdnfinit = 0; taub = 5; tauy = 300; phi = 3;
L = 150; normSTD = 6; bAMP = 5;
T = 10000; 
% two temporal offsets
offsets = [50 , 1500 , -1500];
for xx = 1:length(offsets)
    % set up BTDP stimualtion protocol
    offset = offsets(xx); 
    X = zeros(T , 1); Iext = zeros(T , 1); 
    Tinit = 2000;
    for kk = 0:9
        X(Tinit+kk:100:Tinit+1010 , 1) = 1;
        Iext(Tinit+offset+kk:100:Tinit+offset+1010 ) = 1;
    end

    %set up accumulators
    Y = zeros(T , 1); W = ones(T,1)*Winit; W(1,:) = Winit;
    pBDNF = zeros(T,1); pBDNF(1,:) = pBDNFinit;
    bdnf = zeros(T,1); bdnf(1,:) = bdnfinit;
    mmp9 = zeros(T,1); mmp9(1,:) = mmp9init;
    % iterate neurotrophin model
    for tt = 2:T
        W(tt,1) = max(W(tt-1,1) + (1./tauw)*(alpha*pBDNF(tt-1,1) + beta*bdnf(tt-1,1)) , 0);
        Y(tt) =  Y(tt-1) + (1./tauy)*(-Y(tt-1) +W(tt,:)*X(tt,:)' + bAMP*Iext(tt));
        mmp9(tt,:) = mmp9(tt-1,:) + (1./taum)*(-mmp9(tt-1,:) + phi*X(tt-1,:));
        pBDNF(tt,:) = max(pBDNF(tt-1,:) + (1./taup)*(-pBDNF(tt-1,:) + (1-eta)*Y(tt-1) - mmp9(tt-1,:).*pBDNF(tt-1,:)) , 0);
        bdnf(tt,:) = min(bdnf(tt-1,:) + (1./taub)*(-bdnf(tt-1,:) + eta*Y(tt-1) + mmp9(tt-1,:).*pBDNF(tt-1,:)) , 1);
    end
    %%
    % scale for plotting purposes
    dT = repmat((1:T)/1000 , 6 , 1)';  
    dVars = [X , Iext , mmp9 , Y , 5*(bdnf - pBDNF) + 0.5 , Winit + (W-Winit)*50 ];
    % plot outcome with gramm
    IDS = categorical([ repmat("Pre" , T , 1) ;repmat("Post" , T , 1) ; repmat("MMP9" , T , 1) ; repmat("Calcium" , T , 1) ; repmat("BDNF minus proBDNF" , T , 1) ; repmat("W" , T , 1)]);
    ROWIDS = categorical([repmat("Pre & Post" , 2*T , 1) ;  repmat("MMP9 & Calcium" , 2*T , 1) ; repmat("BDNF minus proBDNF" , T , 1); repmat("W" , T , 1)]);
    figure;
    g = gramm('x' , dT(:) , 'y' , dVars(:) , 'color' , IDS);
    g.facet_grid(ROWIDS , []);
    g.geom_line;
    g.axe_property('XLim' , [0. , 5.5] , 'YLim' , [0 , 1] , 'XTick' , 1:5 , 'YTick' , 0:0.5:1);
    g.set_color_options('map' , 'brewer_dark');
    g.set_names('x' , 'Time (s)' , 'row' , 'var' , 'color' , 'var' , 'y' , 'magnitude');
    g.set_order_options('row' , categorical(["Pre & Post" ; "MMP9 & Calcium" ; "BDNF minus proBDNF"; "W"]));
    g.draw;
end
%%