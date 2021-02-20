%% Figure 3 Panel B
% Computes and displays the correlations between spike trains generated
% from receptive fields at different delta orientations.
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
RW = load('../data/WMSHORT.mat');
RWmov = RW.WM(: , : , 1:10000);
WN = load('../data/WNSHORT.mat');
WNmov = WN.WNmov;
M = size(WNmov , 1);  T = size(WNmov , 3);

%%
a = 0.01; b = 9.4;  % nonlinearity parameters

% biologically plausible RF parameters
ferretSpreadInDegree = 5.3; fullVisualFieldDegree = 120;
muSTD = (2*pi*ferretSpreadInDegree/fullVisualFieldDegree);

S1 = 1/15; S2 = 1/8;

epoches = 25; 

%%
accFR = [];

thetaSpace = linspace(0 , 2*pi , 100); % different values of theta to simulate
[XSPACE , YSPACE] = meshgrid(linspace(-pi , pi , M) , linspace(-pi , pi , M));
% reference receptive field centered at origin
MUs = [0 , 0];
RFsize = 13.4;
refFX = mixGauss([XSPACE(:) , YSPACE(:)] , MUs , pi , S1 , S2 , RFsize);
refFX = reshape(refFX , [M , M]);

totCRW = nan(epoches , length(thetaSpace));
totCWN = nan(epoches , length(thetaSpace));

if exist('../data/Figure3PanelB.mat')
    load('../data/Figure3PanelB.mat')
else
    for xx = 1:epoches

        fprintf('Sweeping: Theta = ')
        % iterate over all rotations
         for ii = 1:length(thetaSpace)
            fprintf('%0.2f, ' , thetaSpace(ii) )
            % random center of rotated RF
            RFDIST = trandn(-0.8*(pi/muSTD)*ones(1,1),0.8*(pi/muSTD)*ones(1,1))*muSTD;
            RFANGLE = rand(1,1)*2*pi - pi;
            [MUx , MUy] = pol2cart(RFANGLE , RFDIST);
            MUs = [MUx , MUy];

            rotFX = mixGauss([XSPACE(:) , YSPACE(:)] , MUs , thetaSpace(ii) , S1 , S2 , RFsize);
            rotFX = reshape(rotFX , [M , M]);   
            % run input movie stimulation
            fxs = cat(3 , refFX , rotFX);
            S = zeros(T , 2);
            for tt = 1:T
                cWM = RWmov(: , : , tt);
                [Sin , ~] = getSin(cWM , fxs , a , b);
                S(tt , :) = Sin';
            end
            % compute correlations between spike trains and store in matrix
            totCor = corrcoef(smoothdata(S , 'movmean' ,  60));
            totCRW(xx , ii) = totCor(1,2);
            % repeat for white noise stimulation
            for tt = 1:T
                cWM = WNmov(: , : , tt);
                [Sin , ~] = getSin(cWM , fxs , a , b);
                S(tt , :) = Sin';
            end
            totCor = corrcoef(smoothdata(S , 'movmean' ,  60));
            totCWN(xx , ii) = totCor(1,2);
        end
        fprintf('\n')
    end
    save('../data/Figure3PanelB.mat')
end

%%
% Plot the correlation against the delta theta
figure;
g = gramm('x' , 180*thetaSpace/pi - 180 , 'y' , [totCRW ; totCWN] , ...
    'color' , categorical([repmat("retinal wave" , epoches , 1);repmat("white noise" , epoches , 1)]));
g.stat_summary('type' , 'std' , 'geom' , 'area')
g.axe_property('XLim' , [-180 , 180] , 'YLim' , [-0.1 , 0.3] , 'XTick' , [-100 , 100] , 'YTick' , [-0.1 , 0 ,0.1 ,0.2 ,0.3])
g.set_names('x' ,'Delta theta (°)' , 'y' , 'correlation' , 'color' , 'input type')
g.draw

