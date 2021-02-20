%% Figure 2 Panel B
% Generates phase plane of expected changes in synaptic efficacy for
% different combinations of synaptic density and homogeneous input
% correlations
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
L = 30; normSTD = 6; bAMP = 5;
rho = (2*eta - 1)/(2 - 2*eta);
XDUR = 50;
tauW = 3000/XDUR * (1/(2*beta*(1 - eta))); 
tauY = tauy/XDUR; tauR = taum/XDUR;
FR = XDUR*15/1000/60;
kappa = (-rho/phi - FR)*(tauY + tauR);
T = 14400; 
% vary the density
denSpace = linspace(1 , 15 , 30)./(sqrt(2*pi)*normSTD);
% vary the pairwise homogeneous correlation
cSpace = fliplr(linspace(0 , 0.4 , 30));
% number of epochs
X = 50;
WUWL = zeros(length(cSpace) , length(denSpace) , X);

if exist('../data/Figure2PanelB.mat')
   load('../data/Figure2PanelB.mat')
else
    for xx3 = 1:X
        for xx2 = 1:length(cSpace)
            for xx = 1:length(denSpace)
                % compute the appropriate number of synapses N=L*nu
                den = denSpace(xx); N = floor(L*den);
                pos = sort(rand(N,1))*L;
                sPos = mod((pos + L/2),L);
                dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
                dMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);
                % set magnitude of correlation
                cMAGNITDE = cSpace(xx2);
                fprintf('Epoch : %d/%d , corr : %0.2f , den : %0.2f \n' , xx3 , X , cMAGNITDE , den)
                % initialize accumulators
                R = zeros(N,T); Y = zeros(N,T); W = zeros(N,T); dW = zeros(N,T); W(:) = Winit;
                % generate input traces
                cMat = ones(N)* FR *cMAGNITDE; 
                cMat = cMat - eye(N).*cMat + eye(N)*FR;
                S =  sampleCovPoisson(ones(N,1)*FR,cMat,T) > 0;
                % iterate minimal model
                for tt = 2:T
                    Sin = S(:,tt);
                    R(: , tt) = R(: , tt-1)*exp(-1./tauR) + phi*Sin*(1 - exp(-1./tauR));
                    Waug = repmat(W(: , tt-1)' , [N , 1]).*dMat;
                    Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin)*(1 - exp(-1./tauY)); 
                    dW(: , tt) =  (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) );
                end


                WUWL(xx2 , xx , xx3) = mean(dW(:)); 
            end
        end
    end
    save('../data/Figure2PanelB.mat')
end
%%
% average over epochs
mWUWL = nanmean(WUWL , 3); 
%%
% scale factor
dWMAX = 0.003;
figure; hold on;
imagesc(denSpace , cSpace , mWUWL(: , : , 1)*tauW/dWMAX);  
axis square;
set(gca,'YDir','normal')
xlabel('synapse density')
ylabel('Pairwise correlation \gamma_{i,j}')
aNspace = linspace(1 , 15 , 100)./(sqrt(2*pi)*normSTD);

% compute analytic contours
KK = 1000; % resolution
aCSpace = linspace(0 , 0.4, KK); sSpace = linspace(1 , 15 , KK); dWAna = zeros(KK , KK );
for jj = 1:KK
    S = sSpace(jj);
    for ii = 1:KK
        c = aCSpace(ii);
        dWAna(ii , jj) = ((S-1)*phi*(c/(tauR + tauY) + FR) + phi*(1/(tauR + tauY) + FR) + rho*S)*FR*Winit;
    end
end
% generate contour plot and asymptote
clSpace =  round(linspace(-1 , 1 , 11) , 4); clSpace(6) = [];
plot(aNspace , (kappa.*aNspace*(sqrt(2*pi)*normSTD) - 1)./(aNspace*(sqrt(2*pi)*normSTD) - 1) , 'LineWidth' , 2 , 'LineStyle' , '-.' ,  'Color' , rgb('black'))
[C , h] = contour(sSpace./(sqrt(2*pi)*normSTD) , aCSpace , dWAna/dWMAX ,...
    clSpace , 'ShowText' , 'on' , 'LabelSpacing' , 250 , 'LineWidth' , 2 , 'LineColor' , rgb('white')); 
clabel(C , h , 'Color' , rgb('white'))
caxis([-1.2 , 1.2])
colormap([flipud(cbrewer('div' , 'RdBu' , 100))]); 
cb = colorbar; cb.Label.String = 'expected dW (a.u.)';
xlim([min(denSpace), max(denSpace)])
ylim([min(cSpace) , max(cSpace)])
