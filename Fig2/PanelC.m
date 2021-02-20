%% Figure 2 Panel C
% Generates traces of the synaptic efficacies for three examples from the
% different dynamic regimes outlined in panel B
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
normSTD = 6; bAMP = 5;
rho = (2*eta - 1)/(2 - 2*eta);
XDUR = 50;
tauW = tauw/XDUR * (1/(2*beta*(1 - eta))); 
tauY = tauy/XDUR; tauR = taum/XDUR;

T = 864000;
cMap = cbrewer('div' , 'RdBu' , 3);

%%
% high correlation, medium density
L = 64; N = 20; % nu = 20/64 = 0.3125
cMAGNITDE = 0.35; % homogeneous correlations = 0.35
pos = (sort(1:N)*L/(N))'; sPos = mod((pos + L/2),L);
dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
SMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);

% Initialize variables and record variables
R = zeros(N,T); Y = zeros(N,T); W = zeros(N,T); W(:,1) = Winit;
FR = 50*15/1000/60;
cMat = ones(N)* FR*cMAGNITDE; cMat = cMat - eye(N).*cMat + eye(N)*FR;
% generate Poisson input
S =  sampleCovPoisson(ones(N,1)*FR,cMat,T) > 0;
% iterate minimal model
for tt = 2:T
    Sin = S(:,tt); 
    R(: , tt) = R(: , tt-1)*exp(-1./tauR) + phi*Sin*(1 - exp(-1./tauR));
    Waug = repmat(W(: , tt-1)' , [N , 1]).*SMat;
    Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin)*(1 - exp(-1./tauY));
    W(: , tt) =  min(max(W(: , tt-1) + (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) ) , 0),1) ; 
end
% plot with gramm
M = 5; % only a subset of synapes
WSUB = W(1:M , 1:200:end)';% WSUB = WSUB(:);
dT = (1:200:T)*50/1000/60/60;% repmat((1:20:T)*50/1000/60/60 , M , 1)'; dT = dT(:);

figure;
subplot(1,3,1)
plot(dT , WSUB , 'Color' , rgb('black'))
xlim([0 , 4]);
ylim([0 , 1]) 
yticks([0 , 0.5 , 1])
xticks([0 , 2 , 4])
axis square

%%
% low correlation, low density
L = 64; N = 8; % nu = 8/64 = 0.125
cMAGNITDE = 0.15; % homogeneous correlations = 0.35
pos = (sort(1:N)*L/(N))'; sPos = mod((pos + L/2),L);
dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
SMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);

% Initialize variables and record variables
R = zeros(N,T); Y = zeros(N,T); W = zeros(N,T); W(:,1) = Winit;
FR = 50*15/1000/60;
cMat = ones(N)* FR*cMAGNITDE; cMat = cMat - eye(N).*cMat + eye(N)*FR;
% generate Poisson input
S =  sampleCovPoisson(ones(N,1)*FR,cMat,T) > 0;
% iterate minimal model
for tt = 2:T
    Sin = S(:,tt); 
    R(: , tt) = R(: , tt-1)*exp(-1./tauR) + phi*Sin*(1 - exp(-1./tauR));
    Waug = repmat(W(: , tt-1)' , [N , 1]).*SMat;
    Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin)*(1 - exp(-1./tauY));
    W(: , tt) =  min(max(W(: , tt-1) + (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) ) , 0),1) ; 
end

M = 5; % only a subset of synapes
WSUB = W(1:M , 1:20:end)';% WSUB = WSUB(:);
dT = (1:20:T)*50/1000/60/60;% repmat((1:20:T)*50/1000/60/60 , M , 1)'; dT = dT(:);

subplot(1,3,2)
plot(dT , WSUB , 'Color' , rgb('black'))
xlim([0 , 4]);
ylim([0 , 1]) 
yticks([0 , 0.5 , 1])
xticks([0 , 2 , 4])
axis square

%%

%%
% low correlation, high density
L = 64; N = 32; % nu = 32/64 = 0.5
cMAGNITDE = 0.15; % homogeneous correlations = 0.35
pos = (sort(1:N)*L/(N))'; sPos = mod((pos + L/2),L);
dMat = min(pdist2(pos,pos) , pdist2(sPos,sPos));
SMat = normpdf(dMat , 0 , normSTD)*(sqrt(2*pi)*normSTD);

% Initialize variables and record variables
R = zeros(N,T); Y = zeros(N,T); W = zeros(N,T); W(:,1) = Winit;
FR = 50*15/1000/60;
cMat = ones(N)* FR*cMAGNITDE; cMat = cMat - eye(N).*cMat + eye(N)*FR;
% generate Poisson input
S =  sampleCovPoisson(ones(N,1)*FR,cMat,T) > 0;
% iterate minimal model
for tt = 2:T
    Sin = S(:,tt); 
    R(: , tt) = R(: , tt-1)*exp(-1./tauR) + phi*Sin*(1 - exp(-1./tauR));
    Waug = repmat(W(: , tt-1)' , [N , 1]).*SMat;
    Y(: , tt) = Y(: , tt-1)*exp(-1./tauY) + (Waug*Sin)*(1 - exp(-1./tauY));
    W(: , tt) =  min(max(W(: , tt-1) + (1./tauW)*( Y(: , tt-1).*  (R(: , tt-1) + rho) ) , 0),1) ; 
end

M = 5;
WSUB = W(1:M , 1:20:end)';
dT = (1:20:T)*50/1000/60/60;

subplot(1,3,3)
plot(dT , WSUB , 'Color' , rgb('black'))
xlim([0 , 4]);
ylim([0 , 1]) 
yticks([0 , 0.5 , 1])
xticks([0 , 2 , 4])
axis square
