%% Figure 3 Panel E
% Computes and displays the number of turnovers per day of the simulation
% as well as the survival and stable fraction
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;

addpath(genpath('../tools'));
close all
%%
%

N = 30; % number of synapses


close all;
if exist('../data/Figure3PanelE.mat')
    load('../data/Figure3PanelE.mat')
else
    % get a list of all the finished simulations
    fListNO = rdir('../sims/Fig3/NO*.mat');
    if isempty( fListNO )
        fprintf('run batch_process_ferret.m first\n')
    end
    % 15 days of the simulation
    KK = 15 + 1;
    dT = linspace(0 , 15 , KK);
    adT = linspace(0 , 15 , 100*KK);
    regIDXsNOBAP = cell(1 , 1);
    INITTIME = []; INITFRAC = []; FINALTIME = []; FINALFRAC = []; 
    % for each simulation
    for ii = 1:length(fListNO)
        fListNO(ii).name
        % load the turnover times
        dat = load(fListNO(ii).name , '-regexp' , '(storeConfig)|(regIDX)');
        regIDXsNOBAP{end + 1} = dat.regIDX; % accumulate turnover times
        % determine turnover times for each individual synapse
        timeToTO = abs(diff(dat.storeConfig.thetas,1,2)) > 0;
        timeToTO = timeToTO * diag(dat.regIDX);
        for kk = length(dat.regIDX)-1:-1:1
            for jj = 1:N
                if timeToTO(jj , kk) == 0
                    timeToTO(jj , kk) = timeToTO(jj , kk + 1);
                end
            end
        end
        % first turnover of one of the initial 30 synapses
        timeToTOInit = timeToTO(: , 1); 
        timeToTOInit(timeToTOInit == 0 ) = []; % don't count synapses that never turnover
        [timeToTOInit  ,~] = sort(timeToTOInit);
        fracAtTime = 1 - (1:length(timeToTOInit))/N; % compute the survival fraction for different time points
        INITTIME = [INITTIME; 0 ; timeToTOInit];  INITFRAC = [INITFRAC ; 1 ; fracAtTime']; % store them
        finalTO = sum(1 - (timeToTO == 0)'); finalTO(finalTO == 0) = []; % same for stable fraction
        FINALTIME = [FINALTIME ;  sort(dat.regIDX(finalTO)')]; FINALFRAC = [FINALFRAC ; fliplr(1 - (0:length(finalTO)-1)/N)'];
    end
    regIDXsNOBAP = regIDXsNOBAP(2:end);
    save('../data/Figure3PanelE.mat')
end
%%
% plotting
figure; hold on;
cMAP2 = cbrewer('seq' , 'Blues' , 5);
mNOBAP = cell2mat(regIDXsNOBAP); 
totTurnoverNO = cellfun(@length , regIDXsNOBAP);

h1 = histcounts(mNOBAP , dT);
bar((dT(1:end-1)) , h1./length(fListNO) , 'BarWidth' , 0.75 , 'FaceColor' , cMAP2(end-2,:) , 'LineStyle' , 'none')
f = fit(dT(1:end-1)' , h1'./length(fListNO) , 'a/(x + b)^c' , 'Lower' , [0 , 0 , 1]);
plot((adT) , f(adT) , 'LineWidth' , 2 , 'Color' , rgb('black'))
xlabel('Simulated time (days)')
ylabel('Turnovers per 6 hours')
xticks([0 , 5 , 10 , 15])
yticks(0:10:40);%[0 , 50 , 100 , 150])
axis square
xlim([-1 , 16])
%%
figure;
g = gramm('x' , ([FINALTIME ; INITTIME ]) , ...
                'y' , [FINALFRAC ; INITFRAC ], ...
                'color' , categorical([repmat("stable fraction" , length(FINALTIME) ,1); ...
                                           repmat("survival fraction" , length(INITTIME),1) ]));
g.set_color_options('map' , 'd3_10' )
g.stat_summary('type' , 'std' ,  'geom' , 'area','bin_in' , 17 , 'bin_edges' , [0 , 0.04 , 0.08 , 0.16 , 0.32 , 0.64 , 1:5:15] )
g.set_names('x' , 'Time (days)' , 'y' , 'Fraction' , 'color' , 'simulation type' , 'linestyle' , 'type of fraction')
g.axe_property('XLim' , [0  ,15] , 'YLim' , [0 , 1] , 'YTick' , [0 , 0.5 , 1] , 'XTick' , [0 , 0.1 , 1 , 10 ] , 'XScale' , 'Log')
g.draw