%% Figure 2 Panel D
% Generates plot of change in synaptic efficacy as a function of
% homogeneous correlations and for two densities (high and low)
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% May 2020;
addpath(genpath('../tools'));
close all
%%
REPS = 100; CLENGTH = 25;
cSpace = linspace(0 , 1.0 , CLENGTH);

accDWHIGH = zeros(CLENGTH , REPS);
accDWLOW = zeros(CLENGTH , REPS);

for xx = 1:length(cSpace)
    cSpace(xx)
    dW = runWithDensityAndCorrelation(0.5 , cSpace(xx) , REPS);
    accDWHIGH(xx , :) = dW;
    dW = runWithDensityAndCorrelation(0.05 , cSpace(xx) , REPS);
    accDWLOW(xx , :) = dW;
end

%%

figure;
g = gramm('x' , repmat(cSpace , [1 , 2*REPS]) , 'y' , [accDWLOW(:)/0.9 ; accDWHIGH(:)/0.9] ,...
    'color' , [zeros(REPS*CLENGTH , 1) ; ones(REPS*CLENGTH , 1)]);
g.set_point_options('base_size' , 3)
g.geom_point('alpha' , 0.25 );
g.stat_summary;
g.axe_property('YLim' , [0 , 1.2] , 'XLim' , [0 , 1] , ...
    'XTick' ,0:0.25:1 , 'YTick' , 0:0.4:1.2 , ...
    'PlotBoxAspectRatio' , [1 , 1 , 1]);
g.set_color_options('map' , [rgb('red') ; rgb('green')]);
g.set_names('x' , 'correlation' , 'y' , 'synaptic efficacy')
g.draw;
