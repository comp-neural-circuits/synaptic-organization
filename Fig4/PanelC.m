%% Figure 4 Panel C
% Generates the scatter plot in the bottom row of Figure 4 panel C
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;
addpath(genpath('../tools'));
close all
% load experimental data
schollFig2H = csvread('../data/SchollFigure2H.csv'); schollFig2H = sortrows(schollFig2H , 1);
iacarausoFig2B = csvread('../data/IacarusoFigure2b.csv'); iacarausoFig2B = sortrows(iacarausoFig2B , 1);

close all;
% plot and fit a Gaussian
figure;

subplot(1,2,1); hold on;
f = fit(schollFig2H(:,1) , schollFig2H(:,2) , 'a*exp(-(x^2/(2*c^2)))');
scatter(schollFig2H(:,1) , schollFig2H(:,2) , 'MarkerFaceColor' , rgb('green') , 'MarkerEdgeColor' , rgb('black'))
plot(linspace(0  , 80) , f(linspace(0 , 80))  , 'Color' , rgb('black'))

xlim([0 , 80]); xticks(0:40:80);
ylim([0 , 0.5]); yticks([0 , 1])
axis square
box off
xlabel('RF offset')
ylabel('pdf')

subplot(1,2,2); hold on;
f2 = fit(iacarausoFig2B(:,1) , iacarausoFig2B(:,2) , 'a*exp(-(x^2/(2*c^2)))' , 'StartPoint', [0.5 , 26]);
scatter(iacarausoFig2B(:,1) , iacarausoFig2B(:,2) , 'MarkerFaceColor' , rgb('purple') , 'MarkerEdgeColor' , rgb('black'))
plot(linspace(0  , 80) , f2(linspace(0 , 80)) , 'Color' , rgb('black'))

xlim([0 , 80]); xticks(0:40:80);
ylim([0 , 0.5]); yticks([0 , 1])
axis square
box off
xlabel('RF offset')
ylabel('pdf')