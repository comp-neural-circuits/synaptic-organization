%% Figure 5 Panel B and C
% Generates the illustration of the bAP attenuation constant and of the
% effect of the bAP on a proximal and a distal synapse
%
% Author: Jan H Kirchner
% email: jan.kirchner@brain.mpg.de
% September 2019;

close all;
figure; hold on;
dX = linspace(0 , 500 , 10000);
cMAP = cbrewer('seq' , 'Oranges' , 6); cMAP = cMAP(2:end-1 , :);
conSpace = 25:50:175;
ps = [];
for ii = 1:length(conSpace)
    p = plot(dX , normpdf(dX , 0 , conSpace(ii))*(sqrt(2*pi)*conSpace(ii)) , 'Color' , cMAP(ii,:));
    ps = [ps , p];
end
legend(ps , num2str(conSpace')); legend boxoff
box off;
axis square;
xlabel('Distance to soma (micron)')
ylabel('Attenuation factor (a.u.)')
xticks(-500:250:500)
yticks(0:0.5:1)

%%

dat = load('../data/BAP_FERRET_TREE_BIOLOGICAL_737573.30009299970697611570.mat');
%%
% pick suitable synapses
somDVec = dat.somDVec;
lCal = find((somDVec > prctile(somDVec , 5)) & (somDVec < prctile(somDVec , 20)) , 1 , 'first');
hCal = find((somDVec > prctile(somDVec , 90)) & (somDVec < prctile(somDVec , 100)) , 1 , 'first');
T = length(dat.Aexcerpt);
figure;
subplot(3,1,1); hold on;
plot((1:T)*50/1000/60 , dat.Aexcerpt)
plot([0 , 8] , [dat.bAP_Thresh , dat.bAP_Thresh])
axis square
xlim([3 , 6]); 
xlabel('Time (min)'); ylabel('som activation (a.u.)')
xticks(3:6)
yticks([0 , 50 , 100])

subplot(3,1,2); hold on;
plot((1:T)*50/1000/60 , dat.Yexcerpt(lCal , :))
ylim([0 , 3])
axis square
xlim([3 , 6]); 
xlabel('Time (min)'); ylabel('calcium (a.u.)')
xticks(3:6)
yticks([0 , 1 , 2 , 4])
dat.compSomDist(dat.pos(lCal))
subplot(3,1,3); hold on;
plot((1:T)*50/1000/60 , dat.Yexcerpt(hCal , :))
ylim([0 , 3])
axis square
xlim([3 , 6]); 
xlabel('Time (min)'); ylabel('calcium (a.u.)')
xticks(3:6)
yticks([0 , 1 , 2 , 4])
dat.compSomDist(dat.pos(hCal))