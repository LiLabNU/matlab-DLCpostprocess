clear all
close all
%%
% parameters
framerate = 30; % frame rate of the videos are recorded at
fixedratio = 3; % fixed ratio of the behavioral task
pre = 2; % how many second to include before sucrose delivery
post = 0; % how many second to include after sucrose port entry
blockDiv = 4; % how many blocks to divide for each session
minSuc = 10; % minimal number of sucrose intakes to be included for the analysis

% groups
group1 = ["19479";"19480";"19481";"19482"]; % Females CeA opto; ChR2 good
group2 = ["19510";"19511";"19512";"19513"]; % Females CeA opto; ChR2 bad
group3 = ["19489";"19490";"19494";"19495";"19496";"19497";"19498"]; % Females CeA opto; Control

group4 = ["19475";"19476";"19478"]; % Males CeA opto; ChR2 good
group5 = ["19474";"19486";"19487";"19488"]; % Males CeA opto; ChR2 bad
group6 = ["19506";"19507";"19508";"19509";"19483";"19484";"19485"]; % Males CeA opto; Control

miceToUse = [group3];
%miceToUse = ["378-1";"378-2";"379-1";"379-2";"380-1";"380-2";"381-1";"381-2"];
toPlot = 0;
[avge, Diste, trialTSe] = batchProcessTrajectories(miceToUse,framerate,fixedratio,pre,post,blockDiv,minSuc,toPlot);


miceToUse = [group6];
toPlot = 0;
[avgc, Distc, trialTSc] = batchProcessTrajectories(miceToUse,framerate,fixedratio,pre,post,blockDiv,minSuc,toPlot);


d = trialTSc.np;
finalL = [];
for i = 1:size(d,2)
    a = d(:,i);
    a = a(~isnan(a));
    l = a(end);
    finalL(i,1) = find(a<=l/4,1,'last');
    finalL(i,2) = find(a<=l/2,1,'last');
    finalL(i,3) = find(a<=l/4*3,1,'last');
    finalL(i,4) = find(a<=l,1,'last');
end

