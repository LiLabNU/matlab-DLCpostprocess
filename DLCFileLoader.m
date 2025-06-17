clear all
close all
clc

%% load MedPC and DLC files with animal ID matched
% 1) Set one or multiple mother folders
folderPaths = {
    '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\PBS\LiPatel_Labs\Personal_Folders\Talia\Behavior\MedPC_Data\DOB\Day8'
    };

% Set offset for the delay of video recording from MedPC
MedPCBin = 100; % number of bins for each second in medpc data. 100 would be 0.01s increments in time
MedPC_program = 'Disc'; % use SA, FRnosepoke, or Disc
framerate = 30; % frame rate of the videos are recorded at 30 frames/second


%% Gather all files from each mother folder
Folders   = {};
FileNames = {};

for f = 1:numel(folderPaths)
    [subFolders, subFileNames] = getDirContents(folderPaths{f});
    Folders   = [Folders  subFolders];
    FileNames = [FileNames; subFileNames];
end

% Convert allFolders to a column cell array
Folders = Folders';
FileNames = FileNames';
% remove .txt medpc files starts with ~$
keepMask = cellfun(@(x) ~startsWith(extractAfter(x, find(x == filesep, 1, 'last')), '~$'), FileNames);
FileNames = FileNames(keepMask);

% get file names
MedPCfiles = FileNames(contains(FileNames,'txt'));
DLCfiles = FileNames(contains(FileNames,'filtered.csv'));
dlc2kinematicsfiles = FileNames(contains(FileNames,'dlc2kinematics.csv'));
% get KPMS file names if applicable (containing a folder with kpms in the folder name
isFolder = cellfun(@isfolder, FileNames);
if contains(FileNames(isFolder),'kpms')
    KPMSfolder = string(FileNames(isFolder));
    KPMSfiles = GetFilesFromFolder(0,'filtered.csv',KPMSfolder);
    KPMSfiles = string(KPMSfiles);
end

% load files
switch MedPC_program
    case 'SA'
        [MedPCdata, trialTS, AnimalIDcell] = SA_MedPC2mat(MedPCfiles, folderPaths);  % FRnosepoke_MedPC2mat SA_MedPC2mat
    case 'FRnosepoke'
        [MedPCdata, trialTS, AnimalIDcell] = FRnosepoke_MedPC2mat(MedPCfiles, folderPaths);
    case 'Disc'
        [MedPCdata, ~, AnimalIDcell, trialTypes, trialTS] = PortEntry_MedPC2mat(MedPCfiles, folderPaths); % Pavlovian Disc task
        trialTypes = trialTypes{1};
        trialTS = trialTS{1};
end
MedPCmice = AnimalIDcell(:,1);
% get DLC file names
DLCfiles = string(DLCfiles);
DLCmice = string(regexp(DLCfiles, '_([^_]*?)DLC', 'tokens'));

%%
% Reorder DLC file names based on the order of MedPC data
[MatchedMice, DLCMatchIdx, MedPCkeepIdx] = matchingMice(MedPCmice, DLCmice);
% Apply MedPCkeepIdx to AnimalIDcell
AnimalIDcell = AnimalIDcell(MedPCkeepIdx, :);
% Apply MedPCkeepIdx and Offset delay to all fields in MedPCdata and trialTS
if ~strcmp(MedPC_program,'Disc')   % applicable for FRnosepoke_MedPC2mat SA_MedPC2mat
    tempFields = fieldnames(MedPCdata);
    for i = 1:numel(tempFields)
        MedPCdata.(tempFields{i}) = MedPCdata.(tempFields{i})(MedPCkeepIdx, :);
    end
    tempFields = fieldnames(trialTS);
    for i = 1:numel(tempFields)
        trialTS.(tempFields{i}) = trialTS.(tempFields{i})(:, MedPCkeepIdx);
        for j = 1:numel(trialTS.(tempFields{i}))
            trialTS.(tempFields{i}){j} = trialTS.(tempFields{i}){j};
        end
    end
else    % applicable for PortEntry_MedPC2mat
    MedPCdata = MedPCdata(MedPCkeepIdx,:);
end

%%%%% MedPC mice are now matched with DLC mice, and both data streams start
%%%%% at the same time

% Apply DLCMatchIdx to DLCfiles
DLCfiles = DLCfiles(DLCMatchIdx);
% Load DLCfiles
accuracyThreshold = 0.9;
for i = 1:numel(DLCfiles)
    % Read the DLC data table (skip the first header row)
    DLCdata = readtable(DLCfiles(i), 'NumHeaderLines', 1);
    % Convert DLCdata into the structured format with interpolations
    bodyPartData = convertDLCtable(DLCdata, accuracyThreshold);
    % Add MatchedMice as the first field
    bodyPartData.Mouse = MatchedMice(i);
    allFields = fieldnames(bodyPartData);
    orderedFields = ['Mouse'; setdiff(allFields, 'Mouse', 'stable')];
    bodyPartData = orderfields(bodyPartData, orderedFields);
    % Store the updated structure
    bodyPartStruct(i) = bodyPartData;
end

if ~isempty(dlc2kinematicsfiles)
    clear temp
    dlc2kinematicsfiles  = dlc2kinematicsfiles(DLCMatchIdx);
    for i = 1:numel(dlc2kinematicsfiles)
        % Read the dlc2kinematics data table (skip the first header row)
        dlc2kinematics = readtable(string(dlc2kinematicsfiles(i)),'VariableNamingRule', 'preserve');
        tableStruct = table2struct(dlc2kinematics, 'ToScalar', true);
        temp(i) = Hao_mergeStructs(bodyPartStruct(i), tableStruct);
    end
    bodyPartStruct = temp;
end

if exist('KPMSfolder', 'var')  % Check if the variable exists
    if numel(KPMSfolder) == numel(Folders)
        KPMSmice = string(regexp(KPMSfiles, '_([^_]*?)DLC', 'tokens'));
        [MatchedMice, KPMSMatchIdx, ~] = matchingMice(MedPCmice, KPMSmice);
        KPMSfiles = KPMSfiles(KPMSMatchIdx);
        for i = 1:numel(KPMSfiles)
            warning off
            % Read the DLC data table (skip the first header row)
            KPMSdata = readtable(fullfile(KPMSfolder,KPMSfiles(i)),'VariableNamingRule', 'preserve');
            tableStruct = table2struct(KPMSdata, 'ToScalar', true);
            temp(i) = Hao_mergeStructs(bodyPartStruct(i), tableStruct);
        end
    else
        warning on
        warning('Some folders do not contain kpms results')
    end
    bodyPartStruct = temp;
end

% offset MedPC data for the delayt of video recording starts
offset = 2.1; % offset the delay between MedPC start and camera recording start, in seconds

%%%%% Important: DLC data is in frames; MedPC data is in 10 milliseconds. To
%%%%% align DLC data with MedPC events, you will need to convert them onto
%%%%% the same clock. For example, to convert both of them onto seconds,
%%%%% you should have DLC data / framerate and MedPC data / MedPCBin (100)

%% for talia
blkName = ["Appetitive1";"Aversive1";"Appetitive2";"Aversive2"];
trialTypeNames = ["RewardCue"; "ShockCue"; "NeutralCue"];
blk = {1:23; 24:46; 47:68; 69:90};
trialTS_aj = round((trialTS - offset)*framerate); % in video frames

fldname = "virtual_body_speed"; % in bodyPartStruct
win = [-5 10];
psthAvezScore = struct;
psthAve = struct;
psthMatrix = struct;
for i = 1: numel(bodyPartStruct)
    for b = 1:numel(blk)
        trialTS_aj_blk = trialTS_aj(blk{b});
        trialTypes_blk = trialTypes(blk{b});
        data = bodyPartStruct(i).(fldname)';
        for ty = 1:numel(trialTypeNames)
            if ty == 3 % neutral cue
                ts = trialTS_aj_blk(trialTypes_blk==4);
            else % others
                ts = trialTS_aj_blk(trialTypes_blk==ty);
            end
            [binEdges, psthAvezScore.(trialTypeNames(ty)).(blkName(b))(i,:), psthMatrix.(trialTypeNames(ty)).(blkName(b))(:,:,i), psthAve.(trialTypeNames(ty)).(blkName(b))(i,:)] ...
                = compute_psth_matrix(data, ts, win, framerate, 'AvezScore'); % data here can either be DLC metrics or medpc port entry. 
            % If you use medPC data here, you will have to convert 'ts' into medpc bins
        end
    end
end

% plot
figure
subplot(1,2,1)
h1 = errorbar_pn_hao(psthAvezScore.NeutralCue.Appetitive1(1:7,:),[0, 0.4470, 0.7410]);
hold on
h2 = errorbar_pn_hao(psthAvezScore.NeutralCue.Aversive2(1:7,:),[0.8500, 0.3250, 0.0980]);
legend([h1(1), h2(1)], {'Appe 1','Aver 2'}); % Only on the first subplot
title('exp')
ylim([-2 6])

subplot(1,2,2)
h1 = errorbar_pn_hao(psthAvezScore.NeutralCue.Appetitive1(8:end,:),[0, 0.4470, 0.7410]);
hold on
h2 = errorbar_pn_hao(psthAvezScore.NeutralCue.Aversive2(8:end,:),[0.8500, 0.3250, 0.0980]);
legend([h1(1), h2(1)], {'Appe 1','Aver 2'}); % Only on the first subplot
title('ctrl')
ylim([-2 6])

figure 
h1 = errorbar_pn_hao(psthAvezScore.NeutralCue.Aversive2(1:7,:),[0, 0.4470, 0.7410]);
hold on
h2 = errorbar_pn_hao(psthAvezScore.NeutralCue.Aversive2(8:end,:),[0.8500, 0.3250, 0.0980]);
legend([h1(1), h2(1)], {'exp','ctrl'}); % Only on the first subplot
title('shock cue aver2')
ylim([-2 6])













%% for valen
%% ploting motion trajectory for each trial between nose poke and sucrose intake
framerate = 30; % frame rate of the videos are recorded at 30 frames/second
numBins = []; % force data into x number of bins. input as empty if you want the originals
fr = 3;

for i = 1:numel(bodyPartStruct)
    firstPE = find_first_portentry_MedPC(trialTS.Sucrose{i}, MedPCdata.PortEntry(i,:));
    sucDel = trialTS.Sucrose{i};

    figure
    for j = 1:numel(firstPE)
        if ~isnan(firstPE(j))
            windowRange = sucDel(j):firstPE(j);
            eventFrames = unique(round((windowRange/MedPCBin-offset)*framerate));

            if ~isempty(numBins)
                nose{i,j} = imresize(bodyPartStruct(i).nose(eventFrames,:),[numBins,2],'bilinear');
            else
                nose{i,j} = bodyPartStruct(i).nose(eventFrames,:);
            end
            sucrosePort = mean(bodyPartStruct(i).sucroseport(eventFrames,:),1);
            topNosePort = mean(bodyPartStruct(i).topnoseport(eventFrames,:),1);
            botNosePort = mean(bodyPartStruct(i).botnoseport(eventFrames,:),1);

%             D = [];
%             for k = 1:size(nose{i,j},1)
%                 D(k) = sqrt((nose{i,j}(k,1) - sucrosePort(1)).^2 + (nose{i,j}(k,2) - sucrosePort(2)).^2);
%             end
%             Dewelling{i,j} = length(D) - find(D <= 20,1,'first');
%             Distance{i,j} = D;
        

            col = 8;
            row = ceil(numel(firstPE)/col);
            subplot(row,col,j)
            title('Trial ' + j)
            scatter(sucrosePort(1), sucrosePort(2), 50, 'k', 'filled');
            hold on
            scatter(topNosePort(1), topNosePort(2), 50, 'k', 'filled');
            scatter(botNosePort(1), botNosePort(2), 50, 'k', 'filled');
            % For the nose data, choose a color value for each point.
            cvals = 1:size(nose{i,j},1);  % or any other vector with the same number of rows as nose
            % Now, scatter plot the nose data:
            scatter(nose{i,j}(:,1), nose{i,j}(:,2), 10, cvals, 'filled', 'MarkerFaceAlpha', 0.25);
            % Set the colormap and display a colorbar
            colormap(jet);
            %colorbar
            hold off
            sgtitle(bodyPartStruct(i).Mouse)

            % if j ~= 1
            %     distance{i,:}(im-1,:) = dtw_distance(nose{i,im},nose{i,im-1});
            % end


        end
        % figure
        % bar(distance{i,:})
    end
end

%%
nbins = 36;
breakpoint = [];
H = {};
for i = 1:numel(bodyPartStruct)
    h = [];
    temp = bodyPartStruct(i).nose;
    minX = min(temp(:,1));
    maxX = max(temp(:,1));
    minY = min(temp(:,2));
    maxY = max(temp(:,2));
    for j = 1:size(nose,2)
        if ~isempty(nose{i,j})
            h(j) = path_entropy(nose{i,j}(:,1), nose{i,j}(:,2), nbins, minX, maxX, minY, maxY);
        end
    end
    H{i} = h;

    if ~isempty(H{i})
        % Detect one significant change point
        idx = findchangepts(H{i}, 'MaxNumChanges', 1, 'Statistic', 'mean');
        if ~isempty(idx)
            breakpoint(i) = idx;
        end
    end
end
disp(breakpoint);

%%
for i = 1:numel(bodyPartStruct)
    %deltaD{i} = diff(cell2mat(Dewelling(i,:)));
    for j = 1:size(Distance,2)
        dist = cell2mat(Distance(i,j));
        if ~isempty(dist)
            d(j) = find(dist <= (max(dist) - min(dist))/2,1,'first');
        else
            d(j) = 0;
        end
    end
 
%        idx{i,:} = findchangepts(d(d~=0), 'MaxNumChanges', 1, 'Statistic', 'mean');
        Dist{i,:} = d;
  
end



group1 = ["19479";"19480";"19481";"19482"]; % Females CeA opto; ChR2 good
group2 = ["19510";"19511";"19512";"19513"]; % Females CeA opto; ChR2 bad
group3 = ["19489";"19490";"19494";"19495";"19496";"19497";"19498"]; % Females CeA opto; Control

group4 = ["19475";"19476";"19478"]; % Males CeA opto; ChR2 good
group5 = ["19474";"19486";"19487";"19488"]; % Males CeA opto; ChR2 bad
group6 = ["19506";"19507";"19508";"19509";"19483";"19484";"19485"]; % Males CeA opto; Control


isInGroup = find(ismember(MedPCmice, group3));
breakpoint(isInGroup)
