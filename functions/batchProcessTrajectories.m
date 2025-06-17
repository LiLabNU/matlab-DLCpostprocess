function [avg, Dist, ts] = batchProcessTrajectories(miceToUse,framerate,fixedratio,pre,post,blockDiv,minSuc,toPlot)

MedPCBin = 100; % number of bins for each second in medpc data

% get files
[MedPCfiles, folderDir] = GetFilesFromFolder(0,'.txt');
[MedPCdata, trialTS, AnimalIDcell] = FRnosepoke_MedPC2mat(MedPCfiles, folderDir);

[DLCfiles, ~] = GetFilesFromFolder(0,'filtered.csv',folderDir{1});
DLCfiles = string(DLCfiles);
pattern = '_([^_]*?)DLC';
matches = string(regexp(DLCfiles, pattern, 'tokens'));

m = 0;
blockTrial = struct();
numRows = max(cell2mat(cellfun(@length,trialTS.Sucrose,"UniformOutput",false)));
ts.suc = nan(numRows,length(miceToUse));
numRows = max(cell2mat(cellfun(@length,trialTS.ActiveNP,"UniformOutput",false)));
ts.np = nan(numRows,length(miceToUse));

for i = 1:size(MedPCfiles,1)
    mouse = (AnimalIDcell(i,1));
    if ismember(mouse,miceToUse)
        m = m+1;
        eventTimestamps = trialTS.Sucrose{i};
        if ~isempty(eventTimestamps)
            [PEdetails, ~] = findBouts(MedPCdata.PortEntry(i,:), MedPCBin, eventTimestamps, AnimalIDcell(i,1), 'mean');
            [FRdetails,~] = findFRduration(trialTS.FirstActive{i}, eventTimestamps, trialTS.ActiveNP{i}, fixedratio);
        end

        ratio = framerate/MedPCBin;
        timestamps.firstNP = ceil(trialTS.FirstActive{i}*ratio); 
        timestamps.firstNPofCompleteFR = ceil(FRdetails.firstCompletedFRTS*ratio);
        timestamps.sucroseDel = ceil(trialTS.Sucrose{i}*ratio);
        timestamps.firstPEntry = ceil(PEdetails.firstBoutOnset*ratio);
        timestamps.firstPExit = ceil((PEdetails.firstBoutOnset + PEdetails.firstBoutDuration)*ratio);
        timestamps.NP = trialTS.ActiveNP{i}*ratio;

        ts.suc(1:length(trialTS.Sucrose{i}),m) = trialTS.Sucrose{i}./MedPCBin;
        ts.np(1:length(trialTS.ActiveNP{i}),m) = trialTS.ActiveNP{i}./MedPCBin;

        % find the DLC file based on the mouse ID
        DLCdata = [];

        match = find(matches == mouse);
     %   mouse = removeSpecialCharacters(mouse);
        fieldname(m,:) = mouse;
        blockTrial.(fieldname(m)) = [];
        if toPlot
            figure
            sgtitle(mouse)
        end

        if ~isempty(match)            
                if length(eventTimestamps)>minSuc
                    DLCfile = DLCfiles(match);
                    if length(match)~=1
                        warning('More than 1 DLC file found for mouse ' +mouse);
                    else
                        disp('Reading DLC file: ' +DLCfile);
                    end                    
                    DLCdata = readtable(DLCfile, 'NumHeaderLines', 1);

                    object.sucrosePort(1) = mean(DLCdata.sucroseport);
                    object.sucrosePort(2) = mean(DLCdata.sucroseport_1);
                    object.topNosePort(1) = mean(DLCdata.topnoseport);
                    object.topNosePort(2) = mean(DLCdata.topnoseport_1);
                    object.botNosePort(1) = mean(DLCdata.botnoseport);
                    object.botNosePort(2) = mean(DLCdata.botnoseport_1);


                    cmap = jet(length(timestamps.sucroseDel)); % Using the 'jet' colormap
                    step = ceil(length(timestamps.sucroseDel)/blockDiv);
                    blocks = 1:step:length(timestamps.sucroseDel);

                    for b = 1:blockDiv
                        if b == 1       
                            timespent(m,b) = timestamps.sucroseDel(blocks(b)+step-1)./framerate;% in seconds;
                            FRduration(m,b) = nanmean(FRdetails.FRDuration(timestamps.firstNPofCompleteFR...
                                <= timestamps.sucroseDel(blocks(b)+step-1)))./MedPCBin; % in seconds;
                            FRcompletion(m,b) = nanmean(FRdetails.FRcompletion(timestamps.firstNP...
                                <= timestamps.sucroseDel(blocks(b)+step-1))); % % binary;
                            numNP(m,b) = nansum(timestamps.NP <= timestamps.sucroseDel(blocks(b)+step-1));
                            numSuc(m,b) = nansum(timestamps.sucroseDel <= timestamps.sucroseDel(blocks(b)+step-1));
                            
                        elseif b==blockDiv
                            timespent(m,b) = 3600; % in seconds;
                            FRduration(m,b) = nanmean(FRdetails.FRDuration(timestamps.firstNPofCompleteFR...
                                >= timestamps.sucroseDel(blocks(b))))./MedPCBin; % in seconds;
                            FRcompletion(m,b) = nanmean(FRdetails.FRcompletion(timestamps.firstNP...
                                >= timestamps.sucroseDel(blocks(b)))); % % binary;
                            numNP(m,b) = nansum(timestamps.NP >= timestamps.sucroseDel(blocks(b)));
                            numSuc(m,b) = nansum(timestamps.sucroseDel >= timestamps.sucroseDel(blocks(b)));                           
                        else
                            timespent(m,b) = timestamps.sucroseDel(blocks(b)+step-1)./framerate;% in seconds;
                            FRduration(m,b) = nanmean(FRdetails.FRDuration(timestamps.firstNPofCompleteFR >= timestamps.sucroseDel(blocks(b))...
                                & timestamps.firstNPofCompleteFR <= timestamps.sucroseDel(blocks(b)+step-1)))./MedPCBin; % in seconds;
                            FRcompletion(m,b) = nanmean(FRdetails.FRcompletion(timestamps.firstNP >= timestamps.sucroseDel(blocks(b))...
                                & timestamps.firstNP <= timestamps.sucroseDel(blocks(b)+step-1))); % % binary;
                            numNP(m,b) = nansum(timestamps.NP >= timestamps.sucroseDel(blocks(b)) & timestamps.NP <= timestamps.sucroseDel(blocks(b)+step-1));
                            numSuc(m,b) = nansum(timestamps.sucroseDel >= timestamps.sucroseDel(blocks(b)) & timestamps.sucroseDel <= timestamps.sucroseDel(blocks(b)+step-1));
                            
                        end
                    end

                  
    

                    b = 1;
                    counter = 0;
                    temp = [];
                    if toPlot
                        subplot(2,ceil(blockDiv/2),b)
                    end
                    for j = 1:length(timestamps.sucroseDel)
                        animal = [];
                        counter = counter+1;
                        %firstNP = timestamps.firstNPofCompleteFR(j);
                        suc = timestamps.sucroseDel(j);
                        Pentry = timestamps.firstPEntry(j);
                        Pexit = timestamps.firstPExit(j);

                        range = (suc-pre*framerate: Pentry+post*framerate);
                        if ~isempty(range) && length(DLCdata.nose) > range(end)
                            animal.nose(:,1) = DLCdata.nose(range);
                            animal.nose(:,2) = DLCdata.nose_1(range);
                            animal.butt(:,1) = DLCdata.tailbase(range);
                            animal.butt(:,2) = DLCdata.tailbase_1(range);

                            latency{m,b}(counter) = (Pentry-suc)/framerate; % in seconds
                            distances{m,b}(counter) = sum(sqrt(diff(animal.nose(:,1)).^2 + diff(animal.nose(:,2)).^2)); % in pixels
                            speed{m,b}(counter) = distances{m,b}(counter)./latency{m,b}(counter); % in pixels/s                         
                           
                            nearSucP{m,b}(counter) = isNear(animal.nose, object.sucrosePort, animal.butt);                            
                            nearNP{m,b}(counter) = isNear(animal.nose, object.botNosePort, animal.butt);                             

                            blockTrial.(fieldname(m)){counter,b}(:,1) = animal.nose(:,1);
                            blockTrial.(fieldname(m)){counter,b}(:,2) = animal.nose(:,2);


                            if ismember(j,blocks(2:end))
                                if toPlot
                                    scatter(object.sucrosePort(1),object.sucrosePort(2), 100, "filled", "k")
                                    scatter(object.topNosePort(1),object.topNosePort(2), 100, "filled", "k")
                                    scatter(object.botNosePort(1),object.botNosePort(2), 100, "filled", "k")
                                    title('Block ' + string(b))
                                end
                                blockData{m,b} = temp;
                                b = b+1;
                                counter = 0;
                                temp = [];
                                if toPlot
                                    subplot(2,ceil(blockDiv/2),b)
                                end
                            end
                            temp = [temp; animal.nose];
                            if toPlot
                                scatter(animal.nose(:,1),animal.nose(:,2), 10, cmap(j,:), "filled","MarkerFaceAlpha","0.25");
                                hold on
                            end
                        end
                    end
                    if j == length(timestamps.sucroseDel)
                        blockData{m,b} = temp;
                        if toPlot
                            title('Block ' + string(b))
                            scatter(object.sucrosePort(1),object.sucrosePort(2), 100, "filled", "k")
                            scatter(object.topNosePort(1),object.topNosePort(2), 100, "filled", "k")
                            scatter(object.botNosePort(1),object.botNosePort(2), 100, "filled", "k")
                        end
                    end
                else
                    disp(['Sucrose intakes less than ' num2str(minSuc) ' for mouse ' num2str(mouse)]);
                end            
        else
            disp('No matched DLC file found for mouse ' +mouse);
        end
    end
end

avg.latency = cellfun(@mean, latency);
avg.distances = cellfun(@mean, distances);
avg.speed = cellfun(@mean, speed);
avg.nearNP = cellfun(@mean, nearNP);
avg.nearSucP = cellfun(@mean, nearSucP);
avg.FRduration = FRduration;
avg.FRcompletion = FRcompletion;
avg.numNP = numNP;
avg.numSuc = numSuc;
avg.timespent = timespent;


% for i = 1:size(MedPCfiles,1)
%     for j = 1:size(blockData,2)
%         [~, R2(i,j)] = fitLinear(blockData{i,j}(:,1), blockData{i,j}(:,2));
%
%         %         % Display the results
%         %         fprintf('Linear Fit Coefficients: \n');
%         %         disp(p);
%         %         fprintf('R^2: %f\n', R2);
%     end
% end

%%
%
% for i = 1:size(MedPCfiles,1)
%     %     reference = blockData{i,1}; % Use the first session as a reference
%     %     standardLength = size(reference,1);
%     %     % Interpolate all trajectories to this length
%     %     for j = 1:blockDiv
%     %         oldX = linspace(0, 1, size(blockData{i,j}, 1));
%     %         newX = linspace(0, 1, standardLength);
%     %         uniformXYData{i,j} = [interp1(oldX, blockData{i,j}(:,1), newX)', interp1(oldX, blockData{i,j}(:,2), newX)'];
%     %         [~, xygrid{i,j}(:,1), xygrid{i,j}(:,2), P{i,j}] = estimateDensity(blockData{i,j}, 100, [1, 1], 0);
%     %     end
%     if ~isempty(blockData{i,1})
%         for j = 1:blockDiv-1
%             [~, jaccardIndex(i,j)] = compareScatterOverlap(blockData{i,j}, blockData{i,j+1}, 200, 0,0);
%         end
%     end
% end

avgTrajectory = {};
for i = 1:size(fieldname,1)
    for j = 1:blockDiv
        if ~isempty(blockTrial.(fieldname(i)))
            [avgTrajectory{i,j}(:,1), avgTrajectory{i,j}(:,2)] = averageTrajectories(blockTrial.(fieldname(i))(:,j));
        else
            avgTrajectory{i,j} = [];
        end
    end
end

for i = 1:size(fieldname,1)
    for j = 1:blockDiv-1
        if ~isempty(avgTrajectory{i,j})
            Dist(i,j)=cdtw(avgTrajectory{i,j},avgTrajectory{i,j+1},0);
        else
            Dist(i,j) = 0;
        end
    end
end

function isNear = isNear(point1, point2, point3)
    % Calculate the threshold as the distance between point1 and point3
    threshold = sqrt((point1(1) - point3(1)).^2 + (point1(2) - point3(2)).^2)/2;
    
    % Calculate the distance between point1 and point2
    distance = sqrt((point1(:,1) - point2(:,1)).^2 + (point1(:,2) - point2(:,2)).^2);
    
    % Check if the distance is less than or equal to the threshold
    isNear = sum(distance <= threshold)/length(distance);
end
end