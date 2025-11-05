%% Initialize
clear;
clc;
close all;

%% Load the data
dataPath = "D:\CAPTURE Project\GUI\Data\";
kdePath = "D:\CAPTURE Project\GUI\EventKDE\BW2_20Raters\KDEs_New\";
files = dir(dataPath+"*.mat");
fileNames = {files.name};
tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
if isempty(tokens) % Check if patient files were found
    disp('No patient file has been found');
    folderPath = uigetdir(pwd, 'Select the folder that contains individual data');
    if folderPath ~= 0 % If the user selects a folder
        files = dir(fullfile(folderPath,'*.mat'));
        fileNames = {files.name};
        tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
        addpath(genpath(folderPath));
    else
        error('Patient files not found');
    end
end
patientNumbers = tokens(:,1)';
walkNumbers = tokens(:,2)';

%% Load Event Data
eventTable = [];
for pwIdx = 1:length(patientNumbers)
    pIdx = patientNumbers(pwIdx);
    wIdx = walkNumbers(pwIdx);
    % try
        dat = load(dataPath+"RWNApp_RW"+pIdx+"_Walk"+wIdx+".mat");
    % catch
    %     continue;
    % end
    kdeDat = readtable(kdePath+"RWN"+pIdx+"_"+wIdx+"_KDE.csv");
    events = dat.evnts_tbl;
    kdeVal = kdeDat.Y;
    kdeTime = kdeDat.X;
    events.time = events.NTP - dat.ntp_kde2(1);
    events.kdeScore = zeros(size(events,1),1);
    events.kdeDistance = zeros(size(events,1),1);
    events = events(events.time>0,:);
    for evIdx = 1:size(events,1)
        tdiff = kdeTime - events.time(evIdx);
        [timeDist, kdeIndex] = min(abs(tdiff));
        events.kdeScore(evIdx) = kdeVal(kdeIndex);
        events.kdeDistance(evIdx) = timeDist;
    end
    events.pId = repmat(pIdx,size(events,1),1);
    events.wId = repmat(wIdx,size(events,1),1);
    eventTable = cat(1,eventTable,events);
        % evCounts = groupsummary(app.RWdat.evTable,"Event");

end
%%
temp = eventTable(eventTable.kdeDistance>.1,:);
eventTable.kdeDistance(eventTable.kdeDistance<=.1) = -1;
eventTable.kdeDistance(eventTable.kdeDistance>.1) = 1./eventTable.kdeDistance(eventTable.kdeDistance>.1); 
eventTable.kdeDistance(eventTable.kdeDistance==-1)=1;
eventTable.Score = eventTable.kdeScore.*eventTable.kdeDistance;
%% 
T = eventTable;
figure;
histogram(T.Score,100)

lowThresh  = .1;prctile(T.Score, 25);
highThresh = .25;prctile(T.Score, 75);

% Classify each score
T.Level = repmat("Mid", height(T), 1);
T.Level(T.Score <= lowThresh)  = "Low";
T.Level(T.Score >= highThresh) = "High";
T.Level = categorical(T.Level, ["Low","Mid","High"]);

levelCounts = countcats(T.Level);
labels = categories(T.Level);

figure;
donutchart(levelCounts, labels);
title("Overall Score Level Distribution Lowthr: "+lowThresh + " Highthr: "+ highThresh);

%% 
T = eventTable;
disp(unique(T.Event))
T = T(ismember(T.Event,["Doorway","Lost Beg","Lost End","Outdoor Beg","Outdoor End","Stop Beg","Stop End"]),:);
figure;
histogram(T.Score,100)

lowThresh  = prctile(T.Score, 25);
highThresh = prctile(T.Score, 75);

% Classify each score
T.Level = repmat("Mid", height(T), 1);
T.Level(T.Score <= lowThresh)  = "Low";
T.Level(T.Score >= highThresh) = "High";
T.Level = categorical(T.Level, ["Low","Mid","High"]);

levelCounts = countcats(T.Level);
labels = categories(T.Level);

figure;
donutchart(levelCounts, labels);
title("Overall Score Level Distribution Lowthr: "+lowThresh + " Highthr: "+ highThresh);

T.Category = T.Event;
cats = unique(T.Category);

figure;
for i = 1:numel(cats)
    thisGroup = T(T.Category == cats(i), :);
    c = countcats(thisGroup.Level);
    subplot(2, ceil(numel(cats)/2), i);
    donutchart(c, categories(thisGroup.Level));
    title(string(cats(i)));
end
%% No Zero
T = eventTable(eventTable.Score>0.01,:);
figure;
histogram(T.Score,100)

lowThresh  = .1;prctile(T.Score, 25);
highThresh = .25;prctile(T.Score, 75);

% Classify each score
T.Level = repmat("Mid", height(T), 1);
T.Level(T.Score <= lowThresh)  = "Low";
T.Level(T.Score >= highThresh) = "High";
T.Level = categorical(T.Level, ["Low","Mid","High"]);

levelCounts = countcats(T.Level);
labels = categories(T.Level);

figure;
donutchart(levelCounts, labels);
title("Overall Score Level Distribution Lowthr: "+lowThresh + " Highthr: "+ highThresh);

%% 
T = eventTable;
disp(unique(T.Event))
T = T(ismember(T.Event,["Doorway","Lost Beg","Lost End","Outdoor Beg","Outdoor End","Stop Beg","Stop End"]),:);
figure;
histogram(T.Score,100)

lowThresh  = .05;prctile(T.Score, 25);
highThresh = .25;prctile(T.Score, 75);

% Classify each score
T.Level = repmat("Mid", height(T), 1);
T.Level(T.Score <= lowThresh)  = "Low";
T.Level(T.Score >= highThresh) = "High";
T.Level = categorical(T.Level, ["Low","Mid","High"]);

levelCounts = countcats(T.Level);
labels = categories(T.Level);

figure;
donutchart(levelCounts, labels);
title("Overall Score Level Distribution Lowthr: "+lowThresh + " Highthr: "+ highThresh);

T.Category = T.Event;
cats = unique(T.Category);

figure;
for i = 1:numel(cats)
    thisGroup = T(T.Category == cats(i), :);
    c = countcats(thisGroup.Level);
    subplot(2, ceil(numel(cats)/2), i);
    donutchart(c, categories(thisGroup.Level));
    title(string(cats(i)));
end
%%
% T.Group = categorical(repmat("Events", height(T), 1));
% % Choose the column you want to group by (e.g., Group)
% groupVar = 'Group';
% 
% % Count number of High/Mid/Low per group
% counts = groupsummary(T, groupVar, "count", "Level");
% 
% % (Optional) If you want to aggregate all groups together:
% % levelCounts = countcats(categorical(T.Level));
% 
% % For pie chart of Level distribution across all groups
% levelCounts = groupcounts(T.Level);
% labels = categories(categorical(T.Level));
% 
% figure;
% pie(levelCounts, labels);
% title('Distribution of Score Levels');