% PowerTestDebug
clear
clc
%%
% Automatically read the patient list and walk from the file
% names based on ..RWX_WalkY.mat
files = dir('*.mat');
fileNames = {files.name};
tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
if isempty(tokens) % Check if patient files were found
    folderPath = 'D:\CAPTURE Project\GUI\Data\';
    files = dir(fullfile(folderPath,'*.mat'));
    fileNames = {files.name};
    tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
    addpath(genpath(folderPath));
end
patientNumbers = tokens(:,1);
walkNumbers = tokens(:,2);
%%
Fs = 250;
winSize = 5;

%%
pDat = struct();
for fIdx = 15:22%1:length(tokens(:,1))
    pIdx = patientNumbers(fIdx);
    wIdx = walkNumbers(fIdx);
    dat = LoadNewData(pIdx,wIdx);
    eventNames = dat.evnts_tbl.Event;
    eventIdx = find(eventNames == "Lost End");

    pDat(pIdx).walk(wIdx).datF = dat.d_np(:,4);
    pDat(pIdx).t = dat.ntp_np;
    pDat(pIdx).eventTs = dat.evnts_tbl.NTP(eventIdx);
    datEpoch = zeros(winSize*2*Fs+1,length(pDat(pIdx).eventTs));
    for evIdx = 1:length(pDat(pIdx).eventTs)
        [~,npIdx] = min(abs(pDat(pIdx).eventTs(evIdx) - pDat(pIdx).t));
        startIdx = npIdx-round(Fs*winSize);
        endIdx = npIdx+round(Fs*winSize);
        datEpoch(:,evIdx) = pDat(pIdx).walk(wIdx).datF(startIdx:endIdx);
    end    
    pDat(pIdx).walk(wIdx).datE=datEpoch;
end


%% Power Calculations
dat = pDat(3).walk;

eeg = [];
for wIdx = 2:8
    d = dat(wIdx).datE;
    eeg = cat(2,eeg,d);
end

[pWelch,f] = pwelch(eeg,round(2*Fs),round(Fs),Fs);

[s,f,t] = spectrogram(eeg(:,1),round(2*Fs),round(Fs),0:.5:85,Fs,'yaxis');
waterplot(10*log10(s),f,t);

figure
spectrogram(eeg(:,1),round(Fs),round(Fs/2),0:.5:85,Fs,'yaxis')

figure
spectrogram(eeg(:,1),[],[],0:.5:85,Fs,'yaxis')


%%
% function PowerTest_Debug()
a=0;
app = startupFcn();
% evUnique = app.TransitionsLB.Items;
% uniqueDescrip = app.DescriptionLB.Items;
% orderStr = app.SortDropDown.Value;
% 
% app.RWdat.evTable.Include = ismember(app.RWdat.evTable.Event,evUnique);
% evTable = app.RWdat.evTable(app.RWdat.evTable.Include,:);
% [evTable.Description, descriptionUnique] = TextProcess.ExtractDescriptions(evTable.Description, '\[([^\]]+)\]');
% 
% for desIdx = 1:length(uniqueDescrip)
%     evTable.(uniqueDescrip{desIdx}) = contains(evTable.Description,uniqueDescrip{desIdx});
% end
% 
% tStart = app.RWdat.evTable.NTP(1);
% tArray = (evTable.NTP-tStart)/60;
% 
% % Find Event Time Series
% evTimeSeries = zeros(length(tArray),length(evUnique));
% for evIdx = 1:length(evUnique)
%     evName = string(evUnique{evIdx});
%     evTimeSeries(:,evIdx) = strcmp(evTable.Event,evName);
% end
% end
%%
function waterplot(s,f,t)
% Waterfall plot of spectrogram
    waterfall(f,t,abs(s)'.^2)
    set(gca,XDir="reverse",View=[30 50])
    xlabel("Frequency (Hz)")
    ylabel("Time (s)")
end
%%
function dat = LoadNewData(pID,wID)
    app.RWdat.pID = pID;
    app.RWdat.wID = wID;
    dat = load("RWNApp_RW"+app.RWdat.pID+"_Walk"+app.RWdat.wID+".mat");
end
%%
function app = startupFcn()
% Automatically read the patient list and walk from the file
% names based on ..RWX_WalkY.mat
files = dir('*.mat');
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
app.RWdat.patientNumbers = tokens(:,1);
app.RWdat.walkNumbers = tokens(:,2);

app.PatientTypeDD.Items = string(unique(app.RWdat.patientNumbers));
app.PatientTypeDD.ValueIndex = 1;
app.WalkTypeDD.Items = string(unique(app.RWdat.walkNumbers(app.RWdat.patientNumbers == app.PatientTypeDD.ValueIndex)));

app = LoadNewData(app,1,1);  
disp("Transitions")
disp(app.TransitionsLB.Items)

end




%%
function app = UpdateTransitionDescription(app)
evCounts = groupsummary(app.RWdat.evTable,"Event");
keepEvents = true(length(evCounts.Event),1);
if(app.RemoveBeginTransitionsCheckBox.Value)
    keepEvents = ~contains(evCounts.Event,"Beg","IgnoreCase",true) & keepEvents;
end
if(app.RemoveEndTransitionsCheckBox.Value)
    keepEvents = ~contains(evCounts.Event,"End","IgnoreCase",true) & keepEvents;
end
evUnique = evCounts.Event(evCounts.GroupCount>=app.ThresholdEditField.Value & keepEvents);

app.TransitionsLB.Items = evUnique;
app.RWdat.evTable.Include = ismember(app.RWdat.evTable.Event,evUnique);

[~, descriptionUnique] = TextProcess.ExtractDescriptions(app.RWdat.evTable.Description(app.RWdat.evTable.Include), '\[([^\]]+)\]');

app.DescriptionLB.Items = descriptionUnique;

app.SelectAllTransitionsCkb.Value = true;
app.SelectAllDescriptionCkb.Value = true;
app.TransitionsLB.Value = app.TransitionsLB.Items;
app.DescriptionLB.Value = app.DescriptionLB.Items;
end

function varNameList = FindDataVariables(app,varNameCode)
allVarStruct = app.RWdat.allVars;
fieldNames = fieldnames(allVarStruct);
fieldNames = fieldNames(contains(fieldNames,'d_'));
varNameList = [];
% Loop through all vars to find those starting with "d_"
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if startsWith(fieldName, varNameCode) % Check if the field name starts with "d_"
        newName = extractAfter(fieldName, varNameCode); % Remove "d_"
        value = allVarStruct.(fieldName); % Access the field's value
        % Check if the variable is empty
        % Currently only tables and arrays are counted
        if (~isempty(value))
            if(istable(value))
                for inerIdx=1:length(value.Properties.VariableNames)
                    varNameList = cat(1,varNameList,string(newName)+"."+string(value.Properties.VariableNames{inerIdx}));    
                end
            elseif(size(value,2)==1) 
                varNameList = cat(1,varNameList,string(newName));
            end
            
        end
    end
end
end

