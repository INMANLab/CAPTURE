clear;
clc;
close all;
%% add Chronux to the path
ChronuX_path = "D:\Toolboxes\chronux_2_12";
addpath(genpath(ChronuX_path));
set(0, 'DefaultTextInterpreter', 'none')
%% find patient .mat files and add the directory to path
% requires +TextProcess Folder
folderPath = 'D:\CAPTURE Project\GUI\Data';%uigetdir(pwd, 'Select the folder that contains individual data');

if folderPath ~= 0 % If the user selects a folder
    files = dir(fullfile(folderPath,'*.mat'));
    fileNames = {files.name};
    tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
    pList = table;
    pList.pID = tokens(:,1);
    pList.wID = tokens(:,2);
    addpath(genpath(folderPath));
else
    error('Patient files not found');
end
%% Parameters
minSampNumToIgnore = 3; % in samples identify the length to not interpolate| 0.06s -> 15 Samples
minSampNumToInclude = 10*250; % Minimum number of samples to process the chunk
% gainAll = [0,.1];
% patientPars =table;
%%
% winSize = movingwin(1);
% nOverlap = movingwin(2);
% skipFlag = true;
% % skipFlag = false;
% 
% tempTab = table(pIdx,wIdx,chIdx,chunk,winSize,nOverlap,TW,k,skipFlag);
% patientPars = cat(1,patientPars,tempTab);
%% Reading the data
for pIdx = 3:3 %unique(pList.pID)'
    for wIdx = 1:4%unique(pList.wID(pList.pID == pIdx))'
        disp("=========================================================")
        disp("=========================================================")
        disp("=========================================================")
        disp("processing data for RW"+pIdx+"_Walk"+wIdx)
        dat = load("RWNApp_RW"+pIdx+"_Walk"+wIdx+".mat");
        if(isfield(dat, "d_np"))
            npDat = dat.d_np;
            dat.d_np_interp = zeros(size(npDat));
            dat.d_np_cleaned = zeros(size(npDat));
            tDat = dat.ntp_np;
            tDat = tDat-tDat(1);
            Fs = dat.fs_np;            
            for chIdx = 1:size(npDat,2)
                channelDat = npDat(:,chIdx);
                %--------------------- Find long empty areas and replace with NaN
                missIdx = channelDat<-500;
                missIdx = [0;missIdx;0]; % Pad with Zeros
                missLabel = diff([0;missIdx]);
                missLabel(missLabel<0)=0;
                missLabel = cumsum(missLabel).*missIdx;
                sampNum = cumsum(missIdx).*(1-missIdx);
                sampNum = diff([cummax(sampNum);max(sampNum)]);
                idxToIgnore = ismember(missLabel,find(sampNum(sampNum>0)>=minSampNumToIgnore));
                % %----------- Graphics for check
                % figure
                % plot(missIdx,'*');
                % hold on
                % stem(missLabel)
                % stem(sampNum)
                % stem(idxToIgnore)
                % %------------|
                idxToIgnore = idxToIgnore(2:end-1);
                channelDat(idxToIgnore) = NaN;
                %---------------------|
                %--------------------- Interpolate the remaining -500 markers
                marksIdx = channelDat<=-500;

                movingwin = [11, 3];
                figure
                subplot(2,1,1)
                plot(tDat,marksIdx)
                title("RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
                subplot(2,1,2)
                pwelch(double(marksIdx),round(movingwin(1)*250),round(movingwin(2)*250),round(movingwin(1)*250),Fs);
                hold on
                xline(62.5)
                a=0;
            end
        else
            warning("Patient"+pIdx+"_Walk"+wIdx+": No np data found!")
        end
    end
end