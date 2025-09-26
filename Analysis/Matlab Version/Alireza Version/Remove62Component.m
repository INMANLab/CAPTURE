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
minSampNumToIgnore = 20; % in samples identify the length to not interpolate| 0.06s -> 15 Samples
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
for pIdx = 2:5%unique(pList.pID)'
    for wIdx = unique(pList.wID(pList.pID == pIdx))'
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
                nanIdx = isnan(channelDat);
                channelDat(marksIdx) = interp1(tDat(~(marksIdx|nanIdx)),channelDat(~(marksIdx|nanIdx)),tDat(marksIdx),'pchip');
                %---------------------|
                %--------------------- Find non-empty sections for Filtering
                nanIdx = isnan(channelDat);
                if(nanIdx(1)==1)
                    nanIdx(1)=0;
                end
                disp("NaN Ranges for RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+":")
                disp(string(round(diff(reshape(find(diff([nanIdx;0])),2,[]))/Fs,2)))
                % disp(diff(reshape(find(diff([0;nanIdx;0])),2,[]))'/Fs);
                epochIdx = reshape([0;find(diff([nanIdx;0]));length(nanIdx)],2,[]);
                epochIdx(1,:) = epochIdx(1,:)+1;
                epochIdx(:,diff(epochIdx)<minSampNumToInclude)=[];
                %---------------------|
                %----------- Graphics for check
                figure('units','normalized','outerposition',[0 0 .5 1])
                hold on
                plot(tDat,npDat(:,chIdx))
                plot(tDat(idxToIgnore),npDat(idxToIgnore,chIdx),'rx')
                gDat = channelDat;
                gDat(~marksIdx) = NaN;
                plot(tDat,gDat,"-+")
                Plot.AnnotateDataRange(epochIdx/Fs,-500,10)
                xline(epochIdx(1,:)/Fs,'g')
                xline(epochIdx(2,:)/Fs,'r')
                legend("raw data","ignore samples", "interpolated")
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
                %------------|
                %--------------------- Filtering
                dat.d_np_interp(:,chIdx) = channelDat;
                for chunk = 1:size(epochIdx,2)
                    % skipFlag = false;
                    % if(skipFlag)
                    %     continue;
                    % end


                    eegDat = channelDat(epochIdx(1,chunk):epochIdx(2,chunk));
                    %----------------- Remove 62.5 hz
                    params.Fs = Fs;
                    movingwin = [5, 2.5];
                    fRes = 1;
                    TW = movingwin(1)*fRes;
                    [e,v] = dpss(movingwin(1)*Fs,TW,10);
                    k = find(v>0.99,1,'last');
                    if(isempty(k))
                        k=1;
                    end
                    % figure
                    % lv = length(v);
                    % stem(1:lv,v,"filled")
                    % title("Proportion of Energy in [-w,w] of k-th Slepian Sequence")
                    params.tapers = [TW k];
                    params.pad = 1;
                    params.trialave = 0;
                    tau = 10;
                  
                    eegDatF= rmlinesmovingwinc(eegDat, movingwin, tau, params, 0,'',62.5,true);
                  

                    % figure
                    % eegDatF= rmlinesmovingwinc(eegDat, movingwin, tau, params, 0,'y',62.5,true);
                    % title("RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+" Chunk:"+chunk)


                    % params.fpass = [62,63];
                    % [S1,f]=mtspectrumc(eegDat,params);

                    % params.Fs = Fs;
                    % movingwin = [5, 2.5];
                    % params.fpass = [62,63];
                    % params.tapers = [3 5];
                    % params.pad = 1;
                    % [S1,f]=mtspectrumc(eegDatF,params);
                    % [datafit,Amps]=fitlinesc(eegDat,params,0,'y',62.5);
                    % eegDatF = eegDat;
                    % eegDatF = eegDatF-2*datafit;
                    % eegDatF= rmlinesc(eegDat, params, 0,'y',62.5);
                    % figure;plot(datafit)
                    
                    % %----------- Graphics for check
                    % figure('units','normalized','outerposition',[0 0 .5 1])
                    % % % hold on
                    % % % eegDat(isnan(eegDat))=0;
                    % % % rmlinesmovingwinc(eegDat, movingwin, tau, params, 0,'y',62.5,false,gain);
                    % % 
                    % pwelch(eegDat,round(movingwin(1)*250),round(movingwin(2)*250),round(movingwin(1)*250),Fs);hold on
                    % pwelch(eegDatF,round(movingwin(1)*250),round(movingwin(2)*250),round(movingwin(1)*250),Fs)
                    % %------------|

                    % %---------------- Notch Filter
                    % gain = gainAll(pIdx);
                    % wo = 62.5/(250/2);  
                    % bw = wo/65;
                    % [b62,a62] = iirnotch(wo, bw,gain); % 60Hz IIR filter
                    % eegDatF = filtfilt(b62,a62,eegDat);
                    % 
                    % %----------- Graphics for check
                    % figure
                    % pwelch(eegDat,round(movingwin(1)*250),round(movingwin(2)*250),round(movingwin(1)*250),Fs);hold on
                    % pwelch(eegDatF,round(movingwin(1)*250),round(movingwin(2)*250),round(movingwin(1)*250),Fs)
                    % %------------|
                    % %----------------------|


                    
                    % ------------- return value
                    channelDat(epochIdx(1,chunk):epochIdx(2,chunk)) = eegDatF;
                end
                dat.d_np_cleaned(:,chIdx) = channelDat;
                %---------------------|
            end
            %----------- Graphics for check
            for chIdx = 1:size(dat.d_np,2)
                figure('units','normalized','outerposition',[0 0 .5 1])

                gDat = dat.d_np(:,chIdx);
                subplot(4,1,1);
                plot((1:length(gDat))/Fs,gDat);
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat(isnan(gDat))=0;
                subplot(4,1,4)
                pwelch(gDat,movingwin(1)*250,movingwin(2)*250,[],Fs)
                title("Power Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat = dat.d_np_interp(:,chIdx);
                subplot(4,1,2)
                plot((1:length(gDat))/Fs,gDat);
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat(isnan(gDat))=0;
                subplot(4,1,4)
                hold on
                pwelch(gDat,movingwin(1)*250,movingwin(2)*250,[],Fs);
                title("Power Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat = dat.d_np_cleaned(:,chIdx);
                subplot(4,1,3)
                plot((1:length(gDat))/Fs,gDat);
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat(isnan(gDat))=0;
                subplot(4,1,4)
                hold on
                pwelch(gDat,movingwin(1)*250,movingwin(2)*250,[],Fs);
                title("Power Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
                legend("raw","interpolated","filtered")
            end
            %------------|
            % close all;
        else
            warning("Patient"+pIdx+"_Walk"+wIdx+": No np data found!")
        end
    end
end


%%

tEnd = 2;
t = 1/250:1/250:tEnd;
t = t*1000;
tN = 1/62.5:1/62.5:tEnd;
tN = tN *1000;
dat = A;

figure
subplot 211
plot(t,dat(1:length(t)))
xlabel("time (ms)")
hold on 
xline(tN)
subplot 212
pwelch(dat(1:length(t)),[],[],[],250)

wo = 62.5/(250/2);  
bw = wo/35;
[b62,a62] = iirnotch(wo, bw); % 60Hz IIR filter
datF = filtfilt(b62,a62,dat);

figure
subplot 211
plot(t,datF(1:length(t),1))
xlabel("time (ms)")
hold on 
xline(tN)
subplot 212
pwelch(datF(1:length(t)),[],[],[],250)

figure
plot(t,dat(1:length(t),1))
hold on
plot(t,datF(1:length(t),1))