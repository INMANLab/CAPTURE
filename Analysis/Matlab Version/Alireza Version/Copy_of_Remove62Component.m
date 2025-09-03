clear;
clc;
%% add Chronux to the path
ChronuX_path = "D:\Toolboxes\chronux_2_12";
addpath(genpath(ChronuX_path));
%% find patient .mat files and add the directory to path
% requires +TextProcess Folder
folderPath = uigetdir(pwd, 'Select the folder that contains individual data');
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
%% Reading the data
for pIdx = unique(pList.pID)'
    for wIdx = unique(pList.wID(pList.pID == pIdx))'
        dat = load("RWNApp_RW"+pIdx+"_Walk"+wIdx+".mat");
        if(isfield(dat, "d_np"))
            npDat = dat.d_np;
            Fs = dat.fs_np;            
            for chIdx = 1:size(npDat,2)
                %---------- Find empty sections
                % tempDat = npDat(:,chIdx);
                % tempDat(isnan(tempDat))=-512;
                % tempDat(tempDat>-500) = 0;
                % missIdx = [1;find(diff(tempDat));length(tempDat)];

                missIdx = isnan(npDat(:,chIdx))|npDat(:,chIdx)<=-500;

                % missIdx = round(rand(50,1));
                % 
                % figure
                % plot(missIdx,'*');
                % A = [diff(missIdx);0];
                % hold on
                % % plot(A)
                % A(A>0)=0;
                % % stem(-A)
                % A = cumsum(missIdx).*(1-missIdx);
                % stem(A)
                % A = [diff(cummax(A));0];
                % stem(A)
                % B = B(end:-1:1);
                % B = movprod(B,[0 1]);
                % stem(B)

                %----------- Temp test Graph
                figure
                plot(missIdx,'*');
                hold on
                epochIdx = reshape(find(diff([0;missIdx;0])),2,[]);
                % epochIdx(2,:) = epochIdx(2,:)-1;
                epochSize = diff(epochIdx);
                epochIdx = epochIdx(:,epochSize>5);

                % missIdx = [1;find(diff(missIdx));length(tempDat)];
                % %----------- Temp test Graph
                % figure
                % plot(npDat(:,chIdx)/std(npDat(:,chIdx),"omitnan"))
                % hold on
                % plot(missIdx(1:end-1),diff(missIdx)/std(diff(missIdx),"omitnan"))
                % %---------------------------|
                % % missIdx = diff(missIdx);
                % epochIdx = reshape(missIdx,2,[]);
                % epochIdx(:,diff(epochIdx)<10*Fs)=[];
                % xline(epochIdx(1,:),'r')
                % xline(epochIdx(2,:),'g')


                

                epochIdx = reshape([0;find(missIdx);length(missIdx)],2,[]);
                epochIdx(1,:) = epochIdx(1,:)+1;
                for chunk = 1:size(epochIdx,2)
                    eegDat = npDat(epochIdx(1,chunk):epochIdx(2,chunk),chIdx);
                    %---- Interpolate -500 markers
                    marksIdx = eegDat<=-500;
                    t = (epochIdx(1,chunk):epochIdx(2,chunk))/Fs;
                    eegDat(marksIdx) = interp1(t(~marksIdx),eegDat(~marksIdx),t(marksIdx),'pchip');
                    %----------------- Remove 62.5 hz
                    params.Fs = Fs;
                    params.tapers = [1 1];
                    params.pad = 1;
                    params.trialave = 0;

                    tau = 10;
                    movingwin = [3, 1.5];

                    eegDatF = rmlinesmovingwinc(eegDat, movingwin, tau, params, [],[],62.5,true);
                    %------------- return value
                    npDat(epochIdx(1,chunk):epochIdx(2,chunk),chIdx) = eegDat;
                end
            end
            dat.d_np_cleaned = npDat;
            % Graphics for check
            gIdx = 1;
            for chIdx = 1:size(dat.d_np,2)
                figure

                gDat = dat.d_np;
                subplot(size(gDat,2),2,gIdx);gIdx = gIdx+1;
                plot((1:size(gDat,1))/Fs,gDat(:,chIdx));
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
                
                subplot(size(gDat,2),2,gIdx);gIdx = gIdx+1;
                pwelch(gDat(:,chIdx),[],[],[],Fs);
                title("Power Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")

                gDat = dat.d_np_cleaned;
                subplot(size(gDat,2),2,gIdx);gIdx = gIdx+1;
                plot((1:size(gDat,1))/Fs,gDat(:,chIdx));
                xlabel("time (sec)")
                ylabel("Amplitude")
                title("Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
                
                subplot(size(gDat,2),2,gIdx);gIdx = gIdx+1;
                pwelch(gDat(:,chIdx),[],[],[],Fs);
                title("Power Raw: RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+"")
            end





        else
            warning("Patient"+pIdx+"_Walk"+wIdx+": No np data found!")
        end
    end
end