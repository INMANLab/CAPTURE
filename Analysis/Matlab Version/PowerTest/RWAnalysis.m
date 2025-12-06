classdef RWAnalysis < handle
    % Specgram Steps (best to use gui):
    % 1) RWA = RWAnalysis;
    % 2) <optional> RWA.loadData('PatientIdx',1,'transtype','Outdoor Beg'); %plotMultTransSpecGramPerm loads data automatically
    % 3) RWA.plotMultTransSpecGramPerm('transtype','Doorway','desctype','in2in','regiontype','AntHipp','walktype','All Walks','permtype','standard','correctiontype','cluster'); 
    %
    % InOut analysis steps:
    % 1) RWA = RWAnalysis;
    % 2) RWA.getMultSegData;
    % 3) RWA.plotMultInOutDiff;
    % 4) RWA.plotMultInOutGLME;
    %
    % Notes:
    % 1) Redo in/out boxplot to incorporate lensky and change marker shape by brain region (Done)
    % 2) Baseline indoor/indoor transitions (Done)
    % 3) Combine hippocampus into one (hipp vs all others) (Done - added | symbol for regions)
    % 4) Need to mark stopping at a closed door vs walking through open door for all doors (Done - keywords added)
    % 6) Incorporate luminance, fixation rate, audio features, cadence
    % 7) eBosc pepisode (2-20Hz)
    % 
    % 8) Add ability to run multiple "Transitions" in the same spectrogram (Correct Turn Beg|Incorrect Turn Beg) (Done)
    % 9) Plot band filtered individual trials color coded by participant (Done)
    % 10) Incorporate the Pupil Eye tracking data in to GLMM model (Done)
    % 11) Incorporate the Real World Nav Event Segmentation into the GLMM Model (Done)
    % 12) Make tool to create event times for RWN Event Segmentation to do event based spectrograms at different levels of participant agreement (>15 people hit the spacebar at these time points, <5 people hit the spacebar at these time points, etc.)
    % 13) Remove outliers before filtering (Done)
    % 14) Adjustable smoothing window for PermSpecGram option in GUI (Done)
    % 15) Save figure checkbox (Done)
    % 16) Separate figure for trials data (Done)
    % 17) Print threshold for clusters**** (Does not make sense)
    % 18) Enable TransitionDD/RegionDD for difference in GUI (Done)
    % 19) Trials plot for other permutation functions (Done)
    % 20) Make new permutation code to use baseline and FDR (Done)
    % 21) Make variable window size (Done)
    % 22) Use sum(cluster values) in the correction instead of just size (Done)
    % 23) Increase freq range for new permutation code (Done)
    % 24) Variable window size for GLM
    % 25) Event related GLM
    % 26) Leave one for Patient drop down in GUI (Done)
    % 27) Boxplot comparison (Done)
    %
    % 1) Add KDE to GLM (Done)
    % 2) Spectrograms around KDE peaks (high, mid, low, trough) (Done)
    % 3) Event types related to KDE data broken down into pie charts
    % 4) Correlation confusion matrix of all bands and all predictors (by channel i.e. 20 for all walks)
    % 5) GLM outputs to text (done)
    % 6) Walk spectrograms with all predictors (done)
    % 7) Boxplot for each of the predictors/freq(pwr) by indoor/outdoor (done)
    % 8) The updated Kiersten data has new events that change the glme results slightly (not sure how...)
    %
    % 1/f tilt calculation across time (or before/after) for each transition specgram (fooof v2.0 donaghue) in gui (done)
    % *Correlation confusion matrix - run stats and combine by region, patient, etc (new gui?)
    % *Add head movements to event alignment (eye tracking IMU?) - add an event for this
    % *Add skin conductance data
    % *Correct/Incorrect turns gamma timing with low freq suppression (phase amp coupling) (GUI?)
    % *Plot average KDE value for each transition (break down doorways - indoor/indoor, indoor/outdoor, outdoor/indoor, outdoor/outdoor, closed, open, etc. -> RWN_annotationbrackets.xls)
    % *Exclude overlap trials in gui
    % *Comparing KDE to expert rater (how much lag?)
    % *eBosc pepisode (2-20Hz, 2-5, 5-8, 8-12) and alignment to events
    % *Leave one predictor out for glme and record change in r2 (backwards step regression)
    % *GLME filtered by region (new gui?)
    % *Add another velocity condition to completely separate change from no change
    %
    % Every KDE peak as an option (done)
    % KDE peaks that don't have a speed change (velocity dropdown accomplishes this) (done)
    % Ratio between before/after in specgram power for gui (add a average value to title) (done)
    % Only pick one electrode per participant in AntHipp to prevent double dipping (add custom dropdown and new window with list of all chans from all patients) (done)
    % Checkbox for including velocity (done)
    % Iterate and save for each category (done)
    % Add tercile for velocity (add to title) (done)
    % Add check for overlapping transitions and include value in figure (percentage of trials and average time of overlap) (done)
    % Add more detail to figure titles for saving (done)
    


    %%%%%%%%%%%%%%%%%%%
    properties
        RootDir; PatientID; DB; WalkNames; WalkNums; DParsed; PatientList; PatientIdx;
        ChanLabels; InOutTimes; DTable; MultInOut; DisabledWalks; FB; Vid; RegionTable;
        MultTrans; StopGoWalks; AnalysisFile; WinSegSec; WinTransSec;
    end

    %%%%%%%%%%%%%%%%%%%
    methods %initialization, calculations, etc

        function obj = RWAnalysis(varargin) %constructor (runs when object is created for the first time)
            obj.parseInputs(varargin{:});

            obj.RootDir = '\\155.100.91.44\D\Tyler\RealWorld';
            obj.AnalysisFile = fullfile(obj.RootDir,'RWAnalysis.mat');
            obj.WinSegSec = 2; %Also smoothing kernel for specgrams in trans analysis
            obj.WinTransSec = 12;

            obj.PatientList = {
                'RW1\Original' %1(some alignment problems, xsens is good), 2(lots of dropouts w/alignment problems, xsens missing), 3(a few alignment problems, xsens missing), 4(good), 5(good), 6(good, small xsens dropout), 7(reverse, dropouts in np and xsens, alignment good) 
                'RW2\Original' %1, 2(xsens drops out), 3, 4, 5(marks/led do not align well), 6(no np marks), 7(no np marks)
                'RW3\Original' %d1(needs events), 2, d3(pupil corrupt - needs events), 4, 5(xsens not synced), 6(xsens missing timestamps), 7, 8(1st xsens not synced, 2nd xsens is corrupt)
                'RW4\Original' %1, 2, d3 (needs events), 4, 5, 6, 7, 8
                'RW5\Original' %1, 2, 3, 4, 5, d6(rp file missing), d7(rp file missing), d8(rp file missing)
            };

            obj.DisabledWalks = {
                []
                [] %[5,6,7]
                [] %[1,3]
                [] %3
                [] %[6,7,8]
            };

            obj.StopGoWalks = {
                [6,5]
                [5,6]
                [7,6]
                [6,7]
                [7,6]
            };

            obj.ChanLabels = [
                {'Left Anterior Hippocampus (Entorhinal Cortex/Dentate Gyrus)','Left Anterior Hippocampus (Lateral-CA1/Perirhinal)','Right Anterior Hippocampus (Ventral - CA1/Perirhinal)','Right Ventrolateral Temporal Cortex'}
                {'Left Entorhinal - Perirhinal','Left Ventrolateral Temporal Cortex','Right Anterior Hippocampus (Ventral - Perirhinal/Subiculum)','Right Ventroateral Temporal Cortex (Anterior Fusiform Gyrus)'}
                {'Left Anterior Hippocampus (Ventral - Entorhinal/Subiculum)','Left Anterior Hippocampus (Ventral - Subiculum/Perirhinal)','Left Amygdala - Anterior Hippocampus (Long Axis)','Left Anterior Hippocampus - Mid Hippocampus (Long Axis)'}
                {'Left Amygdala - Anterior Hippocampus (Long Axis - Dentate Gyrus)','Left Mid Hippocampus - Posterior Hippocampus','Right Amygdala - Anterior Hippocampus (CA1)','Right Mid Hippocampus - Posterior Hippocampus (Long Axis - parahippocampal border)'}
                {'Left Amygdala - Anterior Hippocampus (Long Axis - Dentate Gyrus)','Left Mid Hippocampus - Posterior Hippocampus (Long Axis)','Right Amygdala - Anterior Hippocampus (Long Axis - subiculum)','Right Mid Hippocampus - Posterior Hippocampus (Long Axis- parahippocampal border)'}
            ];

            %Finding regions
            obj.RegionTable = [];
            for m=1:size(obj.ChanLabels,1)
                for k=1:4 %chan
                    chanlabel = obj.ChanLabels{m,k};
                    b1 = contains(chanlabel,("Amygdala"|"Anterior Hippocampus"));
                    b2 = contains(chanlabel,("Temporal"));
                    b3 = contains(chanlabel,("Entorhinal"|"Perirhinal"));
                    if b1
                        llb = 1; %anterior
                        llbname = 'AntHipp';
                    else %posterior
                        if b2 %temporal
                            llb = 2;
                            llbname = 'LatTemp';
                        elseif b3 %entorhinal
                            llb = 3;
                            llbname = 'Ent+Peri';
                        else
                            llb = 4;
                            llbname = 'PostHipp+Para';
                        end
                    end
                    obj.RegionTable = vertcat(obj.RegionTable,table(m,k,{chanlabel},llb,{llbname},'VariableNames',{'Patient','Chan','ChanLabel','Region','RegionLabel'}));
                end
            end
            
        end

        function p = parseInputs(obj,varargin)
            %Parses name-value input pairs and saves the value to the
            %matching class property or exports to p structure.
            if rem(nargin-1,2)
                error('Name-value pairs needed as input!')
            end

            PropNames = properties(obj); 
            InputNames = varargin(1:2:end);
            InputVals = varargin(2:2:end);       

            p = [];
            p.plottype = ''; %?
            p.transtype = ''; %transition type from events table (i.e. 'Outdoor Beg')
            p.regiontype = ''; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            p.customregion = []; %nx2 where cols are patient,chan (i.e. [1,1;1,2;2,1;2,4], used in filterMultTransData)
            p.patienttype = ''; %'All Patients', '1,2,5'
            p.walktype = ''; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '2,4,6'
            p.desctype = ''; %bracketed keyword found in the description (i.e. closed, open, in2in, in2out, etc.)
            p.patient = []; %1,2,3,etc
            p.walknum = []; %walk number (1,2,3,etc)
            p.chan = []; %1,2,3,4
            p.evntnum = []; %1,2,3,etc
            p.fpass = []; 
            p.fidx = []; %index into list of freq bins (plotMultInOutDiff)
            p.saveflag = false;
            p.permtype = ''; %standard, zscore (use with plotMultTransSpecGramPerm)
            p.correctiontype = ''; %cluster, pixel, fdr (use with plotMultTransSpecGramPerm)
            p.pval = 0.05; %pval for specgrams
            p.pvalclust = 0.01; %pval to use for cluster correction
            p.transrng = [-10,10]; %transition window in sec 
            p.normrng = [-10,10]; %normalization window in sec (must be within transrng)
            p.fullwalknorm = false; %each chan is normalized across the entire walk (bypasses normrng above) 
            p.smoothwin = 4; %4sec smoothing window across time for PermSpecGram
            p.clim = []; %typically [-1,1] for dB and [-10,10] for zscore
            p.plottrials = false; %plot individual trials for a specgram for the specified trialsfreqrng
            p.trialsfreqrng = []; %avg freq band in Hz for plotting specgram trials
            p.plotboxcomp = false; %boxplot comparison of mean power across trials for two rectangular regions in a spectrogram
            p.boxcomprng = []; %rows are box1/box2, cols are time/freq limits for each box
            p.predtype = ''; %predictor type -> Vel, Fix, KDE, nAmb (plotMultInOutPred)
            p.veltype = ''; %VelHigh, VelLow, VelHighTercile, VelMidTercile, VelLowTercile
            p.includevel = false; %include velocity overlay in MultTransSpecGram plots

            if ~all(ismember(InputNames,[PropNames;fieldnames(p)]))
                error('Input names are incorrect!')
            end

            for k=1:length(InputNames)
                if any(strcmp(InputNames{k},PropNames))
                    if isa(InputVals{k},'string')
                        obj.(InputNames{k}) = char(InputVals{k});
                    else
                        obj.(InputNames{k}) = InputVals{k};
                    end
                else
                    p.(InputNames{k}) = InputVals{k};
                end
            end
        end %parseInputs

        function loadData(obj,varargin)
            %loadData('PatientIdx',1);
            obj.PatientID = [];
            obj.DB = [];
            obj.DTable = {};
            obj.WalkNames = [];
            obj.WalkNums = [];
            obj.DParsed = [];
            obj.PatientIdx = [];
            obj.InOutTimes = [];
            obj.Vid = [];

            obj.parseInputs(varargin{:});
            if isempty(obj.PatientIdx)
                error('Need to specify PatientIdx!');
            end

            obj.PatientID = regexp(obj.PatientList{obj.PatientIdx},'(RW)\d{1}','match','once');
            obj.WalkNames = dir(obj.RootDir);
            obj.WalkNames = obj.WalkNames(contains({obj.WalkNames.name},obj.PatientID));
            obj.WalkNames = regexp({obj.WalkNames.name},'Walk\d+','match','once');
            obj.WalkNames(cellfun(@isempty,obj.WalkNames)) = [];
            obj.WalkNums = cellfun(@str2double,regexp(obj.WalkNames,'\d+','match','once'));
            TotalWalks = length(obj.WalkNames);

            %[notchB,notchA] = iirnotch(60/125,0.012);

            obj.DTable = cell(TotalWalks,17);
            for k=1:TotalWalks
                fprintf('Loading patient %s walk %d...\n',obj.PatientID,k);

                WalkNum = obj.WalkNums(k);
                if ~ismember(WalkNum,obj.DisabledWalks{obj.PatientIdx})
                    % Fname = fullfile(obj.RootDir,sprintf('Walk%0.0f\\RWNApp_%s_Walk%0.0f.mat',WalkNum,obj.PatientID,WalkNum));
                    Fname = fullfile(obj.RootDir,sprintf('RWNApp_%s_Walk%0.0f.mat',obj.PatientID,WalkNum));
                    if isfile(Fname)
                        SS = load(Fname);

                        %Check data
                        if WalkNum~=k
                            disp('WalkNum is not monotonically increasing!')
                        end

                        if isempty(SS.evnts_tbl)
                            error('Events table is empty for walk %0.0f',WalkNum);
                        end
                       
                        if isempty(SS.d_xs)
                            fprintf('Xsens data is missing for walk %0.0f\n',WalkNum);
                        end

                        if ~isempty(SS.ntp_np)
                            tsamp = floor(SS.ntp_np.*SS.fs_np); tsamp = tsamp - tsamp(1) + 1; %these should be monotonically increasing samples
                            tsec = ((0:size(SS.d_np,1)-1)./SS.fs_np+SS.ntp_np(1))'; %time in sec comparable to ntp_np but using monotonically increasing samples
                            diff_tsamp = diff(tsamp);
                            if any(abs(SS.ntp_np-tsec)>(2/SS.fs_np))
                                fprintf(['Neuropace NTP time vector for %s has some jitter.\nMin/max jitter from expected value of 1 is %0.0f/%0.0f samples.' ...
                                    '\nTotal # jitter values %0.0f/%0.0f.\nMax deviation from monotonic is %0.2e sec.\n'],obj.PatientID,min(diff_tsamp),...
                                    max(diff_tsamp),sum(diff_tsamp~=1),length(diff_tsamp),max(abs(SS.ntp_np-tsec)));
                            end
                        end

                        if any(diff(SS.ntp_np)<=0)
                            disp('Time vector is not monotonically increasing!')
                        end

                        %Verify table values
                        for m=1:size(SS.evnts_tbl,1)
                            [~,idx] = min(abs(SS.ntp_np - SS.evnts_tbl.NTP(m)));
                            if idx~=SS.evnts_tbl.NPSample(m)
                                fprintf('Mismatch between NPSample and NTP times in event table row %d walk %d!\n',m,k);
                            end
                            [~,idx] = min(abs(SS.ntp_pupil - SS.evnts_tbl.NTP(m)));
                            if idx~=SS.evnts_tbl.PupilFrame(m)
                                fprintf('Mismatch between PupilFrame and NTP times in event table row %d walk %d!\n',m,k);
                            end
                            [~,idx] = min(abs(SS.ntp_gp - SS.evnts_tbl.NTP(m)));
                            if idx~=SS.evnts_tbl.GoProFrame(m)
                                fprintf('Mismatch between GoProFrame and NTP times in event table row %d walk %d!\n',m,k);
                            end
                        end
                        
                        %Checking for duplicates in table
                        [uE,~,uEIdx] = unique(SS.evnts_tbl.Event);
                        ridx = [];
                        for m=1:size(uE,1)
                            idx = find(m==uEIdx);
                            if length(idx)>1
                                pf = SS.evnts_tbl.PupilFrame(idx);
                                [uP,~,uPIdx] = unique(pf);
                                if length(uP)<length(pf)
                                    fprintf('Duplicates were found in the event table for walk %d, event "%s"! Removing...\n',k,uE(m));
                                end
                                for n=1:length(uP)
                                    idx2 = (uPIdx==n);
                                    if sum(idx2)>1
                                        idx3 = idx(idx2); %indices in evnts_tbl
                                        SS.evnts_tbl.Description(idx3(1)) = cell2mat(SS.evnts_tbl.Description(idx3)'); %combining descriptions from all duplicates
                                        ridx = [ridx;idx3(2:length(idx3))]; %saving a list of indices to remove (keeping first)
                                    end
                                end
                            end
                        end
                        SS.evnts_tbl(ridx,:) = [];
                        disp('Removing rows...')
                        disp(ridx')

                        %Checking for duplicate Walk Beg and Walk End
                        if sum(contains(SS.evnts_tbl.Event,'Walk Beg'))>1
                            fprintf('Duplicate "Walk Beg" events were found for %s walk %d\n',obj.PatientID,k)
                        end
                        if sum(contains(SS.evnts_tbl.Event,'Walk End'))>1
                            fprintf('Duplicate "Walk End" events were found for %s walk %d\n',obj.PatientID,k)
                        end

                        %SS.d_np = filtfilt(notchB,notchA,SS.d_np);

                        %Interpolate marks
                        warning('off','MATLAB:interp1:NaNstrip');
                        ts = obj.findMarks(SS.d_np);
                        if ~isempty(ts)
                            TS = reshape((ts+(-1:14))',1,[]);
                            t = setdiff(1:size(SS.d_np,1),TS);
                            d_sub = interp1(t,SS.d_np(t,:),TS,'spline','extrap');
                            SS.d_np(TS,:) = d_sub;
                        end

                        %Set outliers to nan
                        SS.d_np(SS.d_np<-500) = nan;

                        obj.DTable(k,1) = {SS.d_np}; %data only for 1 walk
                        obj.DTable(k,2) = {SS.ntp_np}; %time in sec (this does not start at zero)
                        obj.DTable(k,3) = {SS.fs_np};

                        obj.DTable(k,4) = {SS.d_xs};
                        obj.DTable(k,5) = {SS.ntp_xs};
                        obj.DTable(k,6) = {SS.fs_xs};

                        obj.DTable(k,7) = {SS.ntp_pupil};
                        obj.DTable(k,8) = {SS.ntp_gp};

                        if isfield(SS,'ntp_gaze')
                            obj.DTable(k,9) = {SS.d_gaze_fix};
                            obj.DTable(k,10) = {SS.ntp_gaze};
                            obj.DTable(k,11) = {SS.fs_gaze};
                        end

                        if isfield(SS,'ntp_amb')
                            obj.DTable(k,12) = {SS.d_amb};
                            obj.DTable(k,13) = {SS.ntp_amb};
                            obj.DTable(k,14) = {SS.fs_amb};
                        end

                        if isfield(SS,'ntp_kde')
                            %KDEPeakHighTercile, KDEPeakMidTercile, KDEPeakLowTercile,
                            %KDEPeakHigh, KDEPeakLow, KDETrough, KDEPeakAll
                            %(trough is defined as trs<prctile(pks,5))
                            obj.DTable(k,15) = {SS.d_kde};
                            obj.DTable(k,16) = {SS.ntp_kde};
                            obj.DTable(k,17) = {SS.fs_kde};
                            [pks,ploc] = findpeaks(SS.d_kde,'MinPeakProminence',0.5,'NPeaks',50); %peaks
                            [~,pidx] = sort(pks,'descend'); pks = pks(pidx); ploc = ploc(pidx); pks5 = prctile(pks,5);
                            [trs,tloc] = findpeaks(-SS.d_kde,'MinPeakProminence',0.5,'NPeaks',50); trs = -trs; %troughs
                            [~,tidx] = sort(trs,'descend'); trs = trs(tidx); tloc = tloc(tidx);
                            evnt_table = [];
                            for m=1:length(pks) %going highest to lowest
                                ntp_kde = SS.ntp_kde(ploc(m));
                                [~,pupilframe] = min(abs(SS.ntp_pupil - ntp_kde));
                                [~,goproframe] = min(abs(SS.ntp_gp - ntp_kde));
                                [~,npsample] = min(abs(SS.ntp_np - ntp_kde));
                                evnt_table = vertcat(evnt_table,table({"KDEPeakAll"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                if m<floor(length(pks)/2)
                                    evnt_table = vertcat(evnt_table,table({"KDEPeakHigh"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                else
                                    evnt_table = vertcat(evnt_table,table({"KDEPeakLow"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                end
                                if m<floor(length(pks)/3)
                                    evnt_table = vertcat(evnt_table,table({"KDEPeakHighTercile"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                elseif m>=floor(length(pks)/3) && m<2*floor(length(pks)/3)
                                    evnt_table = vertcat(evnt_table,table({"KDEPeakMidTercile"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                else
                                    evnt_table = vertcat(evnt_table,table({"KDEPeakLowTercile"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                end
                            end
                            for m=1:length(trs)
                                if trs(m)<pks5
                                    ntp_kde = SS.ntp_kde(tloc(m));
                                    [~,pupilframe] = min(abs(SS.ntp_pupil - ntp_kde));
                                    [~,goproframe] = min(abs(SS.ntp_gp - ntp_kde));
                                    [~,npsample] = min(abs(SS.ntp_np - ntp_kde));
                                    evnt_table = vertcat(evnt_table,table({"KDETrough"},{""},pupilframe,goproframe,npsample,ntp_kde,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                end
                            end
                            SS.evnts_tbl = vertcat(SS.evnts_tbl,evnt_table);
                            SS.evnts_tbl = sortrows(SS.evnts_tbl,"NTP");
                        end

                        %Full walk wavelet specgram (smoothed, downsampled, median normalized)
                        d = SS.d_np; nan_idx = isnan(d); d(nan_idx) = 0;
                        [f,~,cfs] = morseSpecGram(d,SS.fs_np,[2,120]); %2 to 120Hz
                        cfs(permute(repmat(nan_idx,1,1,length(f)),[1,3,2])) = nan; %time x freq x chan
                        pwr = abs(cfs).^2;
                        pwr = smoothdata(pwr,1,'movmean',250*obj.WinSegSec); %smooth across time (this handles nans better than smooth3)
                        pwr = pwr./median(pwr,1,"omitnan"); %full epoch norm (time x freq x chan)

                        obj.DTable(k,18) = {pwr(1:10:end,:,:)}; %smoothed/downsampled specgram power (factor of 10 so 25Hz)
                        obj.DTable(k,19) = {SS.ntp_np(1:10:end)};
                        obj.DTable(k,20) = {SS.fs_np./10}; %fs=25Hz
                        obj.DTable(k,21) = {f};

                        ptnum = str2double(regexp(obj.PatientID,'\d+$','match','once'));

                        SS.evnts_tbl =  addvars(SS.evnts_tbl,repmat(ptnum,size(SS.evnts_tbl,1),1),'NewVariableNames','Patient','Before','Event');
                        SS.evnts_tbl =  addvars(SS.evnts_tbl,repmat(WalkNum,size(SS.evnts_tbl,1),1),'NewVariableNames','Walk','Before','Event');

                        obj.DB = vertcat(obj.DB,SS.evnts_tbl);
                    end
                end
            end
        
            obj.DTable = cell2table(obj.DTable,'VariableNames',{...
                'd_np','ntp_np','fs_np',...
                'd_xs','ntp_xs','fs_xs',...
                'ntp_pupil','ntp_gp',...
                'd_gaze_fix','ntp_gaze','fs_gaze',...
                'd_amb','ntp_amb','fs_amb',...
                'd_kde','ntp_kde','fs_kde',...
                'd_wav','ntp_wav','fs_wav','f_wav',...
                });

        end %loadData

        function parseTransData(obj,varargin)
            %Parse data into +-12 sec segments centered on transitions.
            %Parses all transistions in the table. Run loadData first.
            %parseTransData;
            DT_np = []; %+-12sec transitions centered on specified transition type (i.e. 'Outdoor Beg')
            DT_np_idx = []; %index into data for center of transition
            DT_xs = [];
            EvntTrans = {}; %event type
            DescTrans = {}; %event description
            WalkNumTrans = []; %walk number for transitions
            NumSampTrans = []; %number of samples in data for each walk
            
            for m=1:size(obj.DTable,1) %by walk
                walknum = obj.WalkNums(m);
                numsamp = size(obj.DTable.d_np{m},1); %number of samples in walk data (needed to remove windows that extend beyond data limits)

                dbb = obj.DB(obj.DB.Walk==walknum,:);

                d_np = obj.DTable.d_np{m};
                ntp_np = obj.DTable.ntp_np{m}; %ntp in sec
                fs_np = obj.DTable.fs_np(m);
                if iscell(fs_np); fs_np = fs_np{1}; end

                d_xs = obj.DTable.d_xs{m};
                ntp_xs = obj.DTable.ntp_xs{m}; %ntp in sec
                fs_xs = obj.DTable.fs_xs(m);
                if iscell(fs_xs); fs_xs = fs_xs{1}; end
                if isempty(fs_xs); fs_xs = 100; end

                %12sec window around transition
                win_trans_sec = obj.WinTransSec; %sec
                ntp = dbb.NTP; %all transitions types
                evnt = dbb.Event;
                desc = dbb.Description;
                win_trans_np = (-win_trans_sec*fs_np:win_trans_sec*fs_np); d_np_trans = nan(length(win_trans_np),length(ntp),4); d_np_trans_idx = nan(length(ntp),1);
                win_trans_xs = (-win_trans_sec*fs_xs:win_trans_sec*fs_xs); d_xs_trans = nan(length(win_trans_xs),length(ntp));
                for k=1:length(ntp)
                    [~,midx] = min(abs(ntp_np-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_np(1))>=1 && (midx+win_trans_np(end))<=length(d_np)
                            d_np_trans_idx(k) = midx;
                            d_np_trans(:,k,:) = d_np(midx+win_trans_np,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of NP data. Filling with NaNs.\n',k,walknum);
                        end
                    end
                    [~,midx] = min(abs(ntp_xs-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_xs(1))>=1 && (midx+win_trans_xs(end))<=length(d_xs)
                            d_xs_trans(:,k) = d_xs(midx+win_trans_xs);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of XS data. Filling with NaNs.\n',k,walknum);
                        end
                    end
                end
                DT_np = cat(2,DT_np,d_np_trans);
                DT_np_idx = cat(1,DT_np_idx,d_np_trans_idx);
                DT_xs = cat(2,DT_xs,d_xs_trans);
                EvntTrans = cat(1,EvntTrans,evnt);
                DescTrans = cat(1,DescTrans,desc);
                WalkNumTrans = cat(1,WalkNumTrans,ones(size(d_np_trans,2),1).*walknum);
                NumSampTrans = cat(1,NumSampTrans,ones(size(d_np_trans,2),1).*numsamp);

            end %for DTable (by walk)

            obj.DParsed.TransSeg.DT_np = DT_np;
            obj.DParsed.TransSeg.DT_np_idx = DT_np_idx;
            obj.DParsed.TransSeg.DT_xs = DT_xs;
            obj.DParsed.TransSeg.TT_np_samp = win_trans_np;
            obj.DParsed.TransSeg.TT_np_sec = win_trans_np/fs_np;
            obj.DParsed.TransSeg.TT_xs_samp = win_trans_xs;
            obj.DParsed.TransSeg.TT_xs_sec = win_trans_xs/fs_xs;
            obj.DParsed.TransSeg.EvntTrans = EvntTrans;
            obj.DParsed.TransSeg.DescTrans = DescTrans;
            obj.DParsed.TransSeg.WalkNumTrans = WalkNumTrans;
            obj.DParsed.TransSeg.NumSampTrans = NumSampTrans;
            obj.DParsed.TransSeg.Outliers = any(any(isnan(DT_np),1),3); %<-500 outliers handled in loadData
            obj.DParsed.TransSeg.WinTransSec = obj.WinTransSec;
        end %parseTransData

        function parseSegData(obj,varargin)
            %Parse data into 4sec indoor/outdoor segments. Run loadData
            %first. This parses based on in/out boundaries. Make sense to
            %chunk up data sequentially from the start and tag segments as
            %in or out. Need to try this at some point...
            %parseSegData;
            DS_np = []; %4 sec segments marked as inside/outside
            DS_xs = []; %velocity data for xsens in 4sec segments to match DS_np
            DS_gaze = [];
            DS_amb = [];
            DS_kde = [];
            InFlag = logical([]); %flag for inside segments
            WalkNumSeg = []; %walk number for segments

            TT = {}; %in/out start/end times table

            for m=1:size(obj.DTable,1) %by walk
                walknum = obj.WalkNums(m);

                dbb = obj.DB(obj.DB.Walk==walknum,:);

                d_np = obj.DTable.d_np{m};
                ntp_np = obj.DTable.ntp_np{m}; %ntp in sec
                fs_np = obj.DTable.fs_np(m);
                if iscell(fs_np); fs_np = fs_np{1}; end

                d_xs = obj.DTable.d_xs{m};
                ntp_xs = obj.DTable.ntp_xs{m}; %ntp in sec
                fs_xs = obj.DTable.fs_xs(m);
                if iscell(fs_xs); fs_xs = fs_xs{1}; end
                if isempty(fs_xs); fs_xs = 100; end

                %%%%%%%%%%%% New Data %%%%%%%%%%%%%%%%%%
                d_gaze = obj.DTable.d_gaze_fix{m};
                ntp_gaze = obj.DTable.ntp_gaze{m}; %ntp in sec
                fs_gaze = obj.DTable.fs_gaze(m);
                if iscell(fs_gaze); fs_gaze = fs_gaze{1}; end
                if isempty(fs_gaze); fs_gaze = 200; end

                d_amb = obj.DTable.d_amb{m};
                ntp_amb = obj.DTable.ntp_amb{m}; %ntp in sec
                fs_amb = obj.DTable.fs_amb(m);
                if iscell(fs_amb); fs_amb = fs_amb{1}; end
                if isempty(fs_amb); fs_amb = 15; end

                d_kde = obj.DTable.d_kde{m};
                ntp_kde = obj.DTable.ntp_kde{m}; %ntp in sec
                fs_kde = obj.DTable.fs_kde(m);
                if iscell(fs_kde); fs_kde = fs_kde{1}; end
                if isempty(fs_kde); fs_kde = 60; end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %Segmenting data into 4sec segments and classifying as
                %indoor/outdoor. Only analyze segments that occur during
                %actual walk
                if any(contains(dbb.Event,'Walk Beg'))
                    npfirst = dbb.NPSample(find(contains(dbb.Event,'Walk Beg'),1));
                    nplast = dbb.NPSample(find(contains(dbb.Event,'Walk End'),1));
                    npfirst_ntp = dbb.NTP(find(contains(dbb.Event,'Walk Beg'),1));
                    nplast_ntp = dbb.NTP(find(contains(dbb.Event,'Walk End'),1));
                else
                    error('Walk Beg and Walk End annotations are missing!');
                end

                %Making table of indoor/outdoor times/samples
                dbb(~contains(dbb.Event,'Outdoor'),:) = [];
                if isempty(dbb)
                    error('Outdoor events are missing!');
                else
                    if ~(find(contains(dbb.Event,'Beg'),1)==1 && sum(contains(dbb.Event,'Beg'))==sum(contains(dbb.Event,'End')))
                        error('Indoor/outdoor events are incorrect!');
                    end
                    s = [[npfirst;dbb.NPSample],[dbb.NPSample;nplast]]; %np samples
                    t = [[npfirst_ntp;dbb.NTP],[dbb.NTP;nplast_ntp]]; %ntp
                    l = repmat({'Indoor';'Outdoor'},ceil(size(s,1)/2),1);
                    T = [l(1:size(s,1)),num2cell(s),num2cell(t)];
                    T = cell2table(T,"VariableNames",{'Location','Start_NP','End_NP','Start_NTP','End_NTP'});
                    T = addvars(T,repmat(walknum,size(T,1),1),'Before','Location','NewVariableNames','Walk');
                end
                TT = cat(1,TT,T);

                % Chunk into 4sec (1000sample) segments for inside/outside
                win_seg_sec = obj.WinSegSec; %default 4sec
                win_seg_np = win_seg_sec*fs_np;
                win_seg_xs = win_seg_sec*fs_xs;
                win_seg_gaze = win_seg_sec*fs_gaze;
                win_seg_amb = win_seg_sec*fs_amb;
                win_seg_kde = win_seg_sec*fs_kde;
                for k=1:size(T,1)
                    if (T.Start_NTP(k)-ntp_np(1))<0
                        error('NP start is after start of walk! This is a problem and needs to be checked.');
                    end

                    if ~isempty(d_xs)
                        if (T.Start_NTP(k)-ntp_xs(1))<0
                            error('XS start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_gaze)
                        if (T.Start_NTP(k)-ntp_gaze(1))<0
                            error('Gaze start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_amb)
                        if (T.Start_NTP(k)-ntp_amb(1))<0
                            error('Amb start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_kde)
                        if (T.Start_NTP(k)-ntp_kde(1))<0
                            error('KDE start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    %np data
                    d_np_seg = d_np(T.Start_NP(k):T.End_NP(k),:);
                    d_np_seg(floor(size(d_np_seg,1)/win_seg_np)*win_seg_np+1:end,:) = [];
                    d_np_seg = reshape(d_np_seg,win_seg_np,[],4);
                    DS_np = cat(2,DS_np,d_np_seg);

                    %xsens data (make the same size as d_np_seg)
                    d_xs_seg = nan(win_seg_xs,size(d_np_seg,2));
                    if ~isempty(d_xs)
                        [~,start_idx] = min(abs(ntp_xs - T.Start_NTP(k))); %start samples for xs
                        [~,end_idx] = min(abs(ntp_xs - T.End_NTP(k)));
                        dd_xs = d_xs(start_idx:end_idx);
                        dd_xs(floor(size(dd_xs,1)/win_seg_xs)*win_seg_xs+1:end,:) = []; %100Hz (4*100=400 -> 4sec)
                        dd_xs = reshape(dd_xs,win_seg_xs,[]);
                        if size(dd_xs,2)>size(d_np_seg,2)
                            d_xs_seg = dd_xs(:,1:size(d_np_seg,2));
                            fprintf('XS data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_xs,2)<size(d_np_seg,2)
                            d_xs_seg(:,1:size(dd_xs,2)) = dd_xs;
                            fprintf('XS data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_xs_seg = dd_xs;
                        end
                    end
                    DS_xs = cat(2,DS_xs,d_xs_seg);

                    %gaze data (make the same size as d_np_seg)
                    d_gaze_seg = nan(win_seg_gaze,size(d_np_seg,2));
                    if ~isempty(d_gaze)
                        [~,start_idx] = min(abs(ntp_gaze - T.Start_NTP(k))); %start samples for gaze
                        [~,end_idx] = min(abs(ntp_gaze - T.End_NTP(k)));
                        dd_gaze = d_gaze(start_idx:end_idx);
                        dd_gaze(floor(size(dd_gaze,1)/win_seg_gaze)*win_seg_gaze+1:end,:) = []; %200Hz (4*200=800 -> 4sec)
                        dd_gaze = reshape(dd_gaze,win_seg_gaze,[]);
                        if size(dd_gaze,2)>size(d_np_seg,2)
                            d_gaze_seg = dd_gaze(:,1:size(d_np_seg,2));
                            fprintf('Gaze data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_gaze,2)<size(d_np_seg,2)
                            d_gaze_seg(:,1:size(dd_gaze,2)) = dd_gaze;
                            fprintf('Gaze data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_gaze_seg = dd_gaze;
                        end
                    end
                    DS_gaze = cat(2,DS_gaze,d_gaze_seg);

                    %ambient light data (make the same size as d_np_seg)
                    d_amb_seg = nan(win_seg_amb,size(d_np_seg,2));
                    if ~isempty(d_amb)
                        [~,start_idx] = min(abs(ntp_amb - T.Start_NTP(k))); %start samples for amb
                        [~,end_idx] = min(abs(ntp_amb - T.End_NTP(k)));
                        dd_amb = d_amb(start_idx:end_idx);
                        dd_amb(floor(size(dd_amb,1)/win_seg_amb)*win_seg_amb+1:end,:) = []; %15Hz (4*15=60 -> 4sec)
                        dd_amb = reshape(dd_amb,win_seg_amb,[]);
                        if size(dd_amb,2)>size(d_np_seg,2)
                            d_amb_seg = dd_amb(:,1:size(d_np_seg,2));
                            fprintf('Ambient data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_amb,2)<size(d_np_seg,2)
                            d_amb_seg(:,1:size(dd_amb,2)) = dd_amb;
                            fprintf('Ambient data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_amb_seg = dd_amb;
                        end
                    end
                    DS_amb = cat(2,DS_amb,d_amb_seg);

                    %kde data (make the same size as d_np_seg)
                    d_kde_seg = nan(win_seg_kde,size(d_np_seg,2));
                    if ~isempty(d_kde)
                        [~,start_idx] = min(abs(ntp_kde - T.Start_NTP(k))); %start samples for kde
                        [~,end_idx] = min(abs(ntp_kde - T.End_NTP(k)));
                        dd_kde = d_kde(start_idx:end_idx);
                        dd_kde(floor(size(dd_kde,1)/win_seg_kde)*win_seg_kde+1:end,:) = []; %60Hz (4*60=240 -> 4sec)
                        dd_kde = reshape(dd_kde,win_seg_kde,[]);
                        if size(dd_kde,2)>size(d_np_seg,2)
                            d_kde_seg = dd_kde(:,1:size(d_np_seg,2));
                            fprintf('KDE data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_kde,2)<size(d_np_seg,2)
                            d_kde_seg(:,1:size(dd_kde,2)) = dd_kde;
                            fprintf('KDE data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_kde_seg = dd_kde;
                        end
                    end
                    DS_kde = cat(2,DS_kde,d_kde_seg);

                    WalkNumSeg = cat(1,WalkNumSeg,ones(size(d_np_seg,2),1).*walknum);
                    switch T.Location{k}
                        case 'Indoor'
                            InFlag = cat(1,InFlag,true(size(d_np_seg,2),1));
                        case 'Outdoor'
                            InFlag = cat(1,InFlag,false(size(d_np_seg,2),1));
                    end

                end %for size(T,1)

            end %for DTable

            obj.InOutTimes = TT;

            obj.DParsed.InOutSeg.DS_np = DS_np;
            obj.DParsed.InOutSeg.DS_xs = DS_xs;
            obj.DParsed.InOutSeg.DS_gaze = DS_gaze;
            obj.DParsed.InOutSeg.DS_amb = DS_amb;
            obj.DParsed.InOutSeg.DS_kde = DS_kde;
            obj.DParsed.InOutSeg.InFlag = InFlag;
            obj.DParsed.InOutSeg.WalkNumSeg = WalkNumSeg;
            obj.DParsed.InOutSeg.Outliers = any(any(isnan(DS_np),1),3); %<-500 outliers handled in loadData (set to nan)
            obj.DParsed.InOutSeg.WinSegSec = obj.WinSegSec;
        end %parseSegData
        
        function copyData(obj,varargin)
            %copies data from \\155.100.91.44\D\Data\RealWorldNavigationCory to
            %\\155.100.91.44\D\Tyler\RealWorld for analysis. This will
            %overwrite current files!
            wait_msg = parfor_wait(length(obj.PatientList));
            for m=1:length(obj.PatientList) %iterate patient
                wait_msg.Send;

                rootdir = '\\155.100.91.44\D\Data\RealWorldNavigationCory';
                patientid = regexp(obj.PatientList{m},'(RW)\d{1}','match','once');
                walknames = dir(fullfile(rootdir,obj.PatientList{m}));
                walknames = regexp({walknames.name},'Walk\d+','match','once');
                walknames(cellfun(@isempty,walknames)) = [];
                walknums = cellfun(@str2double,regexp(walknames,'\d+','match','once'));
                totalwalks = length(walknums);

                for k=1:totalwalks
                    fname = sprintf('RWNApp_%s_Walk%0.0f.mat',patientid,walknums(k));
                    src = fullfile(rootdir,obj.PatientList{m},walknames{k},fname);
                    dst = fullfile('\\155.100.91.44\D\Tyler\RealWorld',fname);
                    if isfile(src) && isfolder(fileparts(dst))
                        [status,msg] = copyfile(src,dst);
                        if ~status
                            disp(msg);
                        end
                    else
                        fprintf('Did not copy: %s\n',src);
                    end
                end
            end
            wait_msg.Destroy;

        end

        function loadVid(obj,varargin)
            %Load pupil video for current walk. If video is already loaded,
            %will not reload to save time.
            RunFlag = false;
            if isempty(obj.Vid)
                RunFlag = true;
            else
                if any(obj.PatientIdx~=obj.Vid.Patient)
                    RunFlag = true;
                end
            end
            if RunFlag
                obj.Vid = [];
                for k=1:length(obj.WalkNums)
                    fprintf('Loading video for walk %0.0f...\n',obj.WalkNums(k));
                    Files = findRWAFiles(fullfile('\\155.100.91.44\D\Data\RealWorldNavigationCory',obj.PatientID,'Original',obj.WalkNames{obj.WalkNums(k)}));
                    obj.Vid = vertcat(obj.Vid,table(obj.PatientIdx,obj.WalkNums(k),{VideoReader(Files.pupil_vid_file)},'variablenames',{'Patient','Walk','VidObj'}));
                end
            end
        end

        function getMultSegData(obj,varargin)
            %Generates a table of multi-taper spectrum data for 4sec
            %segments for all patients.
            SSTrials = []; %keep trials for testing in glme
            for k=1:length(obj.PatientList)
                disp(k);

                obj.loadData('PatientIdx',k);
                obj.parseSegData;

                DS = obj.DParsed.InOutSeg.DS_np;
                DS_walk = obj.DParsed.InOutSeg.WalkNumSeg;
                DS_xs = obj.DParsed.InOutSeg.DS_xs;
                DS_gaze = obj.DParsed.InOutSeg.DS_gaze; %boolean for fixation/no fixation
                DS_amb = obj.DParsed.InOutSeg.DS_amb;
                DS_kde = obj.DParsed.InOutSeg.DS_kde;
                InFlag = obj.DParsed.InOutSeg.InFlag;

                OL = obj.DParsed.InOutSeg.Outliers; %these have nans for <-500 condition
                DS(:,OL,:) = [];
                DS_walk(OL) = [];
                DS_xs(:,OL) = [];
                DS_gaze(:,OL) = [];
                DS_amb(:,OL) = [];
                DS_kde(:,OL) = [];
                InFlag(OL) = [];

                DS_xs(:,any(DS_xs>5)) = nan;

                params.Fs = obj.DTable.fs_np(1); %250
                params.trialave = 0;
                params.tapers = [3,5];
                params.fpass = [2,120];
                for m=1:4 %by channel
                    [S,f] = mtspectrumc(DS(:,:,m),params); S = S';
                    
                    SSTrials = vertcat(SSTrials,table(repmat(k,size(S,1),1),...
                        repmat(m,size(S,1),1),DS_walk,InFlag,S,median(DS_xs)',... %xs and amb are just median across 4sec window
                        (sum(DS_gaze)./size(DS_gaze,1))',median(DS_amb)',median(DS_kde)',... %gaze is percentage of time fixating in 4sec window
                        'VariableNames',...
                        {'Pt','Chan','Walk','InFlag','Pwr','Vel','Fix','Amb','KDE'}));
                end %chan

            end %patient

            %Table of individual 4sec segments for all data streams
            %(xs/fix/amb), NaNs exist for xs/fix/amp. Data normalized for
            %each stream to make normal. Distributions for Vel are a bit
            %skewed but not sure how to normalize.
            [uPtChan,~,uPtChanIdx] = unique(SSTrials(:,{'Pt','Chan'}),'rows');
            [uPtWalk,~,uPtWalkIdx] = unique(SSTrials(:,{'Pt','Walk'}),'rows');
            SSTrials = addvars(SSTrials,uPtChanIdx,'NewVariableNames','uPtChan','Before','Pt');
            SSTrials = addvars(SSTrials,uPtWalkIdx,'NewVariableNames','uPtWalk','Before','Pt');
            SSTrials = addvars(SSTrials,log10(SSTrials.Pwr+eps),'NewVariableNames','nPwr','After','Pwr');
            % SSTrials = addvars(SSTrials,atanh((SSTrials.Fix-0.5)*2),'NewVariableNames','nFix','After','Fix'); %fisher z transform to make normal (0 and 1 become infinite -> not sure how to handle this)
            % SSTrials = addvars(SSTrials,SSTrials.Fix,'NewVariableNames','nFix','After','Fix'); %fisher z transform to make normal (0 and 1 become infinite -> not sure how to handle this)
            SSTrials = addvars(SSTrials,log10(SSTrials.Amb+eps),'NewVariableNames','nAmb','After','Amb'); %still bimodal, consider a different normalization
            % SSTrials = addvars(SSTrials,SSTrials.KDE,'NewVariableNames','nKDE','After','KDE');

            % fH = figure;
            % fH.Name = 'nPwr';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     histogram(SSTrials.nPwr(idx),'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f',SSTrials.Pt(idx(1)),SSTrials.Chan(idx(1))))
            % end
            % 
            % fH = figure;
            % fH.Name = 'Vel';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     histogram(SSTrials.Vel(idx),'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f',SSTrials.Pt(idx(1)),SSTrials.Chan(idx(1))))
            % end
            % 
            % fH = figure;
            % fH.Name = 'nFix';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     histogram(SSTrials.nFix(idx),'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f',SSTrials.Pt(idx(1)),SSTrials.Chan(idx(1))))
            % end
            % 
            % fH = figure;
            % fH.Name = 'nAmb';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     histogram(SSTrials.nAmb(idx),'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f',SSTrials.Pt(idx(1)),SSTrials.Chan(idx(1))))
            % end
            % 
            % fH = figure;
            % fH.Name = 'nKDE';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     histogram(SSTrials.nKDE(idx),'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f',SSTrials.Pt(idx(1)),SSTrials.Chan(idx(1))))
            % end

            obj.MultInOut.SSTrials = SSTrials;
            obj.MultInOut.uPtChan = uPtChan;
            obj.MultInOut.uPtChanIdx = uPtChanIdx;
            obj.MultInOut.uPtWalk = uPtWalk;
            obj.MultInOut.uPtWalkIdx = uPtWalkIdx;
            obj.MultInOut.f = f;
            obj.MultInOut.params = params;
            obj.MultInOut.WinSegSec = obj.WinSegSec; %Also smoothing kernel for specgrams in trans analysis
            obj.MultInOut.WinTransSec = obj.WinTransSec;

        end %getMultSegData

        function getMultTransData(obj,varargin)
            %Generates multi-patient transition data for specgram analysis.
            %Power for all walks is generated here now and the specgram is
            %smoothed for faster processing in plotMultTransSpecGramPerm.
            %getMultTransData;
            DT = []; DT_idx = []; DT_xs = []; Walk = []; Evnt = []; Desc = []; StopWalk = []; 
            GoWalk = []; OL = []; Region = []; RegionLabel = []; NSamp = []; %number of data samples in walk
            Chan = []; ChanLabel = []; Patient = []; 
            MultDTable = cell(length(obj.PatientList),1); 
            MultWTable = cell(length(obj.PatientList),1); %wavelet transform
            for k=1:length(obj.PatientList)
                disp(k);
                obj.loadData('PatientIdx',k);
                obj.parseTransData;
                dt = obj.DParsed.TransSeg.DT_np; %time x trial x chan
                dt_idx = permute(repmat(obj.DParsed.TransSeg.DT_np_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_xs = repmat(obj.DParsed.TransSeg.DT_xs,1,1,4); %time x trial x chan
                walk = permute(repmat(obj.DParsed.TransSeg.WalkNumTrans(:),1,4),[3,1,2]); %1 x trial x chan
                nsamp = permute(repmat(obj.DParsed.TransSeg.NumSampTrans(:),1,4),[3,1,2]); %1 x trial x chan
                evnt = permute(repmat(obj.DParsed.TransSeg.EvntTrans(:),1,4),[3,1,2]); %1 x trial x chan
                desc = permute(repmat(obj.DParsed.TransSeg.DescTrans(:),1,4),[3,1,2]); %1 x trial x chan
                stopwalk = permute(repmat(obj.DParsed.TransSeg.WalkNumTrans(:)==obj.StopGoWalks{k}(1),1,4),[3,1,2]); %1 x trial x chan
                gowalk = permute(repmat(obj.DParsed.TransSeg.WalkNumTrans(:)==obj.StopGoWalks{k}(2),1,4),[3,1,2]); %1 x trial x chan
                ol = permute(repmat(obj.DParsed.TransSeg.Outliers(:),1,4),[3,1,2]); %nan and -500 only (1 x trial x chan)
                regiontable = obj.RegionTable(obj.RegionTable.Patient==k,:);
                region = permute(repmat(regiontable.Region(:)',size(dt,2),1),[3,1,2]); %1 x trial x chan
                regionlabel = permute(repmat(regiontable.RegionLabel(:)',size(dt,2),1),[3,1,2]); %1 x trial x chan
                chan = permute(repmat(regiontable.Chan(:)',size(dt,2),1),[3,1,2]); %1 x trial x chan
                chanlabel = permute(repmat(regiontable.ChanLabel(:)',size(dt,2),1),[3,1,2]); %1 x trial x chan
                pt = permute(repmat(k,size(dt,2),4),[3,1,2]); %1 x trial x chan

                DT = cat(2,DT,dt);
                DT_idx = cat(2,DT_idx,dt_idx);
                DT_xs = cat(2,DT_xs,dt_xs);
                Walk = cat(2,Walk,walk);
                NSamp = cat(2,NSamp,nsamp); %number of data samples in walk
                Evnt = cat(2,Evnt,evnt);
                Desc = cat(2,Desc,desc);
                StopWalk = cat(2,StopWalk,stopwalk);
                GoWalk = cat(2,GoWalk,gowalk);
                OL = cat(2,OL,ol);
                Region = cat(2,Region,region);
                RegionLabel = cat(2,RegionLabel,regionlabel);
                Chan = cat(2,Chan,chan);
                ChanLabel = cat(2,ChanLabel,chanlabel);
                Patient = cat(2,Patient,pt);

                MultDTable{k} = obj.DTable.d_np;
                MultWTable{k} = obj.DTable.d_wav;
                % for m=1:length(obj.DTable.d_np) %by walk
                %     d = obj.DTable.d_np{m}; nan_idx = isnan(d); d(nan_idx) = 0;
                %     [f,~,cfs] = morseSpecGram(d,250,[2,120]); %2 to 120Hz
                %     cfs(permute(repmat(nan_idx,1,1,length(f)),[1,3,2])) = nan; %time x freq x chan
                %     pwr = abs(cfs).^2;
                %     % smoothsize = [1000,1,1]+(rem([1000,1,1],2)==0); %smoothing window is 2sec across time only
                %     % pwr = smooth3(pwr,'box',smoothsize);
                %     pwr = smoothdata(pwr,1,'movmean',250*obj.WinSegSec); %smooth across time (this handles nans better than smooth3)
                %     % pwr = pwr./mean(pwr,1,"omitnan"); %full epoch norm (time x freq x chan)
                %     pwr = pwr./median(pwr,1,"omitnan"); %full epoch norm (time x freq x chan)
                %     MultWTable{k}{m,1} = pwr(1:10:end,:,:); %smoothed/downsampled specgram power (factor of 10 so 25Hz)
                % end

            end %patient

            DT = reshape(DT,size(DT,1),[]); %data (unwrapped by time, then trial, then chan)
            DT_idx = reshape(DT_idx,1,[]);
            DT_xs = reshape(DT_xs,size(DT_xs,1),[]);
            Walk = reshape(Walk,1,[]);
            NSamp = reshape(NSamp,1,[]); %number of data samples in walk
            Evnt = reshape(Evnt,1,[]);
            Desc = reshape(Desc,1,[]);
            StopWalk = reshape(StopWalk,1,[]);
            GoWalk = reshape(GoWalk,1,[]);
            OL = reshape(OL,1,[]);
            Region = reshape(Region,1,[]);
            RegionLabel = reshape(RegionLabel,1,[]);
            Chan = reshape(Chan,1,[]);
            ChanLabel = reshape(ChanLabel,1,[]);
            Patient = reshape(Patient,1,[]);
            
            obj.MultTrans.DT = DT;
            obj.MultTrans.DT_idx = DT_idx;
            obj.MultTrans.DT_xs = DT_xs;
            obj.MultTrans.Walk = Walk;
            obj.MultTrans.NSamp = NSamp;
            obj.MultTrans.Evnt = Evnt;
            obj.MultTrans.Desc = Desc;
            obj.MultTrans.StopWalk = StopWalk;
            obj.MultTrans.GoWalk = GoWalk;
            obj.MultTrans.OL = OL; %only includes <-500 and NaNs
            obj.MultTrans.OL2 = any(isnan(obj.MultTrans.DT),1)|isoutlier(max(abs(zscore(obj.MultTrans.DT)))); %NaNs plus simple IED detection
            obj.MultTrans.Region = Region;
            obj.MultTrans.RegionLabel = RegionLabel;
            obj.MultTrans.Chan = Chan;
            obj.MultTrans.ChanLabel = ChanLabel;
            obj.MultTrans.Patient = Patient;
            obj.MultTrans.TimeSamp = obj.DParsed.TransSeg.TT_np_samp;
            obj.MultTrans.TimeSec = obj.DParsed.TransSeg.TT_np_sec;
            obj.MultTrans.TimeSec_xs = obj.DParsed.TransSeg.TT_xs_sec;
            obj.MultTrans.MultDTable = MultDTable; %for permutation testing
            obj.MultTrans.MultWTable = MultWTable; %wavelet transform for permutation testing
            obj.MultTrans.Freq = obj.DTable.f_wav{1}; %freq vector (same for all datasets)
            obj.MultTrans.WinSegSec = obj.WinSegSec; %Also smoothing kernel for specgrams in trans analysis
            obj.MultTrans.WinTransSec = obj.WinTransSec;

        end %getMultSegData

        function calcfBOSCWalk(obj,varargin)
            %calcfBOSCWalk(1); %walk1
            %This runs fooof in python. If libiomp5md.dll error, add
            %KMP_DUPLICATE_LIB_OK=1 as a system variable in windows.
            walk = varargin{1};

            obj.FB.D = obj.DTable.d_np{walk}; %walk 1
            obj.FB.fs = obj.DTable.fs_np(walk);

            obj.FB.D(obj.FB.D<-500) = 0;
            obj.FB.D(isnan(obj.FB.D)) = 0;

            obj.FB.nt = size(obj.FB.D,1); %channels
            obj.FB.nc = size(obj.FB.D,2); %time
            obj.FB.t = (0:obj.FB.nt-1)/obj.FB.fs;
           
            obj.FB.data.label = {'chan1','chan2','chan3','chan4'};
            obj.FB.data.time = {obj.FB.t};
            obj.FB.data.trial = {obj.FB.D'};
            obj.FB.data.fsample = obj.FB.fs;

            % general fBOSC setup
            obj.FB.cfg.fBOSC.F                 = 2:0.5:20;
            obj.FB.cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
            obj.FB.cfg.fBOSC.fsample           = obj.FB.fs;         % current sampling frequency of EEG data

            % padding
            obj.FB.cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
            obj.FB.cfg.fBOSC.pad.detection_s   = 0.1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
            obj.FB.cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

            % fooof parameters - fit with fixed line or allow a knee
            obj.FB.cfg.fBOSC.fooof.aperiodic_mode    = 'knee'; %old = eBOSC not fooof
            obj.FB.cfg.fBOSC.fooof.version           = 'python'; %'matlab'

            % threshold settings
            obj.FB.cfg.fBOSC.threshold.duration	= repmat(10, 1, numel(obj.FB.cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
            obj.FB.cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold

            % episode post-processing
            obj.FB.cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now

            % general processing settings
            obj.FB.cfg.fBOSC.channel           = []; % select posterior channels (default: all)
            obj.FB.cfg.fBOSC.trial             = []; % select trials (default: all)
            obj.FB.cfg.fBOSC.trial_background  = []; % select trials for background (default: all)

            [obj.FB.fBOSC, obj.FB.cfg] = fBOSC_wrapper(obj.FB.cfg, obj.FB.data);

        end

        function [ts,chan] = findMarks(~, D)
            %Find threshold crossings in data (i.e. marks)
            thresh = -500;
            chan = find(~isoutlier(std(D,0,"omitnan")),1);
            dthr = diff([false;D(:,chan)<thresh]);
            [ts,tsidx] = find(dthr==1);
            if length(ts)<10
                ts = [];
                return;
            end
            tsend = find(dthr==-1);
            if length(ts)>length(tsend)
                while 1
                    idx = find((tsend-ts(1:length(tsend)))<0,1);
                    ts(idx) = [];
                    tsidx(idx) = [];
                    if (length(ts)==length(tsend))
                        break;
                    end
                    if isempty(idx)
                        ts(end) = [];
                        tsidx(end) = [];
                    end
                end
            end
            pk_max = diff([ts,tsend],1,2)>3; %discard wfs that are too broad
            wf_beg = ts<=30; %discard waveforms at the edges of the dataset
            wf_end = ts>=(size(D,1)-30); %discard waveforms at the edges of the dataset

            idx = pk_max | wf_beg | wf_end;

            ts(idx) = [];
            tsidx(idx) = [];
            tsend(idx) = [];

            if length(ts)<10
                ts = [];
                return;
            end

            dts = diff(ts); %pattern 4, 8 is a mark
            idx = ~all([dts(1:end-1)==4,dts(2:end)==8],2);
            idx = [idx;true(2,1)];

            ts(idx) = [];
            tsidx(idx) = [];
            tsend(idx) = [];

            if length(ts)<10
                ts = [];
                return;
            end
        end %findMarks

        function [MT,dd,dx] = filterMultTransData(obj,varargin)
            %Filters transition data based on transtype, regiontype, and
            %walktype.
            %
            %transtype -> 'Outdoor Beg', 'Doorway', etc. from event table
            %regiontype -> 'AntHipp', 'LatTemp', 'Ent+Peri', 'PostHipp+Para', 'All Chans'
            %patienttype -> 'All Patients', '1,2,5'
            %walktype -> 'First Walks', 'Last Walks', 'Stop Walks', 'Go Walks', 'All Walks', '2,3'
            %desctype -> keywords that subdivide a transition (closed, open, in2out&closed, etc.) (optional) (& can be used as an operator)
            %
            %filterMultTransData('transtype','Outdoor Beg','regiontype','PostHipp+Para','walktype','All Walks','patienttype','All Patients');
            %filterMultTransData('transtype','Doorway','regiontype','PostHipp+Para','walktype','All Walks','desctype','closed','patienttype','All Patients');
            p = obj.parseInputs(varargin{:});

            if isempty(p.transtype)
                error('transtype must be specified!');
            end
            if isempty(p.regiontype)
                error('regiontype must be specified!');
            end
            if isempty(p.walktype)
                error('walktype must be specified!');
            end
            if isempty(p.patienttype)
                error('patienttype must be specified!');
            end

            transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            customregion = p.customregion; %nx2 matrix where cols are patient,chan
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', [individual walks]
            desctype = p.desctype; %Bracketed keywords in description (close, open, etc.) (optional)
            patienttype = p.patienttype; %'All Patients', '2,3'
            veltype = p.veltype; %VelHigh, VelLow, VelHighTercile, VelMidTercile, VelLowTercile

            if isempty(obj.MultTrans)
                obj.getMultTransData;
            end

            if contains(transtype,'|')
                transcell = regexp(transtype,'\|','split');
                transidx = false(size(obj.MultTrans.Evnt));
                for k=1:length(transcell)
                    transidx = transidx | contains(obj.MultTrans.Evnt,transcell{k});
                end
            else
                transidx = contains(obj.MultTrans.Evnt,transtype);
            end
            if ~any(transidx)
                error('Transition type was incorrect!');
            end

            regionlist = {'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'}; %1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para (set in loadData, custom=pt/chan set manually)
            regioncell = regexp(regiontype,'\|','split');
            regionidx = false(size(regionlist));
            for k=1:length(regioncell)
                regionidx = regionidx | contains(regionlist,regioncell{k});
            end
            regionidx = find(regionidx);
            if isempty(regionidx)
                error('Region type was incorrect!');
            else
                if all(regionidx<5)
                    if ~all(ismember(regioncell,unique(obj.MultTrans.RegionLabel)))
                        error('Region type was incorrect!');
                    end
                elseif any(regionidx==6) %custom
                    if isempty(customregion)
                        error('Region type was "Custom" but customregion matrix is empty!'); %nx2 matrix where cols are patient,chan
                    end
                end
            end

            walklist = {'First Walks','Last Walks','Stop Walks','Go Walks','All Walks'};
            if str2num(walktype)
                walkidx = -1;
            else
                walkidx = find(strcmp(walklist,walktype),1);
            end
            if isempty(walkidx)
                error('Walk type was incorrect!');
            end

            patientlist = {'All Patients','Leave Out 1','Leave Out 2','Leave Out 3','Leave Out 4','Leave Out 5'};
            if str2num(patienttype)
                patientidx = -1;
            else
                patientidx = find(strcmp(patientlist,patienttype),1);
            end
            if isempty(patientidx)
                error('Patient type was incorrect!');
            end
            
            if isempty(desctype)
                descidx = true(size(obj.MultTrans.Desc));
            else
                desccell = regexp(desctype,'&|\|','split');
                if contains(desctype,'|')
                    descidx = false(size(obj.MultTrans.Desc));
                    for k=1:length(desccell)
                        descidx = descidx | ~cellfun(@isempty,regexp(obj.MultTrans.Desc,['\[',desccell{k},'\]']));
                    end
                else
                    descidx = true(size(obj.MultTrans.Desc));
                    for k=1:length(desccell)
                        descidx = descidx & ~cellfun(@isempty,regexp(obj.MultTrans.Desc,['\[',desccell{k},'\]']));
                    end
                end
                if ~any(descidx)
                    error('Desc type was incorrect!');
                end
            end

            %remove outliers
            % olidx = any(isnan(obj.MultTrans.DT),1)|isoutlier(max(abs(zscore(obj.MultTrans.DT))));
            % transidx = transidx & ~olidx;
            transidx = transidx & ~obj.MultTrans.OL2;

            dd = obj.MultTrans.DT(:,transidx);
            dd_idx = obj.MultTrans.DT_idx(transidx);
            dx = obj.MultTrans.DT_xs(:,transidx);
            pt = obj.MultTrans.Patient(transidx);
            rg = obj.MultTrans.Region(transidx); %region number 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
            wk = obj.MultTrans.Walk(transidx);
            swk = obj.MultTrans.StopWalk(transidx);
            gwk = obj.MultTrans.GoWalk(transidx);
            ch = obj.MultTrans.Chan(transidx);
            lb = obj.MultTrans.ChanLabel(transidx); %region label
            ns = obj.MultTrans.NSamp(transidx); %number of samples in walk
            ev = obj.MultTrans.Evnt(transidx);
            ds = obj.MultTrans.Desc(transidx);
            ds_idx = descidx(transidx);
            MT = table(pt',wk',swk',gwk',rg',lb',ch',dd_idx',ns',ev',ds',ds_idx','VariableNames',{'patient','walk','stopwalk','gowalk','regionnum','regionlabel','chan','dataidx','nsamp','evnt','desc','descidx'});

            %filter by region (amygdala/anterior hipp)
            if all(regionidx<5)
                idx = ~ismember(MT.regionnum,regionidx); %1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];
            elseif any(regionidx==6) %custom regions (nx2 matrix where cols are patient,chan)
                idx = false(size(MT,1),1);
                for k=1:size(customregion,1)
                    idx = idx | all(MT.patient==customregion(k,1) & MT.chan==customregion(k,2),2);
                end
                dd(:,~idx) = [];
                dx(:,~idx) = [];
                MT(~idx,:) = [];
            end

            %filter by walk
            if walkidx<5
                if walkidx==-1 %walktype is a number
                    idx = ~ismember(MT.walk,str2num(walktype)); %filter out all but specified walk number
                elseif walkidx<3
                    idx = (MT.walk>3); %keep 1st 3 walks
                    if walkidx==2 %keep last walks
                        idx = ~idx;
                    end
                elseif walkidx<4
                    idx = (MT.stopwalk~=1); %stop walks
                else
                    idx = (MT.gowalk~=1); %go walks
                end
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];
            end

            %filter by patient
            if patientidx==-1
                idx = ~ismember(MT.patient,str2num(patienttype)); %filter out all but specified walk number
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];
            elseif patientidx>1 %leave one out (indices 2 thru 6)
                switch patienttype
                    case 'Leave Out 1'
                        ptlist = setdiff(1:5,1);
                    case 'Leave Out 2'
                        ptlist = setdiff(1:5,2);
                    case 'Leave Out 3'
                        ptlist = setdiff(1:5,3);
                    case 'Leave Out 4'
                        ptlist = setdiff(1:5,4);
                    case 'Leave Out 5'
                        ptlist = setdiff(1:5,5);
                end
                idx = ~ismember(MT.patient,ptlist); %filter out all but specified walk number
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];
            end

            %filter by description
            idx = ~MT.descidx;
            dd(:,idx) = [];
            dx(:,idx) = [];
            MT(idx,:) = [];

            %filter by velocity
            if ~isempty(veltype)
                idx = isoutlier(mean(dx)) | any(isnan(dx));
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];

                t = obj.MultTrans.TimeSec_xs;
                t_idx = t>-5 & t<5; %only focus on change in vel for +-4sec
                mdx = median(dx);
                pk_min = min(dx(t_idx,:));
                vel_chg = (pk_min-mdx)./mdx; %percent change in velocity relative to median
                [~,sidx] = sort(vel_chg); %larger negative values represent largest decrease from median (sorted -1 to 0)

                mid_len = floor(length(sidx)/2); %middle
                ter_len = floor(length(sidx)/3); %tercile
                full_len = length(sidx);
                full_idx = 1:full_len;
                switch veltype
                    case 'VelHigh' %keep highest negative numbers
                        mid_idx = 1:mid_len;
                        idx = sidx(setdiff(full_idx,mid_idx)); %discard index
                    case 'VelLow'
                        mid_idx = mid_len+1:full_len;
                        idx = sidx(setdiff(full_idx,mid_idx));
                    case 'VelHighTercile'
                        ter_idx = 1:ter_len;
                        idx = sidx(setdiff(full_idx,ter_idx));
                    case 'VelMidTercile'
                        ter_idx = ter_len+1:ter_len*2;
                        idx = sidx(setdiff(full_idx,ter_idx));
                    case 'VelLowTercile'
                        ter_idx = ter_len*2+1:full_len;
                        idx = sidx(setdiff(full_idx,ter_idx));
                    otherwise
                        error('veltype is incorrect!')
                end
                dd(:,idx) = [];
                dx(:,idx) = [];
                MT(idx,:) = [];
            end

            disp('Selected transitions:')
            disp(unique(MT.evnt));
            disp('Selected regions:')
            disp(unique(MT.regionlabel));
            disp('Selected patients:')
            disp(unique(MT.patient)'); %pt
            disp('Selected walks:')
            disp(unique(MT.walk)'); %walk
            disp('Total trials:')
            disp(size(MT,1)); 

        end

        function [pwr,tsec,freq,np,tsec_np] = getFilteredMultTransData(obj,MT,tsamp,varargin)
            %MT is generated by filterMultTransData. tsamp is a vector of
            %samples centered at zero specifying the window size of the
            %specgram. The specgram data has been smoothed using a 4sec
            %boxcar and downsampled by a factor of 10 in getMultTransData.

            tsec = tsamp./25; %Fs=25Hz (downsampled by factor 10)
            freq = obj.MultTrans.Freq;

            tsamp_np = tsamp(1)*10:tsamp(end)*10; %25 to 250Hz upsample
            tsec_np = tsamp_np./250;

            ntime = length(tsamp);
            nfreq = length(freq);

            ntime_np = length(tsamp_np);

            [uPtWkCh,~,uPtWkChIdx] = unique(MT(:,{'patient','walk','chan'}));
            
            np = nan(size(uPtWkChIdx,1),ntime_np); %raw np data for running fooof, etc.
            pwr = nan(size(uPtWkChIdx,1),ntime,nfreq);
            for k=1:size(uPtWkCh,1)
                idx = find(k==uPtWkChIdx);
                mt = MT(idx,:);
                pt = uPtWkCh.patient(k);
                wk = uPtWkCh.walk(k);
                ch = uPtWkCh.chan(k);
                didx = round(mt.dataidx)+tsamp_np; 
                sidx = round(mt.dataidx/10)+tsamp; %factor 10 downsample (getMultTransData)
                ridx = any(didx<1,2)|any(didx>mt.nsamp(1),2); didx(ridx,:) = []; sidx(ridx,:) = []; idx(ridx) = []; %remove any windows that are outside data range
                dd = obj.MultTrans.MultDTable{pt}{wk}(:,ch); %time x chan
                dd = reshape(dd(didx,:),size(didx,1),size(didx,2));
                np(idx,:) = dd; %trial x time
                ss = obj.MultTrans.MultWTable{pt}{wk}(:,:,ch); %time x freq x chan
                ss = reshape(ss(sidx,:),size(sidx,1),size(sidx,2),nfreq);
                pwr(idx,:,:) = ss; %trial x time x freq
            end
            np = permute(np,[2,1]);
            pwr = permute(pwr,[2,3,1]); %time x freq x trial
        end

        function dispCategories(obj,varargin)
            regionlist = {'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'};
            walklist = {'First Walks','Last Walks','Stop Walks','Go Walks','All Walks'};
            permlist = {'standard','zscore'};
            correctionlist = {'cluster','pixel'};
            eventlist = unique(obj.MultTrans.Evnt);
            desclist = regexp(obj.MultTrans.Desc,'\[([a-zA-Z0-9]+)\]','tokens');
            desclist(cellfun(@isempty,desclist)) = [];
            desclist = cellfun(@(x)[x{:}],desclist,'UniformOutput',false);
            desclist = unique(cat(2,desclist{:}));
            disp('regiontype');
            disp(regionlist);
            disp(' ');
            disp('walktype');
            disp(walklist);
            disp(' ');
            disp('permtype');
            disp(permlist);
            disp(' ');
            disp('correctiontype');
            disp(correctionlist);
            disp(' ');
            disp('transtype');
            disp(eventlist);
            disp(' ');
            disp('desctype');
            disp(desclist);
            disp(' ');
        end

        function saveMultData(obj,varargin)
            disp('Saving...')
            multinout = obj.MultInOut;
            multtrans = obj.MultTrans;
            save(obj.AnalysisFile,'multinout','multtrans','-v7.3');
            disp('Finished saving.')
        end

        function loadMultData(obj,varargin)
            if isfile(obj.AnalysisFile) || isfile('RWAnalysis.mat')
                if isfile('RWAnalysis.mat')
                    analysisfile = 'RWAnalysis.mat';
                else
                    analysisfile = obj.AnalysisFile;
                end
                disp('Loading data for multiple patients from file...')
                MD = load(analysisfile,'multinout','multtrans');
                if ~isempty(MD.multinout)
                    obj.MultInOut = MD.multinout;
                end
                if ~isempty(MD.multtrans)
                    obj.MultTrans = MD.multtrans;
                end
            else
                disp('Analysis file not found!')
            end
        end

        end %methods

    %%%%%%%%%%%%%%%%%%%%%
    methods %general plotting

        function plotInOutDiff(obj,varargin)
            %Plots spectrum comparison of inside vs outside
            %plotInOutDiff('walknum',1); only walk 1
            %plotInOutDiff('fpass',[2,85]); different fpass
            %plotInOutDiff; all walks
            p = obj.parseInputs(varargin{:});
            
            DS = obj.DParsed.InOutSeg.DS;
            InFlag = obj.DParsed.InOutSeg.InFlag;
            WalkNumSeg = obj.DParsed.InOutSeg.WalkNumSeg;

%             DS_xs = obj.DParsed.InOutSeg.DS_xs;
%             DS_xs(DS_xs>5) = nan;
%             DS_xs = mean(DS_xs,'omitnan');
% 
%             nan_idx = isnan(DS_xs);
%             DS_xs(nan_idx) = [];
%             DS(:,nan_idx,:) = [];
%             InFlag(nan_idx) = [];
%             WalkNumSeg(nan_idx) = [];
% 
%             InFlag = DS_xs>1; %velocity greater than 1

            if ~isempty(p.walknum)
                idx = WalkNumSeg~=p.walknum;
                DS(:,idx,:) = [];
                InFlag(idx) = [];
            end

            if isempty(p.fpass)
                fpass = [2,20];
            else
                fpass = p.fpass;
            end

%             InFlag = InFlag(randperm(length(InFlag)));
            
            params.Fs = 250;
            params.trialave = 0;
            params.err = [2,0.05];
            params.tapers = [3,5];
            params.fpass = [fpass(1)-1,fpass(2)+1];

            fH = figure('Position',[50,50,1700,500],'Visible','on');
            for k=1:4 %by channel
                disp(k);

                params.trialave = 0;
                [S,f] = mtspectrumc(DS(:,:,k),params);
                [PVal,PMsk] = kstestPerm_mex(single(permute(S,[2,3,1])),InFlag(:)',200,0.05);

                SE = strel('rectangle',[3 1]);
                PMsk = imdilate(PMsk,SE);
                PMsk = imerode(PMsk,SE);
                B = regionprops(PMsk(:)','BoundingBox');
                B = cat(1,B.BoundingBox);
%                 disp(B(:,3))
%                 if ~isempty(B)
%                     B(B(:,3)<3,:) = [];
%                 end

                params.trialave = 1;
                [S_in,~,~,Serr_in] = mtspectrumc(DS(:,InFlag,k),params);
                [S_out,f,~,Serr_out] = mtspectrumc(DS(:,~InFlag,k),params);

                Serr_in = log10(Serr_in);
                Serr_out = log10(Serr_out);

                ylimit = [min(reshape([Serr_in,Serr_out],[],1)),max(reshape([Serr_in,Serr_out],[],1))];
                ylimit = [ylimit(1)-diff(ylimit)*0.1,ylimit(2)+diff(ylimit)*0.1];

                aH = subplot(1,4,k,'parent',fH);
                hold(aH,'on');
                pH = nan(size(B,1),1);
                for m=1:size(B,1)
                    b = B(m,[1,3]);
                    x = floor([b(1),sum(b),sum(b),b(1)]);
                    x(x<1) = 1; x(x>length(f)) = length(f);
                    y = [ylimit(1),ylimit(1),ylimit(2),ylimit(2)];
                    pH(m) = patch(f(x),y,'k','FaceAlpha',0.05,'EdgeColor','none','parent',aH);
                end
                p1 = patch([f,fliplr(f)],([Serr_in(1,:),fliplr(Serr_in(2,:))]),'b','facealpha',0.3,'edgecolor','none','parent',aH);
                p2 = patch([f,fliplr(f)],([Serr_out(1,:),fliplr(Serr_out(2,:))]),'r','facealpha',0.3,'edgecolor','none','parent',aH);
                axis(aH,[fpass,ylimit])
                box(aH,'on');       
                xlabel(aH,'Hz')
                ylabel(aH,'log10(uV2/Hz)')
                title(aH,[{['Chan',num2str(k)]},obj.ChanLabels(obj.PatientIdx,k)],"FontSize",10)
                legend(aH,[p1,p2],{'Inside','Outside'})

            end

            fstr = sprintf('%s_AllValidWalks_Spectrum_%0.0f-%0.0fHz_Inside-Outside_n%0.0f-n%0.0f.png',obj.PatientID,fpass(1),fpass(2),sum(InFlag),sum(~InFlag));
            print(fH,fullfile('C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs',obj.PatientID,fstr),'-dpng','-r300');
            close(fH);
        end

        function plotTransSpecGramGUI(obj,varargin)
            %Plots spectrogram for specified patient/walk/transition along
            %with video and allows user to scroll through video frames with
            %current time marked in spectrogram. Use left/right arrows to
            %scroll. Ctrl+arrow advances by 100 frames. Shift+arrow
            %advances by 10 frames.
            %
            %plotTransSpecGramGUI('patient',1,'walknum',1,'evntnum',1,'transtype','Outdoor Beg');

            p = parseInputs(obj,varargin{:});
            if isempty(p.patient)
                error('Need to specify patient!');
            end
            if isempty(p.walknum)
                error('Need to specify walknum!');
            end
            if isempty(p.evntnum)
                error('Need to specify evntnum!');
            end
            if isempty(p.transtype)
                error('Need to specify transtype!');
            end
            pt = p.patient;
            wlk = p.walknum;
            evnt = p.evntnum;
            trans = p.transtype; %'Outdoor Beg', etc

            if all(pt~=obj.PatientIdx)
                disp('loading new patient data...')
                obj.loadData('PatientIdx',pt,'transtype',trans);
            end

            obj.loadVid; %loads videos for all walks for a patient
            vid = obj.Vid.VidObj{wlk};
                
            db = obj.DB;
            db = db(contains(db.Event,trans),:);
            db = db(db.Walk==wlk,:);
            db = db(evnt,:);

            if any(db.Patient~=pt) || any(db.Walk~=wlk)
                error('Patient or walk do not match!');
            end

            winsec = 32;

            fs_np = obj.DTable.fs_np(wlk);
            nsamp = winsec*fs_np;
            tsamp = (-nsamp:nsamp); %+-32 sec (2sec removed for artifact, so +-30sec)
            tsec = tsamp/fs_np;
            tfull = db.NPSample+tsamp; %+-12sec

            nfr = round(vid.FrameRate*winsec);
            tframes = (-nfr:nfr)+db.PupilFrame; %frames matching tfull above
            tfrsec = (-nfr:nfr)./vid.FrameRate;
            cf = 1;

            fH = figure('Position',[50,50,1700,500],'Visible','on');
            fH.Name = sprintf('Pt%0.0f, Wk%0.0f, Evnt%0.0f',pt,wlk,evnt);
            aH = nan(1,4); pH = nan(1,4);
            for k=1:4 %chan
                aH(k) = subplot(4,1,k,'parent',fH);
                d = obj.DTable.d_np{wlk}(tfull,k);
                d(d<-500) = 0;
                d(isnan(d)) = 0;
                lb = obj.ChanLabels{pt,k};
                PSG = PermSpecGram(d,'Fs',250,'FPass',[2,64],'TimeRng',[tsec(1),tsec(end)],'NormRng',[tsec(1)+2,tsec(end)-2]);
                PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','Smoothing',[1000,5],'ErrPerc',[]);
                PSG.plotSpecGram('ShowRaw',true,'Clim',[-3,3],'aH',aH(k));
                plot(aH(k),[0,0],[0,8],'k'); %center line
                pH(k) = plot(aH(k),[tsec(1)+2,tsec(1)+2],[0,8],'k'); %video line
                xlim(aH(k),[tsec(1)+2,tsec(end)-2])
                title(aH(k),lb)
            end

            fH2 = figure;
            aH2 = axes('parent',fH2);
            iH2 = image(read(vid,tframes(cf)),'parent',aH2);
            axis(aH2,'image')
            set(aH2,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
            title(aH2,sprintf('Frame = %0.0f',tframes(cf)))

            SS = [];
            SS.fH = fH;
            SS.aH = aH;
            SS.pH = pH;
            SS.fH2 = fH2;
            SS.aH2 = aH2;
            SS.iH2 = iH2;
            SS.cf = cf;
            SS.tframes = tframes;
            SS.tfrsec = tfrsec;
            SS.vid = vid;

            set(fH,'KeyPressFcn',@obj.keyEvntFcn)
            setappdata(fH,'SS',SS)
            set(fH2,'KeyPressFcn',@obj.keyEvntFcn)
            setappdata(fH2,'SS',SS)
        end

        function openMultTransSpecGramGUI(obj,varargin)
            MultTransSpecGramGUI(obj);
        end

        function varargout = plotMultTransSpecGram(obj,varargin)
            %%%%%% Phasing this out in GUI %%%%%%
            %Plot specgram of transitions for all patients and channels
            %using PermSpecGram.m
            %plotMultTransSpecGram('transtype','Outdoor Beg','regiontype','PostHipp+Para','walktype','All Walks','patienttype','All Patients');
            %plotMultTransSpecGram('transtype','Outdoor Beg','regiontype','AntHipp','walktype','Stop Walks','patienttype','All Patients','plottrials',true,'trialsfreqrng',[5,8]);

            [MT,dd,dx] = filterMultTransData(obj,varargin{:});

            p = obj.parseInputs(varargin{:});
            
            transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,3'
            patienttype = p.patienttype; %'All Patients', '1,3'
            desctype = p.desctype; %close, open, etc. (optional, will skip if empty)
            pval = p.pval; %default p=0.05
            normrng = p.normrng; %default [-10,10]
            smoothwin = p.smoothwin; %default 4sec
            climit = p.clim;
            fpass = [2,120];
            plottrials = p.plottrials;
            trialsfreqrng = p.trialsfreqrng;

            tsec = obj.MultTrans.TimeSec;
            tsec_xs = obj.MultTrans.TimeSec_xs;
            ntrials = size(MT,1);

            PSG = PermSpecGram(dd,'Fs',250,'FPass',fpass,'TimeRng',[tsec(1),tsec(end)],'NormRng',normrng);
            PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','Smoothing',[smoothwin*250,1],'ErrPerc',pval);
            
            fH = figure('Position',[50,50,1200,800]);
            aH = axes('parent',fH);

            PlotStyle = 'contour'; %'image','contour'
            if isempty(climit)
                CLimit = [-1,1]; %works well for dB
            else
                CLimit = climit;
            end
            PSG.plotSpecGram('MaskFlag',true,'MaskAlpha',0.2,'Clim',CLimit,'aH',aH,'PlotStyle',PlotStyle);

            cb = colorbar(aH);
            dxx = mean(dx(:,~isoutlier(mean(dx))),2,"omitnan");
            switch PlotStyle
                case 'contour'
                    set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'yminortick','off');
                    plot(aH,[0,0],fpass,'k');
                    plot(aH,tsec_xs,(dxx-median(dxx))*9+10,'k','LineWidth',2);
                    if plottrials
                        plot(aH,[tsec(1)+2,tsec(end)-2],repmat(trialsfreqrng,2,1),':k')
                    end
                    cblims = [cb.Limits(1),0,cb.Limits(2)];
                    cb.Ticks = cblims;
                    cb.TickLabels = cblims;
                case 'image'
                    plot(aH,[0,0],log2(fpass),'k');
                    plot(aH,tsec_xs,(dxx-median(dxx))*2.5+3.5,'k','LineWidth',2);
                    if plottrials
                        plot(aH,[tsec(1)+2,tsec(end)-2],repmat(log2(trialsfreqrng),2,1),':k')
                    end
                    cb.Limits = [cb.Limits(1),cb.Limits(2)-1];
                    cblims = [cb.Limits(1),(diff(cb.Limits)-1)/2,cb.Limits(2)-1];
                    cb.Ticks = cblims;
                    cb.TickLabels = [CLimit(1),0,CLimit(2)];
            end
            xlim(aH,[tsec(1)+2,tsec(end)-2])
            ylabel(cb,'dB')

            ptstr = regexprep(num2str(unique([MT.patient])'),'\s+','/');
            ttlstr = sprintf('%s, %s, %s, %s, %s, %s, p<%1.2f',...
                transtype,desctype,regiontype,ptstr,...
                num2str(walktype),num2str(ntrials),pval);
            title(aH,ttlstr); 

            if plottrials
                fH2 = figure('Position',[50,50,1200,800]);
                aH2 = axes('parent',fH2);
                colormap('jet');

                [mt,sidx] = sortrows(MT,{'patient','chan','walk'});
                fidx = PSG.fSpec>=trialsfreqrng(1) & PSG.fSpec<=trialsfreqrng(2);
                fdat = squeeze(mean(PSG.Stb(:,fidx,sidx),2))';

                cl = round(prctile(abs(fdat(:)),99));
                CLimit2 = [-cl,cl];
                II = fix((fdat-CLimit2(1))./(CLimit2(2)-CLimit2(1))*255)+1;
                II(II<1) = 1; II(II>256) = 256;
                iH = imagesc(PSG.tSpec,1:size(fdat,1),II,'parent',aH2);
                iH.CDataMapping = "direct";

                xlim(aH2,[tsec(1)+2,tsec(end)-2])
                cb2 = colorbar(aH2);
                cb2.Limits = [cb2.Limits(1),cb2.Limits(2)-1];
                cblims2 = [cb2.Limits(1),(diff(cb2.Limits)-1)/2,cb2.Limits(2)-1];
                cb2.Ticks = cblims2;
                cb2.TickLabels = [CLimit2(1),0,CLimit2(2)];
                ylabel(cb2,'dB')
                title(aH2,sprintf('Trials %0.0f-%0.0f Hz (%s)',trialsfreqrng(1),trialsfreqrng(2),ttlstr))
                [uPt,~,uPtIdx] = unique(mt.patient);
                mYtick = nan(length(uPt),1);
                hold(aH2,'on');
                for k=1:length(uPt)
                    idx = find(uPtIdx==k);
                    mYtick(k) = mean(idx);
                    plot([tsec(1)+2.3,tsec(1)+2.3],[idx(1),idx(end)],'color',[0.5,0.5,0.5],'LineWidth',10)
                end
                set(aH2,'ytick',mYtick,'yticklabel',uPt)
                ylabel(aH2,'patient')
                xlabel(aH2,'sec')
            end

            if nargout
                varargout{1} = fH;
            end

            if nargout>1
                varargout{2} = fH2;
            end

%             rootdir = 'C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs\Combined\InOut_Transition_Specgram\StopGo';
%             fstr = sprintf('RW1-RW5_%s_%s_%s_Transition_Specgram.png',regexprep(regiontype,'\s+',''),regexprep(walktype,'\s+',''),regexprep(transtype,'\s+',''));
%             print(fH,fullfile(rootdir,fstr),'-dpng','-r300');
%             close(fH);

        end

        function varargout = plotMultTransSpecGramPerm(obj,varargin)
            %Plot specgram of transitions for all patients and channels
            %using smoothed/downsampled specgram data generated in
            %getMultTransData. Run dispCategories to see options for
            %transtype, regiontype, walktype, desctype, permtype, and
            %correctiontype. Baseline is rotational shuffle across trials.
            %
            %plotMultTransSpecGramPerm('transtype','Outdoor Beg','regiontype','PostHipp+Para','walktype','All Walks','permtype','standard','correctiontype','cluster','patienttype','All Patients');
            %plotMultTransSpecGramPerm('transtype','Doorway','regiontype','PostHipp+Para','walktype','All Walks','desctype','closed','permtype','standard','correctiontype','cluster','patienttype','All Patients');
            %plotMultTransSpecGramPerm('transtype','Doorway','desctype','in2out&closed','regiontype','AntHipp','walktype','All Walks','permtype','zscore','correctiontype','cluster','patienttype','All Patients');
            %plotMultTransSpecGramPerm('transtype','Correct Turn Beg','regiontype','AntHipp','walktype','All Walks','permtype','zscore','correctiontype','cluster','patienttype','All Patients');
            %plotMultTransSpecGramPerm('transtype','Correct Turn Beg','regiontype','PostHipp+Para','walktype','All Walks','permtype','zscore','correctiontype','cluster','patienttype','All Patients');
            %plotMultTransSpecGramPerm('transtype','Incorrect Turn Beg','regiontype','PostHipp+Para','walktype','All Walks','permtype','zscore','correctiontype','cluster','patienttype','All Patients');
            p = obj.parseInputs(varargin{:});

            if isempty(p.permtype)
                error('permtype must be specified!');
            end

            warning('off','MATLAB:contour:ConstantData');
            
            transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,5'
            veltype = p.veltype; %'', 'High Change', 'Low Change'
            patienttype = p.patienttype; %'All Patients', '2,4'
            desctype = p.desctype; %close, open, etc. (optional, will skip if empty)
            permtype = p.permtype; %standard, zscore
            correctiontype = p.correctiontype; %cluster, pixel, fdr, fdr+cluster (if empty, cluster correction is skipped)
            if strcmp(correctiontype,'fdr+cluster') && strcmp(permtype,'standard')
                error('fdr+cluster cannot be performed with the standard permtype!');
            end
            transrng = p.transrng; %default [-10,10]
            normrng = p.normrng; %default [-10,10] must be within transrng
            if ~p.fullwalknorm
                if normrng(1)<transrng(1) || normrng(2)>transrng(2) || normrng(2)<=normrng(1)
                    error('normrng is incorrect!');
                end
            end
            pval = p.pval; %default p=0.05
            pvalclust = p.pvalclust; %default 0.01
            plottrials = p.plottrials;
            trialsfreqrng = p.trialsfreqrng;
            plotboxcomp = p.plotboxcomp;
%             plotboxcomp = true;
            boxcomprng = p.boxcomprng; %2x4 matrix where rows are each box and cols are low/high ranges for time/freq
%             boxcomprng = [-4,0,24,40;2,6,24,40];
            if isempty(p.clim)
                switch permtype
                    case 'zscore'
                        climit = [-10,10];
                    case 'standard' %dB
                        climit = [-1,1];
                end
            else
                climit = p.clim;
            end

            %Init some params
            tsamp = round(transrng(1)*25):round(transrng(2)*25); %+-10sec at 25Hz (downsampled by 10 from 250)
            tsec = tsamp./25;
            baseidx = tsec>=normrng(1) & tsec<=normrng(2); %baseline (normalization) indices
            nperm = 1000; %permutations

            %Getting trials and specgram data
            [MT,~,dx] = obj.filterMultTransData(varargin{:}); %table of all trials and corresponding info
            [pwr,tsec,freq,np,tsec_np] = obj.getFilteredMultTransData(MT,tsamp); %raw power for all trials (time x freq x trial)

            nfreq = length(freq);
            ntime = length(tsamp);
            ntrials = size(MT,1);

            %Real
            if p.fullwalknorm %Normalization is done across the entire walk for each chan in getMultTransData, so skip normalization across window
                bpwr = ones(1,nfreq); %Since the entire walk is baseline normalized, set to ones for full window norm to preserve relative power
            else
                bpwr = squeeze(mean(mean(pwr(baseidx,:,:),1,"omitnan"),3,"omitnan")); %mean across time and then trials (1 x freq)
                pwr = pwr./mean(pwr,1,"omitnan"); %full epoch norm (time x freq x trial) 
            end
            mpwr = squeeze(mean(pwr,3,"omitnan")); %trial avg (time x freq)
            npwr = 10*log10(mpwr./bpwr); %normalized power in dB
            npwr_trials = 10*log10(pwr./bpwr); %normalized power for individual trials

            %Permutations
            disp('Running permutations...')
            % PM = calcRWAPerm(pwr(baseidx,:,:),bpwr,0); %trialtype=0 for single specgram (vector of -1 or 1 for all trials when computing specgram difference)
            PM = calcRWAPerm_mex(pwr(baseidx,:,:),bpwr,0); %1.3 times faster

            %Finding percentiles of permuted distributions
            pm_freq = reshape(permute(PM,[1,3,2]),[],nfreq); %(ntime*nperm) x nfreq
            pc = prctile(pm_freq,[pval*100/2,100-pval*100/2],1); %2 x nfreq
            mPM = mean(pm_freq,1); %mean across permutations and time (1 x nfreq)
            sPM = std(pm_freq,0,1); %std across permutations and time (1 x nfreq)
%             %Keep time (similar to the zscore method)
%             pm = reshape(PM,[],nperm); %finding percentiles for all time x freq distributions (instead of collapsing across time like above)
%             pc_low = reshape(prctile(pm,pval*100/2,2),ntime,nfreq); %two tail (across 2nd dim -> 1000 perms)
%             pc_high = reshape(prctile(pm,100-pval*100/2,2),ntime,nfreq); %two tail
%             pc = permute(cat(3,pc_low,pc_high),[3,1,2]); %2 x ntime x nfreq
%             mPM = mean(PM,3); %mean across permutations
%             sPM = std(PM,0,3); %std across permutations

            %Finding clusters/max pixels in permuted data
            max_clust_info = nan(nperm,1);
            max_pixel_pvals = zeros(nperm,2);
            for m=1:nperm
                pm = PM(:,:,m);
                switch permtype
                    case 'standard'
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        pm(pm>squeeze(pc(1,:,:)) & pm<squeeze(pc(2,:,:))) = 0;
                    case 'zscore'
                        pm = (pm-mPM)./sPM;
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        switch correctiontype
                            case 'fdr+cluster'
                                PV = 1-normcdf(abs(pm));
                                PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                                pm(PV>=pval) = 0;
                            otherwise
                                pm(abs(pm)<norminv(1-pval)) = 0;
                        end
                end
                clustinfo = bwconncomp(pm);
                stats = regionprops(clustinfo,pm,'pixelvalues');
                max_clust_info(m) = max([0,abs(cellfun(@sum,{stats.PixelValues}))]); %max cluster size using pixel values for each permutation
                % max_clust_info(m) = max([0,cellfun(@numel,clustinfo.PixelIdxList)]); %max cluster size for each permutation
            end
            
            %Finding percentiles of pixel-level distributions
            pc_pixel = prctile(max_pixel_pvals(:),[pval*100/2,100-pval*100/2]); %two tail

            %Thresholding real
            switch permtype
                case 'standard'
                    npwr_thresh = npwr;
                    switch correctiontype
                        case 'fdr' %results are slightly different than PermSpecGram because this operates on log10 data and does rotating shuffle
                            %Finding pvals based on rank in distribution for FDR correction
                            % PV = calcRWAPVal(npwr,pm_freq);
                            PV = calcRWAPVal_mex(npwr,pm_freq); %~3x faster
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(npwr));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(npwr_thresh>squeeze(pc(1,:,:)) & npwr_thresh<squeeze(pc(2,:,:))) = 0;
                    end
                case 'zscore'
                    npwr = (npwr-mPM)./sPM;
                    npwr_thresh = npwr;
                    npwr_trials = (npwr_trials-mPM)./sPM;
                    switch correctiontype
                        case {'fdr','fdr+cluster'}
                            PV = 1-normcdf(abs(npwr));
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(abs(npwr)<norminv(1-pval))=0; %norminv gives stat value at pval for normal distribution (i.e. p=0.05 is -1.65)
                    end
            end

            %Removing clusters
            if contains(correctiontype,'cluster') && ~isempty(correctiontype)
                if size(PM,1)~=ntime
                    error('Cluster correction cannot be performed if the baseline/normalization window is smaller than the transition window!');
                else
                    clustinfo = bwconncomp(npwr_thresh);
                    stats = regionprops(clustinfo,npwr_thresh,'pixelvalues');
                    clust_info = abs(cellfun(@sum,{stats.PixelValues}));
                    % clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                    clust_threshold = prctile(max_clust_info,100-pvalclust*100);
                    whichclusters2remove = find(clust_info<clust_threshold);
                    for k=1:length(whichclusters2remove)
                        npwr_thresh(clustinfo.PixelIdxList{whichclusters2remove(k)})=0;
                    end
                end
            end

            %Plotting
            fH = figure('Position',[50,50,1200,800]);
            aH = axes('parent',fH);
            colormap("jet");
            contourf(tsec,freq,npwr',100,'linecolor','none','parent',aH);
            hold(aH,"on");
            contour(tsec,freq,(npwr_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
            set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'yminortick','off','clim',climit); %now in units of standard deviation
            plot(aH,[0,0],[2,120],'--k','LineWidth',2);
            dxx = dx(:,~isoutlier(mean(dx)));
            dxx = dxx - mean(dxx);
            dxx_me = mean(dxx,2,"omitnan");
            dxx_se = std(dxx,0,2,'omitnan');
            if p.includevel
                plot(aH,obj.MultTrans.TimeSec_xs,(dxx_me)*9+10,'k','LineWidth',2);
                plot(aH,obj.MultTrans.TimeSec_xs,((dxx_me-dxx_se))*9+10,':k','LineWidth',0.5);
                plot(aH,obj.MultTrans.TimeSec_xs,((dxx_me+dxx_se))*9+10,':k','LineWidth',0.5);
            end
            if plottrials
                plot(aH,[tsec(1),tsec(end)],repmat(trialsfreqrng,2,1),':k')
            end
            if plotboxcomp
                rectangle(aH,'Position',[boxcomprng(1,1),boxcomprng(1,3),diff(boxcomprng(1,1:2)),diff(boxcomprng(1,3:4))],'LineStyle','-.');
                rectangle(aH,'Position',[boxcomprng(2,1),boxcomprng(2,3),diff(boxcomprng(2,1:2)),diff(boxcomprng(2,3:4))],'LineStyle','--');
            end
            cb = colorbar(aH);
            cblims = [cb.Limits(1),0,cb.Limits(2)];
            cb.Ticks = cblims;
            cb.TickLabels = cblims;
            xlim(aH,[tsec(1),tsec(end)])
            ylim(aH,[2,120])
            if contains(permtype,'zscore')
                ylabel(cb,'Zscore')
            else
                ylabel(cb,'dB')
            end
            xlabel(aH,'sec');
            ylabel(aH,'Hz');
            if isempty(str2num(patienttype))
                pttypestr = patienttype;
            else
                pttypestr = 'Sel Patients';
            end
            ptstr = regexprep(num2str(unique([MT.patient])'),'\s+','');
            if isempty(str2num(walktype))
                wktypestr = walktype;
            else
                wktypestr = 'Sel Walks';
            end
            wkstr = regexprep(num2str(unique([MT.walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr = cell2mat(rgcell(unique(MT.regionnum)));
            if p.fullwalknorm
                normstr = 'full';
            else
                normstr = regexprep(num2str(normrng),'\s+','t');
            end
            transstr = regexprep(num2str(transrng),'\s+','t');
            %%%%%%%%%%%%%%%%%% Calculating Overlap %%%%%%%%%%%%%%%%%%%%%
            uPWX = unique(MT(:,{'patient','walk','dataidx'})); %need to remove chans that have the same dataidx for a given patient/walk
            [uPW,~,uPWIdx] = unique(uPWX(:,{'patient','walk'}));
            DiffSec = []; %difference between events in same patient/walk converted to sec
            for k=1:size(uPW,1)
                idx = (k==uPWIdx);
                diffsec = diff(sort(uPWX.dataidx(idx)))./250; %difference between events in same patient/walk converted to sec
                DiffSec = cat(1,DiffSec,diffsec);
            end
            OverlapWinSec = diff(transrng); %size of specgram window in sec
            OverlapIdx = DiffSec<OverlapWinSec;
            PercOverlap = sum(OverlapIdx)./numel(DiffSec)*100; %?? Need to think this through carefully
            AvgOverlapSec = mean(OverlapWinSec-DiffSec(OverlapIdx));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ttlstr = sprintf('%s, %s, %s-%s, %s-%s, %s-%s, %s \n(n=%s, p<%1.2f, %s, %s, b=%s, w=%s, o=%0.0fp-%0.0fs)',...
                transtype,desctype,regiontype,rgstr,pttypestr,ptstr,...
                wktypestr,wkstr,veltype,num2str(ntrials),...
                pval,permtype,correctiontype,normstr,transstr,PercOverlap,AvgOverlapSec);
            title(aH,ttlstr); 

            fH2 = [];
            if plottrials
                fH2 = figure('Position',[50,50,1200,800]);
                aH2 = axes('parent',fH2);
                fH2.Colormap = colormap('jet');

                [mt,sidx] = sortrows(MT,{'patient','chan','walk'});
                fidx = freq>=trialsfreqrng(1) & freq<=trialsfreqrng(2);
                fdat = squeeze(mean(npwr_trials(:,fidx,sidx),2))';

                cl = round(prctile(abs(fdat(:)),99));
                CLimit2 = [-cl,cl];
                II = fix((fdat-CLimit2(1))./(CLimit2(2)-CLimit2(1))*255)+1;
                II(II<1) = 1; II(II>256) = 256;
                iH = imagesc(tsec,1:size(fdat,1),II,'parent',aH2);
                iH.CDataMapping = "direct";

                xlim(aH2,[tsec(1),tsec(end)])
                cb2 = colorbar(aH2);
                cb2.Limits = [cb2.Limits(1),cb2.Limits(2)-1];
                cblims2 = [cb2.Limits(1),(diff(cb2.Limits)-1)/2,cb2.Limits(2)-1];
                cb2.Ticks = cblims2;
                cb2.TickLabels = [CLimit2(1),0,CLimit2(2)];
                if contains(permtype,'zscore')
                    ylabel(cb2,'Zscore')
                else
                    ylabel(cb2,'dB')
                end
                title(aH2,sprintf('Trials %0.0f-%0.0f Hz (%s)',trialsfreqrng(1),trialsfreqrng(2),ttlstr))
                [uPt,~,uPtIdx] = unique(mt.patient);
                mYtick = nan(length(uPt),1);
                hold(aH2,'on');
                for k=1:length(uPt)
                    idx = find(uPtIdx==k);
                    mYtick(k) = mean(idx);
                    plot([tsec(1)+0.3,tsec(1)+0.3],[idx(1),idx(end)],'color',[0.5,0.5,0.5],'LineWidth',10)
                end
                set(aH2,'ytick',mYtick,'yticklabel',uPt)
                ylabel(aH2,'patient')
                xlabel(aH2,'sec')
            end

            fH3 = []; fH4 = [];
            if plotboxcomp
                tidx_box1 = tsec>boxcomprng(1,1) & tsec<boxcomprng(1,2);
                tidx_box2 = tsec>boxcomprng(2,1) & tsec<boxcomprng(2,2);
                fidx_box1 = freq>boxcomprng(1,3) & freq<boxcomprng(1,4);
                fidx_box2 = freq>boxcomprng(2,3) & freq<boxcomprng(2,4);
                
                % pwr_box1 = squeeze(mean(mean(10*log10(pwr(tidx_box1,fidx_box1,:)),1),2));
                % pwr_box2 = squeeze(mean(mean(10*log10(pwr(tidx_box2,fidx_box2,:)),1),2));

                pwr_box1 = squeeze(mean(mean(pwr(tidx_box1,fidx_box1,:),1),2)); %average before converting to dB
                pwr_box2 = squeeze(mean(mean(pwr(tidx_box2,fidx_box2,:),1),2));

                t_xs = obj.MultTrans.TimeSec_xs;
                tidx_xs_box1 = t_xs>boxcomprng(1,1) & t_xs<boxcomprng(1,2);
                tidx_xs_box2 = t_xs>boxcomprng(2,1) & t_xs<boxcomprng(2,2);

                dxx = dx;
                dxx(:,isoutlier(mean(dx))) = nan;
                dx_box1 = mean(dxx(tidx_xs_box1,:))';
                dx_box2 = mean(dxx(tidx_xs_box2,:))';

                [uPtCh,~,uPtChIdx] = unique(MT(:,{'patient','regionnum'}));
                mpwr_box1 = nan(size(uPtCh,1),1);
                mpwr_box2 = nan(size(uPtCh,1),1);
                ccdx_box1 = nan(size(uPtCh,1),1);
                ccdx_box2 = nan(size(uPtCh,1),1);
                ptnum = nan(size(uPtCh,1),1);
                rgnum = nan(size(uPtCh,1),1);
                for k=1:size(uPtCh,1)
                    idx = (uPtChIdx==k);
                    ptnum(k) = uPtCh.patient(k);
                    rgnum(k) = uPtCh.regionnum(k);
                    pwr_b1 = pwr_box1(idx);
                    pwr_b2 = pwr_box2(idx);
                    dx_b1 = dx_box1(idx);
                    dx_b2 = dx_box2(idx);
                    mpwr_box1(k) = mean(pwr_b1);
                    mpwr_box2(k) = mean(pwr_b2);
                    nan_idx = isnan(dx_b1)|isnan(dx_b2);
                    cc_b1 = corrcoef(dx_b1(~nan_idx),10*log10(pwr_b1(~nan_idx))); %convert pwr to dB before correlation
                    cc_b2 = corrcoef(dx_b2(~nan_idx),10*log10(pwr_b2(~nan_idx)));
                    ccdx_box1(k) = cc_b1(1,2);
                    ccdx_box2(k) = cc_b2(1,2);
                end

                mpwr_box12 = median([pwr_box1,pwr_box2]);
                pwr_chg_box12 = (mpwr_box12(2)-mpwr_box12(1))/mpwr_box12(1); %percentage change in normalized power

                mmpwr_box12 = median([mpwr_box1,mpwr_box2]);
                mpwr_chg_box12 = (mmpwr_box12(2)-mmpwr_box12(1))/mmpwr_box12(1); %percentage change in normalized power

                fH3 = figure('Position',[50,50,1200,550]); 
                aH3 = [];
                aH3(1) = subplot(1,2,1,'parent',fH3);
                boxplot(aH3(1),10*log10([pwr_box1,pwr_box2]),{'box1 (-.)','box2 (--)'})
                [h,p,ci,stats] = ttest2(10*log10(pwr_box1),10*log10(pwr_box2));
                ylabel(aH3(1),'Power (dB)')
                title(aH3(1),sprintf('All Trials (p=%0.1e, %0.0f%%)',p,pwr_chg_box12*100))
                aH3(2) = subplot(1,2,2,'parent',fH3);
                boxplot(aH3(2),10*log10([mpwr_box1,mpwr_box2]),{'box1 (-.)','box2 (--)'})
                hold(aH3(2),'on');
                mrk_style = {'o','+','*','x'}; %by region -> 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
                mrk_size = [6,10,10,10];
                mrk_color = {'r','g','b','c','m'}; %by patient -> 1:5
                for k=1:length(mpwr_box1)
                    x1 = (rand-0.5)/3+1;
                    x2 = (rand-0.5)/3+2;
                    plot(aH3(2),x1,10*log10(mpwr_box1(k)),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
                    plot(aH3(2),x2,10*log10(mpwr_box2(k)),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);                    
                    plot(aH3(2),[x1,x2],10*log10([mpwr_box1(k),mpwr_box2(k)]),'LineStyle',':','Color',mrk_color{ptnum(k)})
                end
                [h,p,ci,stats] = ttest2(10*log10(mpwr_box1),10*log10(mpwr_box2));
                ylabel(aH3(2),'Mean Power (dB)')
                title(aH3(2),sprintf('Mean Across Patient/Region (p=%0.1e, %0.0f%%)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p,mpwr_chg_box12*100))

                % fH4 = figure('Position',[50,50,650,550]); 
                % aH4 = axes('Parent',fH4);
                % boxplot(aH4,[ccdx_box1,ccdx_box2],{'box1 (-.)','box2 (--)'})
                % hold(aH4,'on');
                % mrk_style = {'o','+','*','x'}; %by region -> 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
                % mrk_size = [6,10,10,10];
                % mrk_color = {'r','g','b','c','m'}; %by patient -> 1:5
                % for k=1:length(ccdx_box1)
                %     x1 = (rand-0.5)/3+1;
                %     x2 = (rand-0.5)/3+2;
                %     plot(aH4,x1,ccdx_box1(k),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
                %     plot(aH4,x2,ccdx_box2(k),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);                    
                %     plot(aH4,[x1,x2],[ccdx_box1(k),ccdx_box2(k)],'LineStyle',':','Color',mrk_color{ptnum(k)})
                % end
                % [h,p,ci,stats] = ttest2(ccdx_box1,ccdx_box2);
                % ylabel(aH4,'Correlation')
                % title(aH4,sprintf('Vel/Pwr Correlation Across Patient/Region (p=%0.1e)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p))

                %%%%%%%%%%%%%%%%%%% Fooof %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % setenv('KMP_DUPLICATE_LIB_OK', 'TRUE'); %this is permanently set in computer system settings
                fbins = 2:0.2:85;

                tidx_box1_np = tsec_np>boxcomprng(1,1) & tsec_np<boxcomprng(1,2);
                tidx_box2_np = tsec_np>boxcomprng(2,1) & tsec_np<boxcomprng(2,2);

                settings.max_n_peaks = 3; %aperiodic component has a better fit if more peaks are removed
                settings.aperiodic_mode = 'fixed';

                FPLabels = {'peak_freq','peak_height','peak_width','aperiodic_offset','aperiodic_exponent','fit_rsquared','fit_error'};

                FP_box1 = nan(size(np,2),length(FPLabels));
                FP_box2 = nan(size(np,2),length(FPLabels));
                % wait_msg = parfor_wait(size(np,2));
                parfor k=1:size(np,2)
                    % wait_msg.Send;
                    pwr_box1 = mean(abs(calcWavTF(np(tidx_box1_np,k),fbins,250)).^2,2);
                    pwr_box2 = mean(abs(calcWavTF(np(tidx_box2_np,k),fbins,250)).^2,2);
                    fr_box1 = fooof(fbins, pwr_box1, [fbins(1),fbins(end)], settings, 1); %fr.peak_params = [freq, height (aperiodic removed), width]
                    fr_box2 = fooof(fbins, pwr_box2, [fbins(1),fbins(end)], settings, 1); %fr.peak_params = [freq, height (aperiodic removed), width]
                    if isempty(fr_box1.peak_params)
                        pkp = nan(1,3);
                    else
                        [~,midx] = max(fr_box1.peak_params(:,2)); %find the highest peak
                        pkp = fr_box1.peak_params(midx,:);
                    end
                    FP_box1(k,:) = [pkp,fr_box1.aperiodic_params,fr_box1.r_squared,fr_box1.error];
                    if isempty(fr_box2.peak_params)
                        pkp = nan(1,3);
                    else
                        [~,midx] = max(fr_box2.peak_params(:,2)); %find the highest peak
                        pkp = fr_box2.peak_params(midx,:);
                    end
                    FP_box2(k,:) = [pkp,fr_box2.aperiodic_params,fr_box2.r_squared,fr_box2.error];
                end
                % wait_msg.Destroy;

                nFP_box12 = nan(size(FP_box1,1),2); %mean normalized across box1/box2 by patient/chan
                mFP_box12 = nan(size(uPtCh,1),2);
                for k=1:size(uPtCh,1)
                    idx = (uPtChIdx==k);
                    fp_box12 = [FP_box1(idx,5),FP_box2(idx,5)];
                    nfp_box12 = fp_box12./mean(fp_box12,2);
                    nFP_box12(idx,:) = nfp_box12;
                    mFP_box12(k,:) = mean(nfp_box12,1);
                end

                fH4 = figure('Position',[50,50,1200,550]); 
                aH4 = [];
                aH4(1) = subplot(1,2,1,'parent',fH4);
                boxplot(aH4(1),nFP_box12,{'box1 (-.)','box2 (--)'})
                [h,p,ci,stats] = ttest2(nFP_box12(:,1),nFP_box12(:,2));
                ylabel(aH4(1),'Normalized Aperiodic Exponent (Fooof)')
                title(aH4(1),sprintf('All Trials (p=%0.1e)',p))
                aH4(2) = subplot(1,2,2,'parent',fH4);
                boxplot(aH4(2),mFP_box12,{'box1 (-.)','box2 (--)'})
                hold(aH4(2),'on');
                for k=1:size(mFP_box12,1)
                    x1 = (rand-0.5)/3+1;
                    x2 = (rand-0.5)/3+2;
                    plot(aH4(2),x1,mFP_box12(k,1),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);
                    plot(aH4(2),x2,mFP_box12(k,2),'Marker',mrk_style{rgnum(k)},'MarkerFaceColor',mrk_color{ptnum(k)},'MarkerEdgeColor',mrk_color{ptnum(k)},'MarkerSize',mrk_size(rgnum(k)),'LineWidth',2);                    
                    plot(aH4(2),[x1,x2],mFP_box12(k,:),'LineStyle',':','Color',mrk_color{ptnum(k)})
                end
                [h,p,ci,stats] = ttest2(mFP_box12(:,1),mFP_box12(:,2));
                ylabel(aH4(2),'Normalized Aperiodic Exponent (Fooof)')
                title(aH4(2),sprintf('Mean Across Patient/Region (p=%0.1e)\nPatient (r=1, g=2, b=3, c=4, m=5)\nRegion (o=anthip, +=lattemp, *=entperi, x=posthip)',p))
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end

            if nargout
                varargout{1} = fH;
            end

            if nargout>1
                varargout{2} = fH2;
            end

            if nargout>2
                varargout{3} = fH3;
            end

            if nargout>3
                varargout{4} = fH4;
            end

            % rootdir = 'C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs\Combined\InOut_Transition_Specgram\StopGo';
            % fstr = sprintf('RW1-RW5_%s_%s_%s_Transition_Specgram.png',regexprep(regiontype,'\s+',''),regexprep(walktype,'\s+',''),regexprep(transtype,'\s+',''));
            % print(fH,fullfile(rootdir,fstr),'-dpng','-r300');
            % close(fH);
        end

        function varargout = plotMultTransSpecDiffPerm(obj,varargin)
            %Plot specgram difference for two transitions for all patients
            %and channels using smoothed/downsampled specgram data
            %generated in getMultTransData.
            %
            %plotMultTransSpecDiffPerm('transtype',{'Correct Turn Beg','Incorrect Turn Beg'},'regiontype',{'PostHipp+Para','PostHipp+Para'},'walktype',{'All Walks','All Walks'},'patienttype',{'All Patients','All Patients'},'permtype','zscore','correctiontype','cluster');
            %plotMultTransSpecDiffPerm('transtype',{'Doorway','Doorway'},'desctype',{'open','closed'},'regiontype',{'PostHipp+Para','PostHipp+Para'},'walktype',{'All Walks','All Walks'},'patienttype',{'All Patients','All Patients'},'permtype','zscore','correctiontype','cluster');
            p = obj.parseInputs(varargin{:});

            if isempty(p.permtype)
                error('permtype must be specified!');
            end

            warning('off','MATLAB:contour:ConstantData');

            transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            customregion = p.customregion;
            if isempty(customregion)
                customregion = {[],[]};
            end
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,3'
            veltype = p.veltype;
            if isempty(veltype)
                veltype = {'',''};
            end
            patienttype = p.patienttype; %'All Patients', '1,3'
            desctype = p.desctype; %close, open, etc. (optional, will skip if empty)
            if isempty(desctype)
                desctype = {'',''};
            end
            permtype = p.permtype; %standard, zscore
            correctiontype = p.correctiontype; %cluster, pixel (can be empty)
            if strcmp(correctiontype,'fdr+cluster') && strcmp(permtype,'standard')
                error('fdr+cluster cannot be performed with the standard permtype!');
            end
            transrng = p.transrng; %default [-10,10]
            normrng = transrng;
            pval = p.pval; %default p=0.05
            pvalclust = p.pvalclust; %default 0.01
            if isempty(p.clim)
                climit = [-10,10];
            else
                climit = p.clim;
            end

            %Init some params
            tsamp = round(transrng(1)*25):round(transrng(2)*25); %+-10sec at 25Hz (downsampled by 10 from 250)
            nperm = 1000; %permutations

            %Getting trials and specgram data
            [MT1,~,dx1] = obj.filterMultTransData('transtype',transtype{1},'regiontype',regiontype{1},'walktype',walktype{1},'desctype',desctype{1},'patienttype',patienttype{1},'veltype',veltype{1},'customregion',customregion{1}); %table of all trials and corresponding info
            [pwr1,tsec,freq] = obj.getFilteredMultTransData(MT1,tsamp); %raw power for all trials (time x freq x trial)

            [MT2,~,dx2] = obj.filterMultTransData('transtype',transtype{2},'regiontype',regiontype{2},'walktype',walktype{2},'desctype',desctype{2},'patienttype',patienttype{2},'veltype',veltype{2},'customregion',customregion{2}); %table of all trials and corresponding info
            pwr2 = obj.getFilteredMultTransData(MT2,tsamp); %raw power for all trials (time x freq x trial)

            nfreq = length(freq);
            ntime = length(tsamp);
            ntrials1 = size(MT1,1);
            ntrials2 = size(MT2,1);

            trialtype = [-ones(1,ntrials1),ones(1,ntrials2)];
            pwr = cat(3,pwr1,pwr2);
            if ~p.fullwalknorm %Normalization is done across the entire walk for each chan in getMultTransData if true
                pwr = pwr./mean(pwr,1,"omitnan"); %full epoch norm (time x freq x trial) %Now doing this across the entire walk for each chan in getMultTransData
            end
            pwr = log10(pwr); %log10 normalize to make a normal distribution for welch's ttest below!!

            %Real (Do unequal trials violate ttest assumptions? Not for Welch's ttest below.)
            tnum = squeeze(mean(pwr(:,:,trialtype==-1),3,"omitnan")-mean(pwr(:,:,trialtype==1),3,"omitnan")); %trial avg diff (time x freq)
            tdenom = sqrt(std(pwr(:,:,trialtype==-1),0,3,"omitnan").^2./ntrials1+std(pwr(:,:,trialtype==1),0,3,"omitnan").^2./ntrials2);
            npwr = tnum./tdenom;

            %Permutations
            disp('Running permutations...')
            % PM = calcRWAPerm(pwr,nan(1,nfreq),trialtype);
            PM = calcRWAPerm_mex(pwr,nan(1,nfreq),trialtype); %2.5 times faster

            %Finding percentiles of permuted distributions
            pm_freq = reshape(permute(PM,[1,3,2]),[],nfreq); %(ntime*nperm) x nfreq
            pc = prctile(pm_freq,[pval*100/2,100-pval*100/2],1); %2 x nfreq
            mPM = mean(pm_freq,1); %mean across permutations and time (1 x nfreq)
            sPM = std(pm_freq,0,1); %std across permutations and time (1 x nfreq)
         
            %Finding clusters in permuted data
            max_clust_info = nan(nperm,1);
            max_pixel_pvals = zeros(nperm,2);
            for m=1:nperm
                pm = PM(:,:,m);
                switch permtype
                    case 'standard'
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq -> maybe only pool across time?)
                        pm(pm>squeeze(pc(1,:,:)) & pm<squeeze(pc(2,:,:))) = 0;
                    case 'zscore'
                        pm = (pm-mPM)./sPM;
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq -> maybe only pool across time?)
                        switch correctiontype
                            case 'fdr+cluster'
                                PV = 1-normcdf(abs(pm));
                                PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                                pm(PV>=pval) = 0;
                            otherwise
                                pm(abs(pm)<norminv(1-pval)) = 0;
                        end
                end
                clustinfo = bwconncomp(pm);
                stats = regionprops(clustinfo,pm,'pixelvalues');
                max_clust_info(m) = max([0,abs(cellfun(@sum,{stats.PixelValues}))]); %max cluster size using pixel values for each permutation
%                 max_clust_info(m) = max([0,cellfun(@numel,clustinfo.PixelIdxList)]); %max cluster size for each permutation
            end

            %Finding percentiles of pixel-level distributions
            pc_pixel = prctile(max_pixel_pvals(:),[pval*100/2,100-pval*100/2]); %two tail

            %Thresholding real
            switch permtype
                case 'standard'
                    npwr_thresh = npwr;
                    switch correctiontype
                        case 'fdr' %results are slightly different than PermSpecGram because this operates on log10 data and does rotating shuffle
                            %Finding pvals based on rank in distribution for FDR correction
                            % PV = calcRWAPVal(npwr,pm_freq);
                            PV = calcRWAPVal_mex(npwr,pm_freq); %~3x faster
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(npwr));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(npwr_thresh>squeeze(pc(1,:,:)) & npwr_thresh<squeeze(pc(2,:,:))) = 0;
                    end
                case 'zscore'
                    npwr = (npwr-mPM)./sPM;
                    npwr_thresh = npwr;
                    switch correctiontype
                        case {'fdr','fdr+cluster'}
                            PV = 1-normcdf(abs(npwr));
                            PV = reshape(mafdr(PV(:),'BHFDR',true),size(PV));
                            npwr_thresh(PV>=pval) = 0;
                        case 'pixel'
                            npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                        otherwise %cluster or empty
                            npwr_thresh(abs(npwr)<norminv(1-pval))=0; %norminv gives stat value at pval for normal distribution (i.e. p=0.05 is -1.65)
                    end
            end

            %Removing clusters
            if contains(correctiontype,'cluster') && ~isempty(correctiontype)
                clustinfo = bwconncomp(npwr_thresh);
                stats = regionprops(clustinfo,npwr_thresh,'pixelvalues');
                clust_info = abs(cellfun(@sum,{stats.PixelValues}));
%                 clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                clust_threshold = prctile(max_clust_info,100-pvalclust*100);
                whichclusters2remove = find(clust_info<clust_threshold);
                for k=1:length(whichclusters2remove)
                    npwr_thresh(clustinfo.PixelIdxList{whichclusters2remove(k)})=0;
                end
            end

            %Plotting
            fH = figure('Position',[50,50,1200,800]);
            aH = axes('parent',fH);
            colormap("jet");
            contourf(tsec,freq,npwr',100,'linecolor','none','parent',aH);
            hold(aH,"on");
            contour(tsec,freq,(npwr_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
            set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'yminortick','off','clim',climit); %now in units of standard deviation
            plot(aH,[0,0],[2,120],'--k','LineWidth',2);
            dxx = mean(dx1(:,~isoutlier(mean(dx1))),2,"omitnan")-mean(dx2(:,~isoutlier(mean(dx2))),2,"omitnan");
            dxx = dxx - mean(dxx);
            dxx_me = mean(dxx,2,"omitnan");
            dxx_se = std(dxx,0,2,'omitnan');
            if p.includevel
                plot(aH,obj.MultTrans.TimeSec_xs,(dxx_me)*9+10,'k','LineWidth',2);
                plot(aH,obj.MultTrans.TimeSec_xs,((dxx_me-dxx_se))*9+10,':k','LineWidth',0.5);
                plot(aH,obj.MultTrans.TimeSec_xs,((dxx_me+dxx_se))*9+10,':k','LineWidth',0.5);
            end
            cb = colorbar(aH);
            cblims = [cb.Limits(1),0,cb.Limits(2)];
            cb.Ticks = cblims;
            cb.TickLabels = cblims;
            xlim(aH,[tsec(1),tsec(end)])
            ylim(aH,[2,120])
            if contains(permtype,'zscore')
                ylabel(cb,'Zscore')
            else
                ylabel(cb,'Tscore')
            end
            xlabel(aH,'sec');
            ylabel(aH,'Hz');
            if isempty(str2num(patienttype{1}))
                pttypestr1 = patienttype{1};
            else
                pttypestr1 = 'Sel Patients';
            end
            if isempty(str2num(patienttype{2}))
                pttypestr2 = patienttype{2};
            else
                pttypestr2 = 'Sel Patients';
            end
            ptstr1 = regexprep(num2str(unique([MT1.patient])'),'\s+','');
            ptstr2 = regexprep(num2str(unique([MT2.patient])'),'\s+','');
            if isempty(str2num(walktype{1}))
                wktypestr1 = walktype{1};
            else
                wktypestr1 = 'Sel Walks';
            end
            if isempty(str2num(walktype{2}))
                wktypestr2 = walktype{2};
            else
                wktypestr2 = 'Sel Walks';
            end
            wkstr1 = regexprep(num2str(unique([MT1.walk])'),'\s+','');
            wkstr2 = regexprep(num2str(unique([MT2.walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr1 = cell2mat(rgcell(unique(MT1.regionnum)));
            rgstr2 = cell2mat(rgcell(unique(MT2.regionnum)));
            if p.fullwalknorm
                normstr = 'full';
            else
                normstr = regexprep(num2str(normrng),'\s+','t');
            end
            transstr = regexprep(num2str(transrng),'\s+','t');
            title(aH,sprintf('%s, %s, %s, %s, %s, %s \n(n=%s, p<%1.2f, %s, %s, b=%s, w=%s)',...
                [transtype{1},'::',transtype{2}],...
                [desctype{1},'::',desctype{2}],...
                [regiontype{1},'-',rgstr1,'::',regiontype{2},'-',rgstr2],...
                [pttypestr1,'-',ptstr1,'::',pttypestr2,'-',ptstr2],...
                [wktypestr1,'-',wkstr1,'::',wktypestr2,'-',wkstr2],...
                [veltype{1},'::',veltype{2}],...
                [num2str(ntrials1),'::',num2str(ntrials2)],...
                pval,permtype,correctiontype,normstr,transstr));

            if nargout
                varargout{1} = fH;
            end

            % rootdir = 'C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs\Combined\InOut_Transition_Specgram\StopGo';
            % fstr = sprintf('RW1-RW5_%s_%s_%s_Transition_Specgram.png',regexprep(regiontype,'\s+',''),regexprep(walktype,'\s+',''),regexprep(transtype,'\s+',''));
            % print(fH,fullfile(rootdir,fstr),'-dpng','-r300');
            % close(fH);

        end

        function plotMultTransSpecGramPermBL(obj,varargin)
            %Plot specgram of transitions for all patients and channels
            %using smoothed/downsampled specgram data generated in
            %getMultTransData. Uses randomly sampled times during a walk
            %for the baseline. Run dispCategories to see options for
            %transtype, regiontype, walktype, desctype, permtype, and
            %correctiontype.
            %
            %plotMultTransSpecGramPerm('transtype','Outdoor Beg','regiontype','PostHipp+Para','walktype','All Walks','permtype','standard','correctiontype','cluster');
            %plotMultTransSpecGramPerm('transtype','Doorway','regiontype','PostHipp+Para','walktype','All Walks','desctype','closed','permtype','standard','correctiontype','cluster');
            %plotMultTransSpecGramPerm('transtype','Doorway','desctype','in2out&closed','regiontype','AntHipp','walktype','All Walks','permtype','zscore','correctiontype','cluster');
            p = obj.parseInputs(varargin{:});

            if isempty(p.permtype)
                error('permtype must be specified!');
            end

            warning('off','MATLAB:contour:ConstantData');

            transtype = p.transtype; %'Outdoor Beg', 'Outdoor End', 'Doorway'
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', [individual walks]
            desctype = p.desctype; %close, open, etc. (optional, will skip if empty)
            permtype = p.permtype; %standard, zscore
            correctiontype = p.correctiontype; %cluster, pixel

            %Some hard-coded params
            tsamp = -250:250; %+-10sec at 25Hz (downsampled by 10 from 250)
            nperm = 1000; %permutations
            pval = 0.05;

            %Getting trials and specgram data
            [MT,~,dx] = obj.filterMultTransData(varargin{:}); %table of all trials and corresponding info
            [pwr,tsec,freq] = obj.getFilteredMultTransData(MT,tsamp); %raw power for all trials (time x freq x trial)

            nfreq = length(freq);
            ntime = length(tsamp);
            ntrials = size(MT,1);

            %Real
            pwr = pwr./mean(pwr,1,"omitnan"); %full epoch norm (time x freq x trial)
            mpwr = squeeze(mean(pwr,3,"omitnan")); %trial avg (time x freq)

            %Permutations (randomly select windows across entire walk)
            disp('Running permutations...')
            rng('shuffle'); %seed the random stream with clock time
            [uPtWkCh,~,uPtWkChIdx] = unique(MT(:,{'patient','walk','chan'}));
            MWT = obj.MultTrans.MultWTable;
            PM = nan(ntime,nfreq,nperm);
            wait_msg = parfor_wait(nperm);
            parfor (m=1:nperm,6)
                wait_msg.Send;
                pm = nan(size(uPtWkChIdx,1),ntime,nfreq);
                for k=1:size(uPtWkCh,1)
                    idx = find(k==uPtWkChIdx);
                    mt = MT(idx,:);
                    pt = uPtWkCh.patient(k);
                    wk = uPtWkCh.walk(k);
                    ch = uPtWkCh.chan(k);
                    dd = MWT{pt}{wk}(:,:,ch); %time x freq (25Hz downsampled wavelet data for specified patient/walk/chan)
                    nt = size(mt,1); %trials
                    ns = size(dd,1); %samples to choose from
                    r = randsample(tsamp(end)+100:ns-tsamp(end)-100,nt)'; %avoid edges with 4sec border (100 samples)
                    didx = r+tsamp; %trial x time
                    dd = reshape(dd(didx,:),size(didx,1),size(didx,2),nfreq); %trial x time x freq
                    cutpoint = randi(size(dd,2),[1,nt]);
                    for n=1:nt %double randomization (seems to eliminate the edge effect)
                        dd(n,:,:) = circshift(dd(n,:,:),cutpoint(n),2); %trial x time x freq (shift time)
                    end
                    pm(idx,:,:) = dd(1:nt,:,:); %trial x time x freq
                end
                pm = pm./mean(pm,2,"omitnan"); %full epoch norm (trial x time x freq)
                pm = squeeze(mean(pm,1,"omitnan")); %trial avg
                PM(:,:,m) = pm;
            end
            wait_msg.Destroy;

            mBase = mean(PM,3,"omitnan"); %mean across permutations for baseline

            %Using mean of all permutations as baseline
            npwr = 10*log10(mpwr./mBase); %normalized power in dB
            PM = 10*log10(PM./mBase);

            mPM = mean(PM,3,"omitnan"); %mean across permutations
            sPM = std(PM,0,3,"omitnan"); %std across permutations

            %Finding percentiles of permuted distributions
            %Keep time (similar to the zscore method)
            pm = reshape(PM,[],nperm); %finding percentiles for all time x freq distributions (instead of collapsing across time like above)
            pc_low = reshape(prctile(pm,pval*100/2,2),ntime,nfreq); %two tail (across 2nd dim -> 1000 perms)
            pc_high = reshape(prctile(pm,100-pval*100/2,2),ntime,nfreq); %two tail
            pc = permute(cat(3,pc_low,pc_high),[3,1,2]);

            %Finding clusters in permuted data
            max_clust_info = nan(nperm,1);
            max_pixel_pvals = zeros(nperm,2);
            for m=1:nperm
                pm = PM(:,:,m);
                switch permtype
                    case 'standard'
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        pm(pm>squeeze(pc(1,:,:)) & pm<squeeze(pc(2,:,:))) = 0;
                    case 'zscore'
                        pm = (pm-mPM)./sPM;
                        max_pixel_pvals(m,:) = [min(pm(:)),max(pm(:))]; %pixel correction distributions (pooled across time and freq)
                        pm(abs(pm)<norminv(1-pval)) = 0;
                end
                clustinfo = bwconncomp(pm);
                max_clust_info(m) = max([0,cellfun(@numel,clustinfo.PixelIdxList)]); %max cluster size for each permutation
            end

            %Finding percentiles of pixel-level distributions
            pc_pixel = prctile(max_pixel_pvals(:),[pval*100/2,100-pval*100/2]); %two tail

            %Thresholding real
            switch permtype
                case 'standard'
                    npwr_thresh = npwr;
                    if contains(correctiontype,'cluster')
                        npwr_thresh(npwr_thresh>squeeze(pc(1,:,:)) & npwr_thresh<squeeze(pc(2,:,:))) = 0;
                    else
                        npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                    end
                case 'zscore'
                    npwr = (npwr-mPM)./sPM;
                    npwr_thresh = npwr;
                    if contains(correctiontype,'cluster')
                        npwr_thresh(abs(npwr)<norminv(1-pval))=0; %norminv gives stat value at pval for normal distribution (i.e. p=0.05 is -1.65)
                    else
                        npwr_thresh(npwr_thresh>pc_pixel(1) & npwr_thresh<pc_pixel(2)) = 0;
                    end
            end

            %Removing clusters
            if contains(correctiontype,'cluster')
                clustinfo = bwconncomp(npwr_thresh);
                clust_info = cellfun(@numel,clustinfo.PixelIdxList);
                clust_threshold = prctile(max_clust_info,100-pval*100);
                whichclusters2remove = find(clust_info<clust_threshold);
                for k=1:length(whichclusters2remove)
                    npwr_thresh(clustinfo.PixelIdxList{whichclusters2remove(k)})=0;
                end
            end

            %Plotting
            fH = figure('Position',[50,50,1200,800]);
            aH = axes('parent',fH);
            colormap("jet");
            contourf(tsec,freq,npwr',40,'linecolor','none','parent',aH);
            hold(aH,"on");
            contour(tsec,freq,(npwr_thresh~=0)',1,'parent',aH,'linecolor','w','linewidth',2);
            if contains(permtype,'zscore')
                set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'clim',[-10,10]); %now in units of standard deviation
            else
                set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'clim',[-0.5,0.5]);
            end
            plot(aH,[0,0],[2,64],'--k','LineWidth',2);
            dxx = mean(dx(:,~isoutlier(mean(dx))),2,"omitnan");
            plot(aH,obj.MultTrans.TimeSec_xs,(dxx-median(dxx))*9+10,'k','LineWidth',2);
            cb = colorbar(aH);
            cblims = [cb.Limits(1),0,cb.Limits(2)];
            cb.Ticks = cblims;
            cb.TickLabels = cblims;
            xlim(aH,[tsec(1),tsec(end)])
            ylim(aH,[2,64])
            if contains(permtype,'zscore')
                ylabel(cb,'Zscore')
            else
                ylabel(cb,'dB')
            end
            xlabel(aH,'sec');
            ylabel(aH,'Hz');
            ptstr = regexprep(num2str(unique([MT.patient])'),'\s+','/');
            title(aH,sprintf('%s, %s, %s, %s, %s, %s, p<%1.2f',...
                transtype,desctype,regiontype,ptstr,...
                num2str(walktype),num2str(ntrials),pval));

            % rootdir = 'C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs\Combined\InOut_Transition_Specgram\StopGo';
            % fstr = sprintf('RW1-RW5_%s_%s_%s_Transition_Specgram.png',regexprep(regiontype,'\s+',''),regexprep(walktype,'\s+',''),regexprep(transtype,'\s+',''));
            % print(fH,fullfile(rootdir,fstr),'-dpng','-r300');
            % close(fH);

        end
        
        function plotTransSpecGramByWalk(obj,varargin)
            %plots specgrams for all walks for a specified event and channel
            %
            %plotTransSpecGramByWalk('patient',1,'chan',1,'evntnum',1,'transtype','Outdoor Beg');

            p = parseInputs(obj,varargin{:});
            if isempty(p.patient)
                error('Need to specify patient!');
            end
            if isempty(p.chan)
                error('Need to specify chan!');
            end
            if isempty(p.evntnum)
                error('Need to specify evntnum!');
            end
            if isempty(p.transtype)
                error('Need to specify transtype!');
            end
            pt = p.patient;
            chan = p.chan;
            evnt = p.evntnum;
            trans = p.transtype;

            if all(obj.PatientIdx~=pt)
                obj.loadData('PatientIdx',pt);
            end
           
            obj.loadVid; %loads all videos for a patient and saves to obj.Vid

            db = obj.DB;
            db = db(contains(db.Event,trans),:);
            walk_list = unique(db.Walk)'; walk_list(end) = []; %last walk is in reverse
            frame_list = nan(1,length(walk_list));
            np_samp = nan(1,length(walk_list));
            for k=1:length(walk_list)
                idx = find(db.Walk==walk_list(k));
                frame_list(k) = db.PupilFrame(idx(evnt));
                np_samp(k) = db.NPSample(idx(evnt));
            end

            winsec = 32;

            fs_np = obj.DTable.fs_np(1);
            nsamp = winsec*fs_np;
            tsamp = (-nsamp:nsamp); %+-32 sec (2sec removed for artifact, so +-30sec)
            tsec = tsamp/fs_np;

            fH = figure('Position',[50,50,1600,1000],'Visible','on');
            aH = nan(1,length(walk_list)); vH = nan(1,length(walk_list));
            for k=1:length(walk_list)
                aH(k) = subplot(length(walk_list),10,(10*k-9:10*k-1),'parent',fH);

                wlk = walk_list(k);
                frm = frame_list(k);
                tfull = np_samp(k)+tsamp;
                d = obj.DTable.d_np{wlk}(tfull,chan);
                d(d<-500) = 0;
                d(isnan(d)) = 0;
                lb = obj.ChanLabels{pt,chan};
                PSG = PermSpecGram(d,'Fs',250,'FPass',[2,64],'TimeRng',[tsec(1),tsec(end)],'NormRng',[tsec(1)+2,tsec(end)-2]);
                PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','Smoothing',[500,1],'ErrPerc',[]);
                PSG.plotSpecGram('ShowRaw',true,'Clim',[-5,5],'aH',aH(k));
                plot(aH(k),[0,0],[0,8],'k'); %center line
                xlim(aH(k),[tsec(1)+2,tsec(end)-2])
%                 xlim(aH(k),[-20,20])
                title(aH(k),sprintf('Walk %0.0f Chan %0.0f (%s)',wlk,chan,lb))

                vH(k) = subplot(length(walk_list),10,10*k,'parent',fH);
                
                image(read(obj.Vid.VidObj{obj.Vid.Walk==wlk},frm),'parent',vH(k));
                axis(vH(k),'image')
                set(vH(k),'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
            end


        end

        function plotWalkSpecGram(obj,varargin)
            %Plots continuous raw spectrogram for an entire walk and
            %overlays events and other variables like speed. Will save to
            %png if specified. Uses MultWTable.
            %
            %plotWalkSpecGram('patient',1,'walknum',1,'chan',1,'saveflag',true);
            p = parseInputs(obj,varargin{:});
            if isempty(p.patient)
                error('Need to specify patient!');
            end
            if isempty(p.walknum)
                error('Need to specify walknum!');
            end
            if isempty(p.chan)
                error('Need to specify chan!');
            end
            if isempty(obj.MultTrans)
                obj.getMultTransData;
            end
            if ~any(obj.PatientIdx==p.patient)
                obj.loadData('PatientIdx',p.patient);
            end
            if isempty(obj.InOutTimes)
                obj.parseSegData;
            end

            patient = p.patient;
            walknum = p.walknum;
            chan = p.chan;

            dbb = obj.DB(obj.DB.Walk==walknum,:);
            % idx = contains(dbb.Event,'Doorway') & contains(dbb.Description,'[closed]');
            idx = contains(dbb.Event,'Doorway');
            ts_evnt = dbb.NPSample(idx)/250;

            t_np = obj.DTable.ntp_np{walknum}; t_np = t_np-t_np(1); %ntp_np (sec)
            d_np = obj.DTable.d_np{walknum}(:,chan); %contains nan for outlier regions (overlay to identify outliers in specgram)

            T = obj.InOutTimes(obj.InOutTimes.Walk==walknum,:); %these are np samples not time in sec
            T.Start_NP = T.Start_NP./250; %convert to sec from start of np recording
            T.End_NP = T.End_NP./250;

            % fs = 250;
            % D(isnan(D)) = 0;
            % D(D<-500) = 0;
            % D(dbb.NPSample(end)+1:end) = [];
            % tsec = (0:length(D)-1)./fs; %this is slightly different than t above (usually some small jitter less than a sample in time)
            % PSG = PermSpecGram(D,'Fs',fs,'FPass',[2,64],'TimeRng',[0,length(D)-1]./fs);
            % PSG.calcSpecGram('AnalysisType','Pwr','NormType','None','Smoothing',[5000,5]);
            % 
            % S = (PSG.S);
            % f = PSG.fSpec;

            %calculating power with the hilbert instead
%             [B,A] = butter(3,[5,8]/(fs/2),'bandpass');
%             Dp = filtfilt(B,A,D);
%             Dp = hilbert(Dp); %bipolar seeg
%             Dp = abs(Dp).^2;
%             Dp = smoothdata(Dp,'movmean',1000);
%             Dp = Dp - mean(Dp);
% 
%             figure;
%             plot(tsec,Dp)
%             hold on;
%             if ~isempty(ts)
%                 plot([ts,ts],[-200,1000],'k')
%             end

            %Outliers are set to zero for specgram calculation in
            %getMultTransData (maybe interpolate instead?)
            fs_spec = 25; %specgrams are smoothed/downsampled to 25Hz
            S = 10*log10(obj.MultTrans.MultWTable{patient}{walknum}(:,:,chan)); %full epoch norm already done in getMultTransData
            f_spec = obj.MultTrans.Freq;
            t_spec = (0:size(S,1)-1)/fs_spec; %sec

            %plot specgram for full walk with all events and inside/outside
            %marked
            if p.saveflag
                fH = figure('Position',[50,50,1800,850],'Visible','off','Colormap',parula);
            else
                fH = figure('Position',[50,50,1800,850],'Visible','on','Colormap',parula);
            end
            aH = subplot(20,1,1:9,'parent',fH);
            % imagesc(t_spec,log2(f_spec),S','parent',aH);
            contourf(t_spec,f_spec,S',100,'linecolor','none','parent',aH);
            hold(aH,'on')
            % plot(t_np,d_np/25+2.^3.5,'r') %raw data

            cc = bwconncomp(isnan(d_np));
            B = regionprops(cc,'BoundingBox');
            B = cat(1,B.BoundingBox);
            for m=1:size(B,1)
                b = B(m,[2,4]);
                x = floor([b(1),sum(b),sum(b),b(1)]);
                x(x<1) = 1; x(x>length(t_np)) = length(t_np);
                y = [2,2,64,64];
                p4 = patch(t_np(x),y,'r','FaceAlpha',1,'EdgeColor','none','parent',aH);
            end

            %transitions
            if isempty(ts_evnt)
                p1 = plot([nan,nan],[2,64],'k','LineWidth',2);
            else
                p1 = plot([ts_evnt,ts_evnt]',repmat([2,64],length(ts_evnt),1)','k','LineWidth',2);
            end

            idx = T.Location=="Indoor";
            p2 = plot([T.Start_NP(idx)+3,T.End_NP(idx)-3]',[2.2,2.2],'color',[0.5,0.5,0.5],'linewidth',10); %indoor
            p3 = plot([T.Start_NP(~idx)+3,T.End_NP(~idx)-3]',[2.2,2.2],'color',[0,1,0],'linewidth',10); %outdoor

            set(aH,'xtick',(0:2:20)*60,'xticklabel',0:2:20,'xlim',[T.Start_NP(1)-6,T.End_NP(end)+6]);
            set(aH,'yscale','log','YTick',2.^(1:6),'yticklabel',2.^(1:6),'yminortick','off','ylim',[2,64]); %now in units of standard deviation
            set(aH,'clim',[-5,5])

            xlabel(aH,'min')
            ylabel(aH,'Hz')
            title(sprintf('%s, Walk%0.0f, Chan%0.0f, %s',obj.PatientID,walknum,chan,obj.ChanLabels{obj.PatientIdx,chan}),"FontSize",14)
           
            legend([p1(1),p2(1),p3(1),p4(1)],{'Doorway','Inside','Outside','Dropouts'})

            aH2 = subplot(20,1,12:20,'parent',fH);
            hold(aH2,"on");

            %plot velocity from xsens
            ntp_np = obj.DTable.ntp_np{walknum};

            % sm_x = 2; %2sec smoothing
            sm_x = 4; %4sec smoothing

            fs_xs = 100;
            d_xs = smoothdata(obj.DTable.d_xs{walknum},'movmean',fs_xs*sm_x)./2; d_xs(d_xs>5) = nan; d_xs = d_xs-nanmedian(d_xs); 
            t_xs = obj.DTable.ntp_xs{walknum}; t_xs = t_xs - ntp_np(1);
            p_xs = plot(t_xs,d_xs,'b');
            p2_xs = plot(t_xs(1:fs_xs*sm_x:end),d_xs(1:fs_xs*sm_x:end),'.b');

            fs_gz = 200;
            d_gz = smoothdata(obj.DTable.d_gaze_fix{walknum},'movmean',fs_gz*sm_x)./1.5; d_gz = d_gz-nanmedian(d_gz); %this needs a bit more smoothing to look good in plot
            t_gz = obj.DTable.ntp_gaze{walknum}; t_gz = t_gz - ntp_np(1);
            p_gz = plot(t_gz,d_gz-1,'r');
            p2_gz = plot(t_gz(1:fs_gz*sm_x:end),d_gz(1:fs_gz*sm_x:end)-1,'.r');

            fs_kd = 60;
            d_kd = smoothdata(obj.DTable.d_kde{walknum},'movmean',fs_kd*sm_x)./12; d_kd = d_kd-nanmedian(d_kd);
            t_kd = obj.DTable.ntp_kde{walknum}; t_kd = t_kd - ntp_np(1);
            p_kd = plot(t_kd,d_kd-2,'g');
            p2_kd = plot(t_kd(1:fs_kd*sm_x:end),d_kd(1:fs_kd*sm_x:end)-2,'.g');

            fs_am = 15;
            d_am = log10(smoothdata(obj.DTable.d_amb{walknum},'movmean',fs_am*sm_x)+eps); d_am = d_am./max(d_am); d_am = d_am-nanmedian(d_am);
            t_am = obj.DTable.ntp_amb{walknum}; t_am = t_am - ntp_np(1);
            p_am = plot(t_am,d_am-3,'m');
            p2_am = plot(t_am(1:fs_am*sm_x:end),d_am(1:fs_am*sm_x:end)-3,'.m');
          
            set(aH2,'xtick',(0:2:20)*60,'xticklabel',0:2:20,'xlim',[T.Start_NP(1)-6,T.End_NP(end)+6]);
            set(aH2,'ytick',-3:0,'yticklabel',{'amb','kde','fix','vel'},'ylim',[-3.5,0.5])
            xlabel(aH2,'min')

            linkaxes([aH,aH2],'x')

            if p.saveflag
                fstr = sprintf('%s_Walk%0.0f_Chan%0.0f_WalkSpecGram.png',obj.PatientID,walknum,chan);
                print(fH,fullfile('D:\RWN\WalkSpecGrams',fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function batchWalkSpecGram(obj,varargin)
            %batchWalkSpecGram('patient',1);
            p = parseInputs(obj,varargin{:});
            if isempty(p.patient)
                error('Need to specify patient!');
            end
            obj.loadData('PatientIdx',p.patient);
            obj.parseSegData;
            walknums = unique(obj.DB.Walk);
            for k=1:4 %chan
                for m=1:length(walknums) %walk
                    try
                        obj.plotWalkSpecGram('patient',p.patient,'chan',k,'walknum',walknums(m),'saveflag',true);
                        fprintf('Success: chan%0.0f, walk%0.0f\n',k,walknums(m))
                    catch
                        fprintf('Error: chan%0.0f, walk%0.0f\n',k,walknums(m))
                    end
                end
            end
        end

        function plotFeatConfusionMatrix(obj,varargin)
            %Plots continuous raw spectrogram for an entire walk and
            %overlays events and other variables like speed. Will save to
            %png if specified. Uses MultWTable.
            %
            %plotWalkSpecGram('patient',1,'walknum',1,'chan',1,'saveflag',true);
            p = parseInputs(obj,varargin{:});
            if isempty(p.patient)
                error('Need to specify patient!');
            end
            if isempty(p.walknum)
                error('Need to specify walknum!');
            end
            if isempty(p.chan)
                error('Need to specify chan!');
            end
            if isempty(obj.MultTrans)
                obj.getMultTransData;
            end
            if ~any(obj.PatientIdx==p.patient)
                obj.loadData('PatientIdx',p.patient);
            end
            if isempty(obj.InOutTimes)
                obj.parseSegData;
            end

            patient = p.patient;
            walknum = p.walknum;
            chan = p.chan;

            dbb = obj.DB(obj.DB.Walk==walknum,:);
            ntp_wk_beg = dbb.NTP(contains(dbb.Event,'Walk Beg'));
            ntp_wk_end = dbb.NTP(contains(dbb.Event,'Walk End'));

            ntp_np = obj.DTable.ntp_np{walknum}; 
            ntp_spec = ntp_np(1:10:end);
            
            %Outliers are set to zero for specgram calculation in
            %getMultTransData (maybe interpolate instead?)
            fs_spec = 25; %specgrams are smoothed/downsampled to 25Hz
            S = obj.MultTrans.MultWTable{patient}{walknum}(:,:,chan); %full epoch norm already done in getMultTransData
            f_spec = obj.MultTrans.Freq;

            if length(ntp_spec)~=size(S,1)
                error('Downsampled ntp times and length of specgram data do not match!');
            end

            fidx = false(length(f_spec),5);
            fidx(:,1) = (f_spec>=4 & f_spec<8); %theta
            fidx(:,2) = (f_spec>=8 & f_spec<12); %alpha
            fidx(:,3) = (f_spec>=12 & f_spec<30); %beta
            fidx(:,4) = (f_spec>=30 & f_spec<60); %gamma
            fidx(:,5) = (f_spec>=60 & f_spec<85); %high gamma

            [~,tidx1] = min(abs(ntp_spec-ntp_wk_beg));
            [~,tidx2] = min(abs(ntp_spec-ntp_wk_end));
            tidx = tidx1:tidx2;

            ss = nan(length(tidx),size(fidx,2));
            for k=1:size(fidx,2)
                s = 10*log10(mean(S(tidx,fidx(:,k)),2)); 
                s(isoutlier(s)) = nan; 
                nanidx = isnan(s); 
                nanval = interp1(tidx(~nanidx),s(~nanidx),tidx(nanidx),'linear','extrap'); 
                s(nanidx) = nanval;
                ss(:,k) = s;
            end

            sm_x = obj.MultTrans.WinSegSec; %2/4sec smoothing

            %xsens
            fs_xs = 100;
            d_xs = smoothdata(obj.DTable.d_xs{walknum},'movmean',fs_xs*sm_x); d_xs(d_xs>5) = nan; 
            ntp_xs = obj.DTable.ntp_xs{walknum}; 
            [~,tidx1] = min(abs(ntp_xs-ntp_wk_beg));
            tidx = tidx1:length(d_xs);
            d_xs = d_xs(tidx);
            nanidx = isnan(d_xs);
            nanval = interp1(tidx(~nanidx),d_xs(~nanidx),tidx(nanidx),'linear','extrap');
            d_xs(nanidx) = nanval;
            d_xs = resample(d_xs,fs_spec,fs_xs);
            d_xs = d_xs(1:size(ss,1));
           
            %gaze
            fs_gz = 200;
            d_gz = smoothdata(obj.DTable.d_gaze_fix{walknum},'movmean',fs_gz*sm_x); 
            ntp_gz = obj.DTable.ntp_gaze{walknum}; 
            [~,tidx1] = min(abs(ntp_gz-ntp_wk_beg));
            tidx = tidx1:length(d_gz);
            d_gz = d_gz(tidx);
            nanidx = isnan(d_gz);
            nanval = interp1(tidx(~nanidx),d_gz(~nanidx),tidx(nanidx),'linear','extrap');
            d_gz(nanidx) = nanval;
            d_gz = resample(d_gz,fs_spec,fs_gz);
            d_gz = d_gz(1:size(ss,1));
           
            %kde
            fs_kd = 60;
            d_kd = smoothdata(obj.DTable.d_kde{walknum},'movmean',fs_kd*sm_x); 
            ntp_kd = obj.DTable.ntp_kde{walknum}; 
            [~,tidx1] = min(abs(ntp_kd-ntp_wk_beg));
            tidx = tidx1:length(d_kd);
            d_kd = d_kd(tidx);
            nanidx = isnan(d_kd);
            nanval = interp1(tidx(~nanidx),d_kd(~nanidx),tidx(nanidx),'linear','extrap');
            d_kd(nanidx) = nanval;
            d_kd = resample(d_kd,fs_spec,fs_kd);
            d_kd = d_kd(1:size(ss,1));
          
            fs_am = 15;
            d_am = log10(smoothdata(obj.DTable.d_amb{walknum},'movmean',fs_am*sm_x)+eps);
            ntp_am = obj.DTable.ntp_amb{walknum}; 
            [~,tidx1] = min(abs(ntp_am-ntp_wk_beg));
            tidx = tidx1:length(d_am);
            d_am = d_am(tidx);
            nanidx = isnan(d_am);
            nanval = interp1(tidx(~nanidx),d_am(~nanidx),tidx(nanidx),'linear','extrap');
            d_am(nanidx) = nanval;
            d_am = resample(d_am,fs_spec,fs_am);
            d_am = d_am(1:size(ss,1));
           
            C = corrcoef([ss,d_xs,d_gz,d_kd,d_am]);
            Clabels = {'theta','alpha','beta','gamma','hg','vel','gaze','kde','amb'};
          
            fH = figure;
            aH = axes('Parent',fH);
            imagesc(C,'Parent',aH);
            clim(aH,[-1,1])
            set(aH,'xtick',1:length(Clabels),'xticklabel',Clabels,'ytick',1:length(Clabels),'yticklabel',Clabels)
            title(sprintf('patient: %0.0f, walk: %0.0f, chan: %0.0f',patient,walknum,chan));
            cb = colorbar(aH);
            ylabel(cb,'correlation')
           
            if p.saveflag
                fstr = sprintf('%s_Walk%0.0f_Chan%0.0f_WalkConfMatrix.png',obj.PatientID,walknum,chan);
                print(fH,fullfile('D:\RWN\ConfMatrix',fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function plotMultInOutFreq(obj,varargin)
            %Generates a boxplot of power comparing in/out segments across
            %all patients and channels. Also runs glme for comparison.
            %
            %plotMultInOutDiff('fidx',1,'saveflag',true);
            p = obj.parseInputs(varargin{:});
            if isempty(obj.MultInOut)
                obj.getMultSegData;
            end
            if isempty(p.fidx)
                error('fidx is missing!')
            end

            SSTrials = obj.MultInOut.SSTrials; %TD: added Fix and Amb 20240418 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)

            % flist = [[2,5];[5,8];[8,12];[12,30];[30,60];[60,85];[85,120]];
            flist = [5,12];
            frng = flist(p.fidx,:);
            f = obj.MultInOut.f;
            fidx = f>frng(1)&f<=frng(2);

            %Mean across freq (5-8Hz) in log10 space
            fstr = ['Pwr',regexprep(num2str(frng),'\s+','')];
            SSTrials = addvars(SSTrials,mean(SSTrials.nPwr(:,fidx),2),'NewVariableNames',fstr,'Before','Pwr');

            %Creating table that is mean across uPtChan for glme
            SS = []; %side-by-side comparison of in/out
            SSmdl = []; %in/out is inline for glme
            for k=1:size(obj.MultInOut.uPtChan,1)
                idx = (obj.MultInOut.uPtChanIdx==k);

                pt = obj.MultInOut.uPtChan.Pt(k);
                chan = obj.MultInOut.uPtChan.Chan(k);
                pwr = SSTrials.nPwr(idx,:); %log10 already applied
                pwr_rng = SSTrials.(fstr)(idx); %log10 already applied
                inflag = SSTrials.InFlag(idx);

                in_n = sum(inflag);
                out_n = sum(~inflag);

                in_pwr = median(pwr(inflag,:));
                out_pwr = median(pwr(~inflag,:));

                in_pwr_rng = median(pwr_rng(inflag));
                out_pwr_rng = median(pwr_rng(~inflag));

                SS = vertcat(SS,table(...
                        k,pt,chan,...
                        in_n,out_n,...
                        in_pwr_rng,out_pwr_rng,...
                        in_pwr,out_pwr,...
                        'VariableNames',{...
                        'uPtChan','Pt','Chan',...
                        'In_N','Out_N',...
                        ['In_',fstr],['Out_',fstr],...
                        'In_Pwr','Out_Pwr',...
                        }));

                SSmdl = vertcat(SSmdl,table(...
                        [k;k],[pt;pt],[chan;chan],...
                        [true;false],[in_n;out_n],...
                        [in_pwr_rng;out_pwr_rng],[in_pwr;out_pwr],...
                        'VariableNames',{...
                        'uPtChan','Pt','Chan',...
                        'InFlag','InFlag_N',...
                        fstr,'Pwr',...
                        }));
            end

            % glme_inflag_ptchan = fitglme(SSmdl,[fstr,' ~ InFlag + (1|uPtChan)']);
            % Adjusted = calcAdjRespGLMEFit(SSmdl,glme_inflag_ptchan,'InFlag');
            % 
            % y = Adjusted.Response;
            % x = logical(SSmdl.InFlag);
            % M = nan(length(x),2);
            % y_in = y(x);
            % y_out = y(~x);
            % M(1:length(y_in),1) = y_in;
            % M(1:length(y_out),2) = y_out;
            % 
            % x = [1+(rand(size(M,1),1)*2-1)/8,2+(rand(size(M,1),1)*2-1)/8];
            % 
            % figure;
            % boxplot(M,{'Inside','Outside'});
            % hold on;
            % plot(x(:,1),M(:,1),'ob');
            % plot(x(:,2),M(:,2),'ob');
            % plot(x',M',':b');
            % 
            % [p,tbl,stats] = kruskalwallis(M,{'Inside','Outside'},'on'); %use raw responses here since these are stats by subject
            % figure;
            % [c,m] = multcompare(stats,'alpha',0.05,'display','on');
% 
            % %Significance for average unique pt/chan response (similar to
            % %boxplot ttest)
            % P = []; T = []; MD = []; SE = []; 
            % for n=1:length(f)
            %     ss = 10.^([SS.In_Pwr(:,n),SS.Out_Pwr(:,n)]);
            %     ss = ss./mean(ss,2);
            %     [h,p,ci,stats] = ttest2(ss(:,1),ss(:,2));
            %     P(n) = p;
            %     T(n) = stats.tstat;
            %     MD(n,:) = mean(ss);
            %     SE(n,:) = std(ss)./sqrt(size(ss,1));
            % end
            % T(P.*length(f)>0.05) = 0;
            % 
            % se = strel('rectangle',[3 1]);
            % PMsk = imdilate(reshape(T~=0,[],1),se);
            % PMsk = imerode(PMsk,se);
            % B = regionprops(PMsk(:)','BoundingBox');
            % I = regionprops(PMsk(:)','PixelIdxList');
            % B = cat(1,B.BoundingBox);
            % disp(B(:,3))
            % if ~isempty(B)
            %     idx = B(:,3)<3;
            %     B(idx,:) = [];
            %     rm_idxs = I(idx,:);
            %     rm_idxs = cat(1,rm_idxs.PixelIdxList);
            %     T(rm_idxs) = 0;
            % end
            % 
            % figure
            % plot(f,T)
            % xlabel('Hz')
            % ylabel('tstat')
            % title('Inside/Outside Difference Across Freq (all pts, all chans, t-test, bonferonni corr)')
            % 
            % fH = figure;
            % aH = axes('Parent',fH);
            % mInPwr = smoothdata(MD(:,1)'); sInPwr = smoothdata(SE(:,1))'.*2;
            % mOutPwr = smoothdata(MD(:,2)'); sOutPwr = smoothdata(SE(:,2))'.*2;
            % hold on;
            % ylimit = [0.96,1.04];
            % pH = nan(size(B,1),1);
            % for m=1:size(B,1)
            %     b = B(m,[1,3]);
            %     x = floor([b(1),sum(b),sum(b),b(1)]);
            %     x(x<1) = 1; x(x>length(f)) = length(f);
            %     y = [ylimit(1),ylimit(1),ylimit(2),ylimit(2)];
            %     pH(m) = patch(f(x),y,'k','FaceAlpha',0.05,'EdgeColor','none','parent',aH);
            % end
            % plot(f,mInPwr,'b');
            % plot(f,mOutPwr,'r');
            % pH1 = patch([f,fliplr(f)],([mInPwr+sInPwr,fliplr(mInPwr-sInPwr)]),'b','facealpha',0.3,'edgecolor','none');
            % pH2 = patch([f,fliplr(f)],([mOutPwr+sOutPwr,fliplr(mOutPwr-sOutPwr)]),'r','facealpha',0.3,'edgecolor','none');
            % legend([pH1,pH2],{'Inside','Outside'})
            % xlabel('Hz')
            % ylabel('Normalized Power')
            % title('In/Out Spectrum Comparison (All Pts, All Chans, All Walks, 95% CI)')
            % xlim(round(f([1,end])))
            % ylim(ylimit)
            % 
            % fH = figure;
            % aH = axes('Parent',fH);
            % mOutInDiff = mOutPwr - mInPwr;
            % sOutInDiff = sqrt(sInPwr.^2 + sOutPwr.^2);
            % hold on;
            % ylimit = [-0.02,0.07];
            % pH = nan(size(B,1),1);
            % for m=1:size(B,1)
            %     b = B(m,[1,3]);
            %     x = floor([b(1),sum(b),sum(b),b(1)]);
            %     x(x<1) = 1; x(x>length(f)) = length(f);
            %     y = [ylimit(1),ylimit(1),ylimit(2),ylimit(2)];
            %     pH(m) = patch(f(x),y,'k','FaceAlpha',0.05,'EdgeColor','none','parent',aH);
            % end
            % plot(f,mOutInDiff,'b');
            % plot(round(f([1,end])),[0,0],'k')
            % patch([f,fliplr(f)],([mOutInDiff+sOutInDiff,fliplr(mOutInDiff-sOutInDiff)]),'b','facealpha',0.3,'edgecolor','none');
            % xlabel('Hz')
            % ylabel('Normalized Power Difference (Outside-Inside)')
            % title('Spectrum Difference (All Pts, All Chans, All Walks, 95% CI)')
            % xlim(round(f([1,end])))
            % ylim(ylimit)

            %Standard boxplot with connecting lines
            ss = [SS.(['In_',fstr]),SS.(['Out_',fstr])]; %log10 already applied
            ss = ss./mean(ss,2);
            [h,pval,ci,stats] = ttest2(ss(:,1),ss(:,2));

            InOutTable = removevars(SS,{'In_Pwr','Out_Pwr',['In_',fstr],['Out_',fstr]});
            InOutTable = addvars(InOutTable,ss(:,1),ss(:,2),'NewVariableNames',{['In_',fstr],['Out_',fstr]});

            labels = {'Inside','Outside'};

            fH = figure('position',[50,50,850,650]);
            aH = axes('Parent',fH);
            boxplot(aH,ss,labels);
            ylimit = [min(reshape(cell2mat(get(findobj('Tag','Lower Whisker'),'ydata')),[],1)),max(reshape(cell2mat(get(findobj('Tag','Upper Whisker'),'ydata')),[],1))];
            ydiff = diff(ylimit)*0.1;
            ylimit = [ylimit(1)-ydiff,ylimit(2)+ydiff];
            hold(aH,'on');
            x = [1+(rand(size(ss,1),1)*2-1)/8,2+(rand(size(ss,1),1)*2-1)/8];
            pt_color = {'r','g','b','c','m'};
            rg_style = {'o','x','*','+'}; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para'
            rg_size = [7,10,10,10];
            for m=1:size(ss,1)
                rg_idx = obj.RegionTable.Region(obj.RegionTable.Patient==SS.Pt(m) & obj.RegionTable.Chan==SS.Chan(m));
                plot(aH,x(m,:),ss(m,:),[rg_style{rg_idx},pt_color{SS.Pt(m)}],'MarkerSize',rg_size(rg_idx),'LineWidth',1.5,'LineStyle','none');
                plot(aH,x(m,:),ss(m,:),pt_color{SS.Pt(m)},'MarkerSize',rg_size(rg_idx),'LineWidth',0.5,'LineStyle',':');
            end
            pH_pt = [];
            for m=1:length(pt_color)
                pH_pt(m) = plot(aH,[nan,nan],[nan,nan],pt_color{m},'LineWidth',1.5);
            end
            pH_rg = [];
            for m=1:length(rg_style)
                pH_rg(m) = plot(aH,nan,nan,[rg_style{m},'k'],'LineWidth',1.5,'MarkerSize',rg_size(m));
            end
            fstr2 = [regexprep(num2str(frng),'\s+','-'),' Hz'];
            ylabel(aH,['Normalized Power (',fstr2,')'])
            set(aH,'FontSize',14);
            legend([pH_pt,pH_rg],{'Patient1','Patient2','Patient3','Patient4','Patient5','AntHipp','LatTemp','Ent+Peri','PostHipp'},'location','northwest')

%             [p,tbl,stats] = kruskalwallis(ss,labels,'on'); %use raw responses here since these are stats by subject
%             figure;
%             [c,m] = multcompare(stats,'alpha',0.05,'display','on');
            
            title(aH,sprintf('Inside vs Outside Power (%s, All Chans, p=%0.1e)',fstr2,pval))
            ylim(aH,ylimit)

            if p.saveflag
                fstr = sprintf('InOut_Comparison_%d-%dHz_wAllLabels.png',frng(1),frng(2));
                fdir = fullfile('D:\RWN\InOutFreq');
                if ~isfolder(fdir)
                    mkdir(fdir);
                end
                print(fH,fullfile(fdir,fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function plotMultInOutPred(obj,varargin)
            %Generates a boxplot of power comparing in/out segments across
            %all patients and channels for specified predictor.
            %
            %plotMultInOutPred('predtype','vel','saveflag',true); %Vel, Fix, KDE, nAmb (this is the only one that is log10 normalized)
            p = obj.parseInputs(varargin{:});
            if isempty(obj.MultInOut)
                obj.getMultSegData;
            end
            if isempty(p.predtype)
                error('predtype is missing!')
            end

            SSTrials = obj.MultInOut.SSTrials; %TD: added Fix and Amb 20240418 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)

            %Creating table that is mean across uPtChan for glme
            SS = []; %side-by-side comparison of in/out
            SSmdl = []; %in/out is inline for glme
            for k=1:size(obj.MultInOut.uPtWalk,1)
                idx = (obj.MultInOut.uPtWalkIdx==k);

                sstrials = SSTrials(idx,:);
                sstrials = sstrials(sstrials.Chan==1,:);

                pt = obj.MultInOut.uPtWalk.Pt(k);
                walk = obj.MultInOut.uPtWalk.Walk(k);

                pred = sstrials.(p.predtype); 
                inflag = sstrials.InFlag;

                in_n = sum(inflag);
                out_n = sum(~inflag);

                in_pred = median(pred(inflag),'omitnan');
                out_pred = median(pred(~inflag),'omitnan');

                SS = vertcat(SS,table(...
                        k,pt,walk,...
                        in_n,out_n,...
                        in_pred,out_pred,...
                        'VariableNames',{...
                        'uPtWalk','Pt','Walk',...
                        'In_N','Out_N',...
                        'In_Pred','Out_Pred',...
                        }));

                SSmdl = vertcat(SSmdl,table(...
                        [k;k],[pt;pt],[walk;walk],...
                        [true;false],[in_n;out_n],...
                        [in_pred;out_pred],...
                        'VariableNames',{...
                        'uPtWalk','Pt','Walk',...
                        'InFlag','InFlag_N',...
                        'Pred',...
                        }));
            end

            %Standard boxplot with connecting lines
            ss = [SS.In_Pred,SS.Out_Pred]; 
            ss = ss./mean(ss,2);
            [h,pval,ci,stats] = ttest2(ss(:,1),ss(:,2));

            InOutTable = removevars(SS,{'In_Pred','Out_Pred'});
            InOutTable = addvars(InOutTable,ss(:,1),ss(:,2),'NewVariableNames',{'In_Pred','Out_Pred'});

            labels = {'Inside','Outside'};

            fH = figure('position',[50,50,850,650]);
            aH = axes('Parent',fH);
            boxplot(aH,ss,labels);
            ylimit = [min(reshape(cell2mat(get(findobj('Tag','Lower Whisker'),'ydata')),[],1)),max(reshape(cell2mat(get(findobj('Tag','Upper Whisker'),'ydata')),[],1))];
            ydiff = diff(ylimit)*0.1;
            ylimit = [ylimit(1)-ydiff,ylimit(2)+ydiff];
            hold(aH,'on');
            x = [1+(rand(size(ss,1),1)*2-1)/8,2+(rand(size(ss,1),1)*2-1)/8];
            pt_color = {'r','g','b','c','m'};
            for m=1:size(ss,1)
                plot(aH,x(m,:),ss(m,:),['.',pt_color{SS.Pt(m)}],'MarkerSize',20,'LineWidth',1.5,'LineStyle','none');
                plot(aH,x(m,:),ss(m,:),pt_color{SS.Pt(m)},'MarkerSize',20,'LineWidth',0.5,'LineStyle',':');
            end
            pH_pt = [];
            for m=1:length(pt_color)
                pH_pt(m) = plot(aH,[nan,nan],[nan,nan],['.',pt_color{m}],'LineWidth',1.5,'MarkerSize',20,'LineStyle','none');
            end
            ylabel(aH,['Normalized ',p.predtype])
            legend(pH_pt,{'Patient1','Patient2','Patient3','Patient4','Patient5'},'location','northwest')

%             [p,tbl,stats] = kruskalwallis(ss,labels,'on'); %use raw responses here since these are stats by subject
%             figure;
%             [c,m] = multcompare(stats,'alpha',0.05,'display','on');
            
            title(aH,sprintf('Inside vs Outside %s (All Walks, p=%0.1e)',p.predtype,pval))
            ylim(aH,ylimit)

            if p.saveflag
                fstr = sprintf('InOut_Comparison_%s_wAllLabels.png',p.predtype);
                fdir = fullfile('D:\RWN\InOutPred');
                if ~isfolder(fdir)
                    mkdir(fdir);
                end
                print(fH,fullfile(fdir,fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function plotMultInOutGLME(obj,varargin)
            %Generates a boxplot of power comparing in/out segments across
            %all patients and channels. Also runs glme for comparison.
            %
            %plotMultInOutDiff;
            %plotMultInOutDiff('saveflag',true);
            if isempty(obj.MultInOut)
                obj.getMultSegData;
            end

            p = obj.parseInputs(varargin{:});

            SSTrials = obj.MultInOut.SSTrials; %TD: added Fix and Amb 20240418 and KDE 20240513 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)
            f = obj.MultInOut.f;

            %Model list
            pred_list = {'InFlag','Vel','nAmb','Fix','KDE'};
            mdl_template = {};
            for k=1:length(pred_list)
                idx = nchoosek(1:length(pred_list),k);
                for n=1:size(idx,1)
                    mdl = ['nPwr ~ ',cell2mat(cellfun(@(x)[x,' + '],pred_list(idx(n,:)),'UniformOutput',false)),'(1|uPtChan)'];
                    mdl_template = cat(1,mdl_template,mdl);
                end
            end
            load('GLME_sidx.mat','sidx');

            flist = [[2,5];[5,8];[8,12];[12,30];[30,60];[60,85];[85,120]];
            LL = nan(length(mdl_template),size(flist,1));
            AC = nan(length(mdl_template),size(flist,1));
            TS = cell(length(mdl_template),size(flist,1));
            for m=[2,1,3,4,5,6,7]%1:size(flist,1)
                frng = flist(m,:);
                fidx = f>frng(1)&f<=frng(2);

                %Mean across freq (5-8Hz) in log10 space
                resp_var = ['nPwr',regexprep(num2str(frng),'\s+','')];
                SSTrials = addvars(SSTrials,mean(SSTrials.nPwr(:,fidx),2),'NewVariableNames',resp_var,'Before','Pwr'); %nPwr is log10

                %Model list
                mdl_list = regexprep(mdl_template,'nPwr',resp_var);

                for k=1:length(mdl_list)
                    disp(mdl_list{k});
                    glme = fitglme(SSTrials,mdl_list{k});
                    if p.saveflag
                        fstr = sprintf('GLME_ModelPrintout_%s_%d-%dHz.txt',regexprep(regexprep(mdl_list{k},'\s+',''),'\+\(1\|uPtChan\)',''),frng(1),frng(2));
                        fdir = fullfile('D:\RWN\GLME\ModelPrintouts');
                        if ~isfolder(fdir)
                            mkdir(fdir);
                        end
                        fid = fopen(fullfile(fdir,fstr),'w');
                        out_str = regexprep(evalc('disp(glme)'),'\%|<strong>|</strong>','');
                        out_str = sprintf('%s\nordR2 = %0.4f\nadjR2 = %0.4f\n',out_str,glme.Rsquared.Ordinary,glme.Rsquared.Adjusted);
                        fprintf(fid,out_str);
                        fclose(fid);
                    end
                    LL(k,m) = glme.LogLikelihood;
                    AC(k,m) = glme.ModelCriterion.AIC;
                    TS(k,m) = {[glme.Coefficients.Name,num2cell(glme.Coefficients.tStat),num2cell(glme.Coefficients.pValue)]};
                end

                % if m==2
                %     [~,sidx] = sort(AC(:,m),'ascend');
                % end
                LL_sorted = LL(sidx,m);
                AC_sorted = AC(sidx,m);
                TS_sorted = TS(sidx,m);
                mdl_list_sorted = mdl_list(sidx);

                for k=1:length(mdl_list_sorted)
                    ts = TS_sorted{k}(2:end,:);
                    mdl = mdl_list_sorted{k};
                    for n=1:size(ts,1)
                        pval = ts{n,3};
                        pred = regexprep(ts{n,1},'_1','');
                        if pval<0.05
                            mdl = regexprep(mdl,pred,['\\bf ',pred,' \\rm ']);
                        end
                    end
                    mdl_list_sorted{k} = mdl;
                end

                if p.saveflag
                    fH = figure('position',[50,50,1800,1000],'Visible','off');
                else
                    fH = figure('position',[50,50,1800,1000],'Visible','on');
                end
                aH = axes('Parent',fH);
                plot(aH,AC_sorted);
                xpos = length(mdl_list_sorted)/2.5;
                ylimit = get(aH,'ylim');
                ypos = ylimit(2)-diff(ylimit)/5;
                text(xpos,ypos,[{'PRED'};TS_sorted{1}(:,1)],'parent',aH,'Interpreter','none','HorizontalAlignment','center');
                text(xpos+3,ypos,[{'TSTAT'};TS_sorted{1}(:,2)],'parent',aH,'HorizontalAlignment','center');
                text(xpos+6,ypos,[{'PVAL'};TS_sorted{1}(:,3)],'parent',aH,'HorizontalAlignment','center');
                set(aH,'xtick',1:length(mdl_list_sorted),'xticklabel',mdl_list_sorted)
                ylabel(aH,'AIC')
                xlim(aH,[0,length(mdl_list_sorted)+1])

                if p.saveflag
                    fstr = sprintf('GLME_AllCombos_FixedSort_%d-%dHz.png',frng(1),frng(2));
                    fdir = fullfile('D:\RWN\GLME');
                    if ~isfolder(fdir)
                        mkdir(fdir);
                    end
                    print(fH,fullfile(fdir,fstr),'-dpng','-r300');
                    close(fH);
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%
    methods %figures, etc.

        function figureRW1SpecGram(obj,varargin)
            %plots specgrams for all walks for RW1 for event where
            %pt chooses left/right. Chan 4 shows most consistent resp.
            % Walk1 = 18531 (frame) 
            % Walk2 = 27870 (frame) 
            % Walk3 = 20253 (frame) 
            % Walk4 = 19291 (frame) 
            % Walk5 = 17403 (frame) 
            % Walk6 = 18463 (frame) 

            pt = 1;
            chan = 4;

            walk_list = (1:6);
%             frame_list = [18531,27870,20253,19291,17403,18463]; %aligned to choice?
            frame_list = [18191,27580,19973,18931,17123,18193]; %aligned to position for left/right choice
%             frame_list = [4593,13299,6022,5105,4087,4022]; %aligned with overhead shower station

            winsec = 32;

            fs_np = obj.DTable.fs_np(1);
            nsamp = winsec*fs_np;
            tsamp = (-nsamp:nsamp); %+-32 sec (2sec removed for artifact, so +-30sec)
            tsec = tsamp/fs_np;

            fH = figure('Position',[50,50,1700,500],'Visible','on');
            aH = []; 
            for k=1:length(walk_list)
                aH(k) = subplot(length(walk_list),1,k,'parent',fH);

                wlk = walk_list(k);
                frm = frame_list(k);
                [~,idx] = min(abs(obj.DTable.ntp_pupil{wlk}(frm)-obj.DTable.ntp_np{wlk}));
                tfull = idx+tsamp; 
                d = obj.DTable.d_np{wlk}(tfull,4);
                d(d<-500) = 0;
                d(isnan(d)) = 0;
                lb = obj.ChanLabels{pt,chan};
                PSG = PermSpecGram(d,'Fs',250,'FPass',[2,64],'TimeRng',[tsec(1),tsec(end)],'NormRng',[tsec(1)+2,tsec(end)-2]);
%                 PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','Smoothing',[1000,5],'ErrPerc',[]);
                PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','Smoothing',[500,1],'ErrPerc',[]);
                PSG.plotSpecGram('ShowRaw',true,'Clim',[-5,5],'aH',aH(k));
                plot(aH(k),[0,0],[0,8],'k'); %center line
%                 xlim(aH(k),[tsec(1)+2,tsec(end)-2])
                xlim(aH(k),[-20,20])
                title(aH(k),lb)
            end


        end

    end

    %%%%%%%%%%%%%%%%%
    methods %callbacks, helper fcns

        function keyEvntFcn(~,varargin)
            ky = varargin{2}.Key;
            md = cell2mat(varargin{2}.Modifier);

            if strcmp(ky,md)
                return;
            end

            if ~isempty(md)
                ky = [md,'+',ky]; %adds modifier to key value (i.e. control+s)
            end

            SS = getappdata(varargin{1},'SS');
            switch ky
                case 'rightarrow'
                    SS.cf = SS.cf + 1;
                case 'shift+rightarrow'
                    SS.cf = SS.cf + 10;
                case 'control+rightarrow'
                    SS.cf = SS.cf + 100;
                case 'leftarrow'
                    SS.cf = SS.cf - 1;
                case 'shift+leftarrow'
                    SS.cf = SS.cf - 10;
                case 'control+leftarrow'
                    SS.cf = SS.cf - 100;
            end
            SS.cf(SS.cf<1) = 1; SS.cf(SS.cf>length(SS.tframes)) = length(SS.tframes);
            SS.iH2.CData = read(SS.vid,SS.tframes(SS.cf));
            title(SS.aH2,sprintf('Frame = %0.0f',SS.tframes(SS.cf)))
            for k=1:4
                set(SS.pH(k),'xdata',[SS.tfrsec(SS.cf),SS.tfrsec(SS.cf)])
            end
            setappdata(SS.fH,'SS',SS)
            setappdata(SS.fH2,'SS',SS)
        end

    end

end %class





