classdef RWAnalysis2 < handle
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
    % X Redo in/out boxplot to incorporate lensky and change marker shape by brain region (Done)
    % X Baseline indoor/indoor transitions (Done)
    % X Combine hippocampus into one (hipp vs all others) (Done - added | symbol for regions)
    % X Need to mark stopping at a closed door vs walking through open door for all doors (Done - keywords added)
    % X Incorporate luminance, fixation rate (done)
    % 
    % X Add ability to run multiple "Transitions" in the same spectrogram (Correct Turn Beg|Incorrect Turn Beg) (Done)
    % X Plot band filtered individual trials color coded by participant (Done)
    % X Incorporate the Pupil Eye tracking data in to GLMM model (Done)
    % X Incorporate the Real World Nav Event Segmentation into the GLMM Model (Done)
    % X Make tool to create event times for RWN Event Segmentation to do event based spectrograms at different levels of participant agreement (>15 people hit the spacebar at these time points, <5 people hit the spacebar at these time points, etc.)
    % X Remove outliers before filtering (Done)
    % X Adjustable smoothing window for PermSpecGram option in GUI (Done)
    % X Save figure checkbox (Done)
    % X Separate figure for trials data (Done)
    % X Print threshold for clusters**** (Does not make sense)
    % X Enable TransitionDD/RegionDD for difference in GUI (Done)
    % X Trials plot for other permutation functions (Done)
    % X Make new permutation code to use baseline and FDR (Done)
    % X Make variable window size (Done)
    % X Use sum(cluster values) in the correction instead of just size (Done)
    % X Increase freq range for new permutation code (Done)
    % X Variable window size for GLM (no need to vary, using 2sec default)
    % X Event related GLM (?)
    % X Leave one for Patient drop down in GUI (Done)
    % X Boxplot comparison (Done)
    %
    % X Add KDE to GLM (Done)
    % X Spectrograms around KDE peaks (high, mid, low, trough) (Done)
    % X Event types related to KDE data broken down into pie charts (student project?)
    % X Correlation confusion matrix of all bands and all predictors (by channel i.e. 20 for all walks) (done)
    % X GLM outputs to text (done)
    % X Walk spectrograms with all predictors (done)
    % X Boxplot for each of the predictors/freq(pwr) by indoor/outdoor (done)
    % X The updated Kiersten data has new events that change the glme results slightly (not sure how...)
    % 
    % X Every KDE peak as an option (done)
    % X KDE peaks that don't have a speed change (velocity dropdown accomplishes this) (done)
    % X Ratio between before/after in specgram power for gui (add a average value to title) (done)
    % X Only pick one electrode per participant in AntHipp to prevent double dipping (add custom dropdown and new window with list of all chans from all patients) (done)
    % X Checkbox for including velocity (done)
    % X Iterate and save for each category (done)
    % X Add tercile for velocity (add to title) (done)
    % X Add check for overlapping transitions and include value in figure (percentage of trials and average time of overlap) (done)
    % X Add more detail to figure titles for saving (done)
    %
    % X 1/f tilt calculation across time (or before/after) for each transition specgram (fooof v2.0 donaghue) in gui (done)
    % X Correlation confusion matrix - run stats and combine by region, patient, etc (new gui?) (done)
    % X Plot average KDE value for each transition (break down doorways - indoor/indoor, indoor/outdoor, outdoor/indoor, outdoor/outdoor,closed, open, etc. -> RWN_annotationbrackets.xls) (done)
    % X Comparing KDE to expert rater (how much lag?) (done)
    % 
    % X Add another velocity condition to completely separate change from no change (done)
    % X Checkbox to turn off specgram and move predictor axis to left (done)
    % X Save raw data as mat or xls (done)
    % X Add skin conductance (cleaned raw, phasic - plot both) (x,y,z,ppg sp02 signal,phasic eda,resp,ecg,raw cleaned eda, 1000Hz) (done)
    % X Fix plotcomp (done)
    % X Add predictor in filename for saving in gui (done)
    % X Iterate through plot predictors (done)
    %
    % X Exclude overlap trials in gui checkbox (done)
    % X GLME filtered by region (new gui?) (done)
    % X Align to peaks in eda (done)
    % X Change ylim for predictors in GUI -> vel:0-1.5, gaze:0-1, kde:1-7,amb: 0-500, eda:-0.2-0.5 (done)
    % X Velocity change low/high 10% (done)
    % X Save 1/f plot during automatic saving (done)
    %
    % X add individual R2 values for each predictor in glme result (leave one predictor out, backwards step regression) (done)
    % X add multitaper back into table for comparison (done)
    %
    % X fooof and pepisode added to 2sec full walk table
    % X eBosc pepisode (2-20Hz, 2-5, 5-8, 8-12) and alignment to events
    % X calc velocity change by patient instead of entire group (done)
    % X Recalculate indoor/outdoor spectrum plots from long ago (done)
    % *Pepisode for 12-16Hz, fooof offset in same plot as exponent
    % *Try out new dynamic fooof (specparam)
    % *mtAlpha ~ Vel + Fix + Amb + KDE + EDA + (1|uPtChan) + (1|Pt) + (Vel + Fix + Amb + KDE + EDA|Walk) (also fix code)
    % *Add theta(4-8)/gamma(30-85) ratio to boxplot to compare with fooof
    % *Add heart rate, respiration rate, *eye tracking, *head movements, etc. 
    % *add in/out and conf matrix plots to gui and boxplot comparison with pepisode
    % Correct/Incorrect turns gamma timing with low freq suppression (phase amp coupling) (GUI?) (captured by fooof?)
    % Add "Body Turn Beg" and "Body Turn End" to RWNGUI.
    % 5/2, 5/6, 5/8 are missing drift table for chest phone

    %%%%%%%%%%%%%%%%%%%
    properties
        RootDir; PatientID; DB; WalkNames; WalkNums; DParsed; PatientList; PatientIdx;
        ChanLabels; InOutTimes; DTable; DisabledWalks; FB; Vid; RegionTable;
        StopGoWalks; AnalysisFile; WinSegSec; WinTransSec; 
        MultTable; MultSeg; MultTrans; 
    end

    %%%%%%%%%%%%%%%%%%%
    methods %initialization, calculations, etc

        function obj = RWAnalysis2(varargin) %constructor (runs when object is created for the first time)
            obj.parseInputs(varargin{:});

            % obj.RootDir = '\\155.100.91.44\D\Tyler\RealWorld';
            obj.RootDir = 'D:\RWN';
            obj.WinSegSec = 2; %Also smoothing kernel for specgrams in trans analysis
            obj.WinTransSec = 12;
            obj.AnalysisFile = fullfile(obj.RootDir,['RWAnalysis2_MedianNorm_',num2str(obj.WinSegSec,'%d'),'sec.mat']);

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
            p.fidx = []; %index into list of freq bins
            p.ftype = ''; %Theta, Alpha, Beta, Gamma, HG (plotMultInOutFreq)
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
            p.plotpowercomp = false; %boxplot comparison of mean power across trials for two rectangular regions in a spectrogram
            p.plotfooofcomp = false; %boxplot comparison of fooof exponenet across trials for x-values of the two rectangular regions (y-value or freq is not constrained)
            p.boxcomprng = []; %rows are box1/box2, cols are time/freq limits for each box
            p.predtype = ''; %predictor type -> Vel, Fix, KDE, Amb, EDA (plotMultInOutPred)
            p.veltype = ''; %VelHigh, VelLow, VelHighTercile, VelMidTercile, VelLowTercile
            p.rmoverlap = false; %remove overlapping trials
            p.rmoutliers = false; %remove outliers in in/out boxplot data
            p.glmemdl = ''; %glme model (i.e. nTheta ~ InFlag + (1|uPtChan)

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
            obj.FB = [];

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
            
            obj.DTable = cell(TotalWalks,31);
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
                        if isempty(SS.fs_xs); obj.DTable(k,6) = {100}; else; obj.DTable(k,6) = {SS.fs_xs}; end

                        obj.DTable(k,7) = {SS.ntp_pupil};
                        obj.DTable(k,8) = {SS.ntp_gp};

                        if isfield(SS,'ntp_gaze')
                            obj.DTable(k,9) = {SS.d_gaze_fix};
                            % obj.DTable(k,9) = {smoothdata(SS.d_gaze_fix,1,'movmean',200*obj.WinSegSec)}; %smooth to match?
                            obj.DTable(k,10) = {SS.ntp_gaze};
                            if isempty(SS.fs_gaze); obj.DTable(k,11) = {200}; else; obj.DTable(k,11) = {SS.fs_gaze}; end
                        end

                        if isfield(SS,'ntp_amb')
                            obj.DTable(k,12) = {SS.d_amb};
                            obj.DTable(k,13) = {SS.ntp_amb};
                            if isempty(SS.fs_amb); obj.DTable(k,14) = {15}; else; obj.DTable(k,14) = {SS.fs_amb}; end
                        end

                        if isfield(SS,'ntp_kde')
                            %KDEPeakHighTercile, KDEPeakMidTercile, KDEPeakLowTercile,
                            %KDEPeakHigh, KDEPeakLow, KDETrough, KDEPeakAll
                            %(trough is defined as trs<prctile(pks,5))
                            obj.DTable(k,15) = {SS.d_kde};
                            obj.DTable(k,16) = {SS.ntp_kde};
                            if isempty(SS.fs_kde); obj.DTable(k,17) = {60}; else; obj.DTable(k,17) = {SS.fs_kde}; end
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
                        fprintf('Calculating wavelet transform for %s walk %d\n',obj.PatientID,k)
                        d = SS.d_np; nan_idx = isnan(d); d(nan_idx) = 0;
                        [f,~,cfs] = morseSpecGram(d,SS.fs_np,[2,120]); %2 to 120Hz
                        cfs(permute(repmat(nan_idx,1,1,length(f)),[1,3,2])) = nan; %time x freq x chan
                        pwr = abs(cfs).^2;
                        pwr = smoothdata(pwr,1,'movmean',250*obj.WinSegSec); %smooth across time (this handles nans better than smooth3) (~0.22Hz -> 0.443/2sec cutoff for 2sec window -> might want to decrease this!)
                        pwr = pwr./median(pwr,1,"omitnan"); %full epoch norm (time x freq x chan)

                        %Realistically only need a 12.5Hz lowpass antialias
                        %filter for this downsample factor, not 0.22Hz.
                        %12.5Hz would be a boxcar filter of length
                        %0.443/12.5 = 0.0354sec or 0.0354*250 ~ 9samples.
                        obj.DTable(k,18) = {pwr(1:10:end,:,:)}; %smoothed/downsampled specgram power (factor of 10 so 25Hz) 
                        obj.DTable(k,19) = {SS.ntp_np(1:10:end)};
                        obj.DTable(k,20) = {SS.fs_np./10}; %fs=25Hz
                        obj.DTable(k,21) = {f};

                        if isfield(SS,'ntp_bio')
                            obj.DTable(k,22) = {SS.d_bio};
                            obj.DTable(k,23) = {SS.ntp_bio};
                            if isempty(SS.fs_bio); obj.DTable(k,24) = {100}; else; obj.DTable(k,24) = {SS.fs_bio}; end
                            if ~isempty(SS.ntp_bio)
                                %finding peaks
                                [pks,ploc] = findpeaks(SS.d_bio.PhasicEDA,'MinPeakProminence',0.05,'NPeaks',100); %peaks
                                [~,pidx] = sort(pks,'descend'); pks = pks(pidx); ploc = ploc(pidx); pks5 = prctile(pks,5);
                                [trs,tloc] = findpeaks(-SS.d_bio.PhasicEDA,'MinPeakProminence',0.05,'NPeaks',100); trs = -trs; %troughs
                                [~,tidx] = sort(trs,'descend'); trs = trs(tidx); tloc = tloc(tidx);
                                evnt_table = [];
                                for m=1:length(pks) %going highest to lowest
                                    ntp_bio = SS.ntp_bio(ploc(m));
                                    [~,pupilframe] = min(abs(SS.ntp_pupil - ntp_bio));
                                    [~,goproframe] = min(abs(SS.ntp_gp - ntp_bio));
                                    [~,npsample] = min(abs(SS.ntp_np - ntp_bio));
                                    evnt_table = vertcat(evnt_table,table({"EDAPeakAll"},{""},pupilframe,goproframe,npsample,ntp_bio,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    if m<floor(length(pks)/2)
                                        evnt_table = vertcat(evnt_table,table({"EDAPeakHigh"},{""},pupilframe,goproframe,npsample,ntp_bio,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    else
                                        evnt_table = vertcat(evnt_table,table({"EDAPeakLow"},{""},pupilframe,goproframe,npsample,ntp_bio,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    end
                                end
                                for m=1:length(trs)
                                    if trs(m)<pks5
                                        ntp_bio = SS.ntp_bio(tloc(m));
                                        [~,pupilframe] = min(abs(SS.ntp_pupil - ntp_bio));
                                        [~,goproframe] = min(abs(SS.ntp_gp - ntp_bio));
                                        [~,npsample] = min(abs(SS.ntp_np - ntp_bio));
                                        evnt_table = vertcat(evnt_table,table({"EDATrough"},{""},pupilframe,goproframe,npsample,ntp_bio,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    end
                                end
                                SS.evnts_tbl = vertcat(SS.evnts_tbl,evnt_table);
                                SS.evnts_tbl = sortrows(SS.evnts_tbl,"NTP");
                            end
                        end %if bio

                        %fBosc on full walk (stored as obj.FB)
                        if isempty(obj.MultTable)
                            fprintf('Calculating fBosc for %s walk %d\n',obj.PatientID,k)
                            obj.calcfBOSCWalk(SS.d_np,SS.fs_np);
                            obj.DTable(k,25) = {obj.FB.PEP(1:10:end,:,:)};
                            obj.DTable(k,26) = {SS.ntp_np(1:10:end)};
                            obj.DTable(k,27) = {SS.fs_np./10}; %fs=25Hz
                            obj.DTable(k,28) = {obj.FB.f_bands};
                        else
                            fprintf('Loading fBosc from MultTable for %s walk %d to save time\n',obj.PatientID,k)
                            dtable = obj.MultTable{obj.PatientIdx};
                            obj.DTable(k,25) = dtable.d_pep(WalkNum);
                            obj.DTable(k,26) = dtable.ntp_pep(WalkNum);
                            obj.DTable(k,27) = {dtable.fs_pep(WalkNum)};
                            obj.DTable(k,28) = dtable.f_pep(WalkNum);
                        end

                        %Pupil IMU
                        if isfield(SS,'ntp_imu')
                            obj.DTable(k,29) = {SS.d_imu};
                            obj.DTable(k,30) = {SS.ntp_imu};
                            if isempty(SS.fs_imu); obj.DTable(k,31) = {200}; else; obj.DTable(k,31) = {SS.fs_imu}; end
                            if ~isempty(SS.ntp_imu)
                                %finding peaks
                                [pks,ploc] = findpeaks(smoothdata(SS.d_imu.gyroY,'movmean',100),'MinPeakProminence',100,'NPeaks',200); %peaks for left/right head turn (left is negative, right is positive)
                                [trs,tloc] = findpeaks(-smoothdata(SS.d_imu.gyroY,'movmean',100),'MinPeakProminence',100,'NPeaks',200); trs = -trs; %troughs
                                rldir = [zeros(length(pks),1);ones(length(trs),1)]; pks = [pks;-trs]; ploc = [ploc;tloc]; %combine right(pos)/left(neg) head turns and make them all positive
                                [~,pidx] = sort(pks,'descend'); pks = pks(pidx); ploc = ploc(pidx); rldir = rldir(pidx); %0=right, 1=left
                                evnt_table = [];
                                for m=1:length(pks) %going highest to lowest
                                    ntp_imu = SS.ntp_imu(ploc(m));
                                    [~,pupilframe] = min(abs(SS.ntp_pupil - ntp_imu));
                                    [~,goproframe] = min(abs(SS.ntp_gp - ntp_imu));
                                    [~,npsample] = min(abs(SS.ntp_np - ntp_imu));
                                    evnt_table = vertcat(evnt_table,table({"HeadRotPeakAll"},{""},pupilframe,goproframe,npsample,ntp_imu,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    if m<floor(length(pks)/2)
                                        evnt_table = vertcat(evnt_table,table({"HeadRotPeakHigh"},{""},pupilframe,goproframe,npsample,ntp_imu,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    else
                                        evnt_table = vertcat(evnt_table,table({"HeadRotPeakLow"},{""},pupilframe,goproframe,npsample,ntp_imu,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    end
                                    if rldir(m) %left
                                        evnt_table = vertcat(evnt_table,table({"HeadRotPeakLeft"},{""},pupilframe,goproframe,npsample,ntp_imu,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    else %right
                                        evnt_table = vertcat(evnt_table,table({"HeadRotPeakRight"},{""},pupilframe,goproframe,npsample,ntp_imu,'VariableNames',{'Event','Description','PupilFrame','GoProFrame','NPSample','NTP'}));
                                    end
                                end
                                SS.evnts_tbl = vertcat(SS.evnts_tbl,evnt_table);
                                SS.evnts_tbl = sortrows(SS.evnts_tbl,"NTP");
                            end
                        end %if IMU

                        ptnum = str2double(regexp(obj.PatientID,'\d+$','match','once'));

                        SS.evnts_tbl =  addvars(SS.evnts_tbl,repmat(ptnum,size(SS.evnts_tbl,1),1),'NewVariableNames','Patient','Before','Event');
                        SS.evnts_tbl =  addvars(SS.evnts_tbl,repmat(WalkNum,size(SS.evnts_tbl,1),1),'NewVariableNames','Walk','Before','Event');

                        obj.DB = vertcat(obj.DB,SS.evnts_tbl);
                    end
                end
            end
            
            % This does not decrease comp time for some reason...
            % dtable = obj.DTable(:,1:3);
            % pcell = cell(TotalWalks,1);
            % fprintf('Calculating fBosc for all walks in parfor loop...\n');
            % parfor k=1:TotalWalks
            %     fB = calcfBOSCWalkSA(dtable{k,1},dtable{k,3});
            %     pcell{k} = fB.PEP(1:10:end,:,:);
            % end
            % toc;
            % for k=1:TotalWalks
            %     obj.DTable(k,25) = pcell(k);
            %     obj.DTable(k,26) = {dtable{k,2}(1:10:end)};
            %     obj.DTable(k,27) = {dtable{k,3}./10}; %fs=25Hz
            %     obj.DTable(k,28) = {[4,8;8,12;12,30]}; %Theta,Alpha,Beta;
            % end
        
            obj.DTable = cell2table(obj.DTable,'VariableNames',{...
                'd_np','ntp_np','fs_np',...
                'd_xs','ntp_xs','fs_xs',...
                'ntp_pupil','ntp_gp',...
                'd_gaze_fix','ntp_gaze','fs_gaze',...
                'd_amb','ntp_amb','fs_amb',...
                'd_kde','ntp_kde','fs_kde',...
                'd_wav','ntp_wav','fs_wav','f_wav',...
                'd_bio','ntp_bio','fs_bio',...
                'd_pep','ntp_pep','fs_pep','f_pep',...
                'd_imu','ntp_imu','fs_imu',...
                });

        end %loadData

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

        function parseTransData(obj,varargin)
            %Parse data into +-12 sec segments centered on transitions.
            %Parses all transistions in the table. Run loadData first.
            %parseTransData;
            DT_np = []; DT_np_idx = []; %index into data for center of transition
            DT_xs = []; DT_xs_idx = []; DT_xs_vchg = []; %used to filter trials based on percent change in velocity
            DT_gz = []; DT_gz_idx = [];
            DT_am = []; DT_am_idx = [];
            DT_kd = []; DT_kd_idx = [];
            DT_wv = []; DT_wv_idx = [];
            DT_bi = []; DT_bi_idx = []; %biopac
            DT_pe = []; DT_pe_idx = []; %bosc pepisode
            DT_im = []; DT_im_idx = []; %imu
            EvntTrans = {}; %event type
            DescTrans = {}; %event description
            WalkNumTrans = []; %walk number for transitions
            NumSampTrans = []; %number of samples in data for each walk
            
            for m=1:size(obj.DTable,1) %by walk
                walknum = obj.WalkNums(m);
                numsamp = size(obj.DTable.d_np{m},1); %number of samples in walk data (needed to remove windows that extend beyond data limits)

                dbb = obj.DB(obj.DB.Walk==walknum,:);

                %NP
                d_np = obj.DTable.d_np{m};
                ntp_np = obj.DTable.ntp_np{m}; %ntp in sec
                fs_np = obj.DTable.fs_np(m);

                %XS
                d_xs = obj.DTable.d_xs{m};
                ntp_xs = obj.DTable.ntp_xs{m}; %ntp in sec
                fs_xs = obj.DTable.fs_xs(m);
                
                %Gaze
                d_gz = obj.DTable.d_gaze_fix{m};
                ntp_gz = obj.DTable.ntp_gaze{m}; %ntp in sec
                fs_gz = obj.DTable.fs_gaze(m);

                %Ambient
                d_am = obj.DTable.d_amb{m};
                ntp_am = obj.DTable.ntp_amb{m}; %ntp in sec
                fs_am = obj.DTable.fs_amb(m);

                %KDE
                d_kd = obj.DTable.d_kde{m};
                ntp_kd = obj.DTable.ntp_kde{m}; %ntp in sec
                fs_kd = obj.DTable.fs_kde(m);

                %Wavelet
                d_wv = obj.DTable.d_wav{m};
                ntp_wv = obj.DTable.ntp_wav{m}; %ntp in sec
                fs_wv = obj.DTable.fs_wav(m);
                f_wv = obj.DTable.f_wav{m};

                %Biopac
                if isempty(obj.DTable.d_bio{m})
                    d_bi = [];
                else
                    d_bi = table2array(obj.DTable.d_bio{m});
                end
                ntp_bi = obj.DTable.ntp_bio{m}; %ntp in sec
                fs_bi = obj.DTable.fs_bio(m);
                f_bi = obj.DTable.d_bio{find(~cellfun(@isempty,obj.DTable.d_bio),1)}.Properties.VariableNames;

                %Bosc pepisode
                d_pe = obj.DTable.d_pep{m};
                ntp_pe = obj.DTable.ntp_pep{m}; %ntp in sec
                fs_pe = obj.DTable.fs_pep(m);
                f_pe = obj.DTable.f_pep{m};

                %IMU
                d_im = struct2array(obj.DTable.d_imu(m));
                ntp_im = obj.DTable.ntp_imu{m}; %ntp in sec
                fs_im = obj.DTable.fs_imu(m);
                f_im = fieldnames(obj.DTable.d_imu(m));

                %12sec window around transition
                win_trans_sec = obj.WinTransSec; %sec
                ntp = dbb.NTP; %all transitions types
                evnt = dbb.Event;
                desc = dbb.Description;
                win_trans_np = (-win_trans_sec*fs_np:win_trans_sec*fs_np); d_np_trans = nan(length(win_trans_np),length(ntp),4); d_np_trans_idx = nan(length(ntp),1);
                win_trans_xs = (-win_trans_sec*fs_xs:win_trans_sec*fs_xs); d_xs_trans = nan(length(win_trans_xs),length(ntp)); d_xs_trans_idx = nan(length(ntp),1); 
                win_trans_gz = (-win_trans_sec*fs_gz:win_trans_sec*fs_gz); d_gz_trans = nan(length(win_trans_gz),length(ntp)); d_gz_trans_idx = nan(length(ntp),1);
                win_trans_am = (-win_trans_sec*fs_am:win_trans_sec*fs_am); d_am_trans = nan(length(win_trans_am),length(ntp)); d_am_trans_idx = nan(length(ntp),1);
                win_trans_kd = (-win_trans_sec*fs_kd:win_trans_sec*fs_kd); d_kd_trans = nan(length(win_trans_kd),length(ntp)); d_kd_trans_idx = nan(length(ntp),1);
                win_trans_wv = (-win_trans_sec*fs_wv:win_trans_sec*fs_wv); d_wv_trans = nan(length(win_trans_wv),length(ntp),length(f_wv),4); d_wv_trans_idx = nan(length(ntp),1);
                win_trans_bi = (-win_trans_sec*fs_bi:win_trans_sec*fs_bi); d_bi_trans = nan(length(win_trans_bi),length(ntp),8); d_bi_trans_idx = nan(length(ntp),1);
                win_trans_pe = (-win_trans_sec*fs_pe:win_trans_sec*fs_pe); d_pe_trans = nan(length(win_trans_pe),length(ntp),size(f_pe,1),4); d_pe_trans_idx = nan(length(ntp),1);
                win_trans_im = (-win_trans_sec*fs_im:win_trans_sec*fs_im); d_im_trans = nan(length(win_trans_im),length(ntp),length(f_im)); d_im_trans_idx = nan(length(ntp),1);
                t_xs = win_trans_xs./fs_xs; %for velocity change calculation
                t_xs_idx = t_xs>-5 & t_xs<5; %only focus on change in vel for +-5sec
                d_xs_vchg = nan(length(ntp),1); %percent velocity change
                for k=1:length(ntp)
                    %NP
                    [~,midx] = min(abs(ntp_np-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_np(1))>=1 && (midx+win_trans_np(end))<=length(d_np)
                            d_np_trans_idx(k) = midx;
                            d_np_trans(:,k,:) = d_np(midx+win_trans_np,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of NP data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %XS
                    [~,midx] = min(abs(ntp_xs-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_xs(1))>=1 && (midx+win_trans_xs(end))<=length(d_xs)
                            d_xs_trans_idx(k) = midx;
                            dx = d_xs(midx+win_trans_xs);
                            d_xs_trans(:,k) = dx;
                            mdx = median(dx);
                            d_xs_vchg(k) = (min(dx(t_xs_idx))-mdx)./mdx; %percent change in velocity relative to median
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of XS data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %Gaze
                    [~,midx] = min(abs(ntp_gz-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_gz(1))>=1 && (midx+win_trans_gz(end))<=length(d_gz)
                            d_gz_trans_idx(k) = midx;
                            d_gz_trans(:,k) = d_gz(midx+win_trans_gz);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of gaze data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %Ambient
                    [~,midx] = min(abs(ntp_am-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_am(1))>=1 && (midx+win_trans_am(end))<=length(d_am)
                            d_am_trans_idx(k) = midx;
                            d_am_trans(:,k) = d_am(midx+win_trans_am);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of ambient data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %KDE
                    [~,midx] = min(abs(ntp_kd-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_kd(1))>=1 && (midx+win_trans_kd(end))<=length(d_kd)
                            d_kd_trans_idx(k) = midx;
                            d_kd_trans(:,k) = d_kd(midx+win_trans_kd);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of kde data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %Wavelet
                    [~,midx] = min(abs(ntp_wv-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_wv(1))>=1 && (midx+win_trans_wv(end))<=size(d_wv,1)
                            d_wv_trans_idx(k) = midx;
                            d_wv_trans(:,k,:,:) = d_wv(midx+win_trans_wv,:,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of wavelet data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %Biopac
                    [~,midx] = min(abs(ntp_bi-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_bi(1))>=1 && (midx+win_trans_bi(end))<=size(d_bi,1)
                            d_bi_trans_idx(k) = midx;
                            d_bi_trans(:,k,:) = d_bi(midx+win_trans_bi,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of biopac data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %Bosc pepisode
                    [~,midx] = min(abs(ntp_pe-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_pe(1))>=1 && (midx+win_trans_pe(end))<=size(d_pe,1)
                            d_pe_trans_idx(k) = midx;
                            d_pe_trans(:,k,:,:) = d_pe(midx+win_trans_pe,:,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of bosc pepisode data. Filling with NaNs.\n',k,walknum);
                        end
                    end

                    %IMU
                    [~,midx] = min(abs(ntp_im-ntp(k)));
                    if ~isempty(midx)
                        if (midx+win_trans_im(1))>=1 && (midx+win_trans_im(end))<=size(d_im,1)
                            d_im_trans_idx(k) = midx;
                            d_im_trans(:,k,:) = d_im(midx+win_trans_im,:);
                        else
                            fprintf('Trans segment %0.0f for walk %0.0f is out of range of imu data. Filling with NaNs.\n',k,walknum);
                        end
                    end
                end %ntp loop
                DT_np = cat(2,DT_np,d_np_trans);
                DT_np_idx = cat(1,DT_np_idx,d_np_trans_idx);
                DT_xs = cat(2,DT_xs,d_xs_trans);
                DT_xs_idx = cat(1,DT_xs_idx,d_xs_trans_idx);
                DT_xs_vchg = cat(1,DT_xs_vchg,d_xs_vchg); %used to filter trials based on percent change in velocity
                DT_gz = cat(2,DT_gz,d_gz_trans);
                DT_gz_idx = cat(1,DT_gz_idx,d_gz_trans_idx);
                DT_am = cat(2,DT_am,d_am_trans);
                DT_am_idx = cat(1,DT_am_idx,d_am_trans_idx);
                DT_kd = cat(2,DT_kd,d_kd_trans);
                DT_kd_idx = cat(1,DT_kd_idx,d_kd_trans_idx);
                DT_wv = cat(2,DT_wv,d_wv_trans);
                DT_wv_idx = cat(1,DT_wv_idx,d_wv_trans_idx);
                DT_bi = cat(2,DT_bi,d_bi_trans);
                DT_bi_idx = cat(1,DT_bi_idx,d_bi_trans_idx);
                DT_pe = cat(2,DT_pe,d_pe_trans);
                DT_pe_idx = cat(1,DT_pe_idx,d_pe_trans_idx);
                DT_im = cat(2,DT_im,d_im_trans);
                DT_im_idx = cat(1,DT_im_idx,d_im_trans_idx);
                EvntTrans = cat(1,EvntTrans,evnt);
                DescTrans = cat(1,DescTrans,desc);
                WalkNumTrans = cat(1,WalkNumTrans,ones(size(d_np_trans,2),1).*walknum);
                NumSampTrans = cat(1,NumSampTrans,ones(size(d_np_trans,2),1).*numsamp);

            end %for DTable (by walk)

            obj.DParsed.Trans.DT_np = DT_np;
            obj.DParsed.Trans.DT_np_idx = DT_np_idx;
            obj.DParsed.Trans.DT_xs = DT_xs;
            obj.DParsed.Trans.DT_xs_idx = DT_xs_idx;
            obj.DParsed.Trans.DT_xs_vchg = DT_xs_vchg;
            obj.DParsed.Trans.DT_gz = DT_gz;
            obj.DParsed.Trans.DT_gz_idx = DT_gz_idx;
            obj.DParsed.Trans.DT_am = DT_am;
            obj.DParsed.Trans.DT_am_idx = DT_am_idx;
            obj.DParsed.Trans.DT_kd = DT_kd;
            obj.DParsed.Trans.DT_kd_idx = DT_kd_idx;
            obj.DParsed.Trans.DT_wv = DT_wv;
            obj.DParsed.Trans.DT_wv_idx = DT_wv_idx;
            obj.DParsed.Trans.DT_bi = DT_bi;
            obj.DParsed.Trans.DT_bi_idx = DT_bi_idx;
            obj.DParsed.Trans.DT_pe = DT_pe;
            obj.DParsed.Trans.DT_pe_idx = DT_pe_idx;
            obj.DParsed.Trans.DT_im = DT_im;
            obj.DParsed.Trans.DT_im_idx = DT_im_idx;
            obj.DParsed.Trans.WT_np_samp = win_trans_np;
            obj.DParsed.Trans.WT_np_sec = win_trans_np/fs_np;
            obj.DParsed.Trans.WT_xs_samp = win_trans_xs;
            obj.DParsed.Trans.WT_xs_sec = win_trans_xs/fs_xs;
            obj.DParsed.Trans.WT_gz_samp = win_trans_gz;
            obj.DParsed.Trans.WT_gz_sec = win_trans_gz/fs_gz;
            obj.DParsed.Trans.WT_am_samp = win_trans_am;
            obj.DParsed.Trans.WT_am_sec = win_trans_am/fs_am;
            obj.DParsed.Trans.WT_kd_samp = win_trans_kd;
            obj.DParsed.Trans.WT_kd_sec = win_trans_kd/fs_kd;
            obj.DParsed.Trans.WT_wv_samp = win_trans_wv;
            obj.DParsed.Trans.WT_wv_sec = win_trans_wv/fs_wv;
            obj.DParsed.Trans.WT_bi_samp = win_trans_bi;
            obj.DParsed.Trans.WT_bi_sec = win_trans_bi/fs_bi;
            obj.DParsed.Trans.WT_pe_samp = win_trans_pe;
            obj.DParsed.Trans.WT_pe_sec = win_trans_pe/fs_pe;
            obj.DParsed.Trans.WT_im_samp = win_trans_im;
            obj.DParsed.Trans.WT_im_sec = win_trans_im/fs_im;
            obj.DParsed.Trans.FS_np = fs_np;
            obj.DParsed.Trans.FS_xs = fs_xs;
            obj.DParsed.Trans.FS_gz = fs_gz;
            obj.DParsed.Trans.FS_am = fs_am;
            obj.DParsed.Trans.FS_kd = fs_kd;
            obj.DParsed.Trans.FS_wv = fs_wv;
            obj.DParsed.Trans.F_wv = f_wv; %freq vector for wavelet
            obj.DParsed.Trans.FS_bi = fs_bi;
            obj.DParsed.Trans.F_bi = f_bi; %field names for biopac data
            obj.DParsed.Trans.FS_pe = fs_pe;
            obj.DParsed.Trans.F_pe = f_pe; %freq vector for pepisode
            obj.DParsed.Trans.FS_im = fs_im;
            obj.DParsed.Trans.F_im = f_im; %field names for biopac data
            obj.DParsed.Trans.EvntTrans = EvntTrans;
            obj.DParsed.Trans.DescTrans = DescTrans;
            obj.DParsed.Trans.WalkNumTrans = WalkNumTrans;
            obj.DParsed.Trans.NumSampTrans = NumSampTrans;
            obj.DParsed.Trans.Outliers = any(any(isnan(DT_np),1),3); %<-500 outliers handled in loadData
            obj.DParsed.Trans.WinTransSec = obj.WinTransSec;
            obj.DParsed.Seg.WinSegSec = obj.WinSegSec;

        end %parseTransData

        function getMultData(obj,varargin)
            %Generates multi-patient transition data for specgram analysis
            %and 2sec segment data for glme analysis.
            %getMultData;

            %%%%%%%%%%%%%%%%% Transitions %%%%%%%%%%%%%%%%
            DT_np = []; DT_np_idx = []; DT_xs_idx = []; DT_xs_vchg = [];
            DT_gz_idx = []; DT_am_idx = []; DT_kd_idx = []; DT_wv_idx = [];
            DT_bi_idx = []; DT_bi_flds = {}; DT_pe_idx = []; DT_im_idx = []; DT_im_flds = {};
            DT_Walk = []; DT_Evnt = []; DT_Desc = []; DT_StopWalk = []; 
            DT_GoWalk = []; DT_OL = []; DT_Region = []; DT_RegionLabel = []; DT_NSamp = []; %number of data samples in walk
            DT_Chan = []; DT_ChanLabel = []; DT_Patient = []; 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%% Segments %%%%%%%%%%%%%%%%%%%
            SegTable = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            MTable = cell(length(obj.PatientList),1); 
            for k=1:length(obj.PatientList)
                disp(k);
                if ~any(obj.PatientIdx==k)
                    obj.loadData('PatientIdx',k); %this now generates downsampled/normalized wavelet data
                end

                %%%%%%%%%%%%%%%%% Transitions %%%%%%%%%%%%%%%%
                obj.parseTransData;

                dt_np = obj.DParsed.Trans.DT_np; %time x trial x chan
                dt_np_idx = permute(repmat(obj.DParsed.Trans.DT_np_idx(:),1,4),[3,1,2]); %1 x trial x chan (center index of window in np data)
                dt_xs_idx = permute(repmat(obj.DParsed.Trans.DT_xs_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_xs_vchg = permute(repmat(obj.DParsed.Trans.DT_xs_vchg(:),1,4),[3,1,2]); %1 x trial x chan (used to filter trials by percentage velocity change)
                dt_gz_idx = permute(repmat(obj.DParsed.Trans.DT_gz_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_am_idx = permute(repmat(obj.DParsed.Trans.DT_am_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_kd_idx = permute(repmat(obj.DParsed.Trans.DT_kd_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_wv_idx = permute(repmat(obj.DParsed.Trans.DT_wv_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_bi_idx = permute(repmat(obj.DParsed.Trans.DT_bi_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_pe_idx = permute(repmat(obj.DParsed.Trans.DT_pe_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_im_idx = permute(repmat(obj.DParsed.Trans.DT_im_idx(:),1,4),[3,1,2]); %1 x trial x chan
                dt_walk = permute(repmat(obj.DParsed.Trans.WalkNumTrans(:),1,4),[3,1,2]); %1 x trial x chan
                dt_nsamp = permute(repmat(obj.DParsed.Trans.NumSampTrans(:),1,4),[3,1,2]); %1 x trial x chan
                dt_evnt = permute(repmat(obj.DParsed.Trans.EvntTrans(:),1,4),[3,1,2]); %1 x trial x chan
                dt_desc = permute(repmat(obj.DParsed.Trans.DescTrans(:),1,4),[3,1,2]); %1 x trial x chan
                dt_stopwalk = permute(repmat(obj.DParsed.Trans.WalkNumTrans(:)==obj.StopGoWalks{k}(1),1,4),[3,1,2]); %1 x trial x chan
                dt_gowalk = permute(repmat(obj.DParsed.Trans.WalkNumTrans(:)==obj.StopGoWalks{k}(2),1,4),[3,1,2]); %1 x trial x chan
                dt_ol = permute(repmat(obj.DParsed.Trans.Outliers(:),1,4),[3,1,2]); %nan and -500 only (1 x trial x chan)
                dt_regiontable = obj.RegionTable(obj.RegionTable.Patient==k,:);
                dt_region = permute(repmat(dt_regiontable.Region(:)',size(dt_np,2),1),[3,1,2]); %1 x trial x chan
                dt_regionlabel = permute(repmat(dt_regiontable.RegionLabel(:)',size(dt_np,2),1),[3,1,2]); %1 x trial x chan
                dt_chan = permute(repmat(dt_regiontable.Chan(:)',size(dt_np,2),1),[3,1,2]); %1 x trial x chan
                dt_chanlabel = permute(repmat(dt_regiontable.ChanLabel(:)',size(dt_np,2),1),[3,1,2]); %1 x trial x chan
                dt_pt = permute(repmat(k,size(dt_np,2),4),[3,1,2]); %1 x trial x chan

                DT_np = cat(2,DT_np,dt_np);
                DT_np_idx = cat(2,DT_np_idx,dt_np_idx);
                DT_xs_idx = cat(2,DT_xs_idx,dt_xs_idx);
                DT_xs_vchg = cat(2,DT_xs_vchg,dt_xs_vchg);
                DT_gz_idx = cat(2,DT_gz_idx,dt_gz_idx);
                DT_am_idx = cat(2,DT_am_idx,dt_am_idx);
                DT_kd_idx = cat(2,DT_kd_idx,dt_kd_idx);
                DT_wv_idx = cat(2,DT_wv_idx,dt_wv_idx);
                DT_bi_idx = cat(2,DT_bi_idx,dt_bi_idx);
                DT_bi_flds = cat(1,DT_bi_flds,obj.DParsed.Trans.F_bi);
                DT_pe_idx = cat(2,DT_pe_idx,dt_pe_idx);
                DT_im_idx = cat(2,DT_im_idx,dt_im_idx);
                DT_im_flds = cat(1,DT_im_flds,obj.DParsed.Trans.F_im);
                DT_Walk = cat(2,DT_Walk,dt_walk);
                DT_NSamp = cat(2,DT_NSamp,dt_nsamp); %number of data samples in walk
                DT_Evnt = cat(2,DT_Evnt,dt_evnt);
                DT_Desc = cat(2,DT_Desc,dt_desc);
                DT_StopWalk = cat(2,DT_StopWalk,dt_stopwalk);
                DT_GoWalk = cat(2,DT_GoWalk,dt_gowalk);
                DT_OL = cat(2,DT_OL,dt_ol);
                DT_Region = cat(2,DT_Region,dt_region);
                DT_RegionLabel = cat(2,DT_RegionLabel,dt_regionlabel);
                DT_Chan = cat(2,DT_Chan,dt_chan);
                DT_ChanLabel = cat(2,DT_ChanLabel,dt_chanlabel);
                DT_Patient = cat(2,DT_Patient,dt_pt);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%% Segments %%%%%%%%%%%%%%%%%%%
                obj.parseSegData;

                DS_np = obj.DParsed.Seg.DS_np;
                DS_xs = obj.DParsed.Seg.DS_xs;
                DS_gz = obj.DParsed.Seg.DS_gz; %boolean for fixation/no fixation
                DS_am = obj.DParsed.Seg.DS_am;
                DS_kd = obj.DParsed.Seg.DS_kd;
                DS_wv = obj.DParsed.Seg.DS_wv;
                DS_ed = obj.DParsed.Seg.DS_bi(:,:,contains(obj.DParsed.Seg.F_bi,'PhasicEDA'));
                DS_pe = obj.DParsed.Seg.DS_pe;
                DS_gy = obj.DParsed.Seg.DS_im(:,:,contains(obj.DParsed.Seg.F_im,'gyroY'));
                DS_wk = obj.DParsed.Seg.WalkNumSeg;
                DS_in = obj.DParsed.Seg.InFlag;

                DS_ol = obj.DParsed.Seg.Outliers; %these have nans for <-500 condition
                DS_np(:,DS_ol,:) = [];
                DS_xs(:,DS_ol) = [];
                DS_gz(:,DS_ol) = [];
                DS_am(:,DS_ol) = [];
                DS_kd(:,DS_ol) = [];
                DS_wv(:,DS_ol,:,:) = [];
                DS_pe(:,DS_ol,:,:) = [];
                DS_ed(:,DS_ol) = [];
                DS_gy(:,DS_ol) = [];
                DS_wk(DS_ol) = [];
                DS_in(DS_ol) = [];

                DS_stopwk = (DS_wk == obj.StopGoWalks{k}(1));
                DS_gowk = (DS_wk == obj.StopGoWalks{k}(2));

                %raw data from a sensor, so need to reject outliers
                DS_xs(:,any(DS_xs>5)) = nan;
                ds_xs = median(DS_xs)';
                ds_xs_norm = ds_xs;
                idx = find(isoutlier(ds_xs_norm)); 
                idx2 = ds_xs(idx)>median(ds_xs,'omitnan'); %only discard outliers if they are greater than the median walking speed (i.e. don't want to discard stops)
                % ds_xs(idx(idx2)) = nan;
                ds_xs_norm(idx(idx2)) = nan;

                %outliers should not exist since it is a percentage (normalization does not seem to work well for this data)
                ds_gz = mean(DS_gz)'; %must be mean to calculate percentage time fixating
                ds_gz_norm = ds_gz;
                % ds_gz(ds_gz==0) = min(ds_gz(ds_gz>0)); %a way to remove tails of the distribution but not sure if this is correct
                % ds_gz(ds_gz==1) = max(ds_gz(ds_gz<1));
                % ds_gz_norm = norminv(ds_gz); 
                % ds_gz_norm = atanh((ds_gz-0.5)*2); %basically does the same as norminv

                %raw data from a sensor, so need to reject outliers
                ds_am = median(DS_am)';
                ds_am_norm = log10(ds_am+eps);
                % idx = isoutlier(ds_am_norm); %no values are removed here
                % ds_am(idx) = nan;
                % ds_am_norm(idx) = nan;

                %outliers should not exist since this is a calculated value
                ds_kd = median(DS_kd)';
                ds_kd_norm = ds_kd;
                
                %this is already processed data (zscore? outliers?)
                DS_ed(:,any(abs(DS_ed)>3)) = nan; %this appears to be zscore data, so discarding >3std since artifact does appear to occur in this data
                ds_ed = median(DS_ed)';
                ds_ed_norm = ds_ed;
                % idx = isoutlier(ds_ed_norm,'ThresholdFactor',10); %most values are centered at zero so high threshold is needed (still seems to throw away actual data - disabling for now)
                % ds_ed(idx) = nan;
                % ds_ed_norm(idx) = nan;

                %IMU gyroY (outliers might be important since this
                %represents head turns and mean/median might not be the
                %best approch for head turns in 2sec segment)
                ds_gy = mean(DS_gy)'; %mean might capture a quick head turn in a 2sec segment better than median
                ds_gy_norm = ds_gy;
                % idx = isoutlier(ds_gy_norm); %no values are removed here
                % ds_gy(idx) = nan;
                % ds_gy_norm(idx) = nan;

                %freq bands of interest
                ds_f_spec = obj.DParsed.Seg.F_wv;
                ds_f_band = [4,8;8,12;12,30;30,70;70,120]; %Theta,Alpha,Beta,Gamma,HG
                ds_nband = size(ds_f_band,1);
                ds_idx60 = (ds_f_spec>=59 & ds_f_spec<=61);
                ds_fidx = false(length(ds_f_spec),5);
                for n=1:ds_nband
                    ds_fidx(:,n) = (ds_f_spec>=ds_f_band(n,1) & ds_f_spec<ds_f_band(n,2)) & ~ds_idx60;
                end

                %multitaper params
                ds_params.Fs = obj.DTable.fs_np(1); %250
                ds_params.trialave = 0;
                ds_params.tapers = [5,9];
                ds_params.fpass = [2,120];
                ds_params.pad = 0;
                [~,ds_f_spec_mt] = mtspectrumc(DS_np(:,:,1),ds_params);
                ds_idx60_mt = (ds_f_spec_mt>=59 & ds_f_spec_mt<=61);
                ds_fidx_mt = false(length(ds_f_spec_mt),5);
                for n=1:ds_nband
                    ds_fidx_mt(:,n) = (ds_f_spec_mt>=ds_f_band(n,1) & ds_f_spec_mt<ds_f_band(n,2)) & ~ds_idx60_mt;
                end

                %fooof params
                ds_settings.max_n_peaks = 3; %aperiodic component has a better fit if more peaks are removed
                ds_settings.aperiodic_mode = 'fixed';
                ds_settings.peak_width_limits = [1,12];
                ds_flabels = {'peak_freq','peak_height','peak_width','aperiodic_offset','aperiodic_exponent','fit_rsquared','fit_error'};

                ds_nseg = size(DS_np,2); %number of segments after outlier removal
                ds_regiontable = obj.RegionTable(obj.RegionTable.Patient==k,:);
                for m=1:4 %by channel
                    % ied_ol = isoutlier(max(abs(DS_np(:,:,m))))'; %simple ied detection

                    ds_S = mtspectrumc(DS_np(:,:,m),ds_params); ds_S = ds_S'; %nseg x freq

                    %fooof
                    ds_FR = nan(ds_nseg,length(ds_flabels));
                    parfor n=1:ds_nseg
                        ds_fr = fooof(ds_f_spec_mt, ds_S(n,:), [ds_f_spec_mt(1),85], ds_settings, 1); %fr.peak_params = [freq, height (aperiodic removed), width] (85Hz matches upper lim in plotMultTransSpecGramPerm)
                        if isempty(ds_fr.peak_params)
                            ds_pkp = nan(1,3);
                        else
                            [~,midx] = max(ds_fr.peak_params(:,2)); %find the highest peak
                            ds_pkp = ds_fr.peak_params(midx,:);
                        end
                        ds_FR(n,:) = [ds_pkp,ds_fr.aperiodic_params,ds_fr.r_squared,ds_fr.error];
                    end

                    ds_wv = nan(ds_nseg,ds_nband);
                    ds_wv_norm = nan(ds_nseg,ds_nband);
                    ds_mt = nan(ds_nseg,ds_nband); %multitaper
                    ds_mt_norm = nan(ds_nseg,ds_nband);
                    for n=1:ds_nband
                        ds_wv(:,n) = median(mean(DS_wv(:,:,ds_fidx(:,n),m),3),1)'; %mean across band, median across time
                        ds_wv_norm(:,n) = log10(ds_wv(:,n)+eps);
                        idx = isoutlier(ds_wv_norm(:,n));
                        ds_wv(idx,n) = nan;
                        ds_wv_norm(idx,n) = nan;

                        ds_mt(:,n) = mean(ds_S(:,ds_fidx_mt(:,n)),2);
                        ds_mt_norm(:,n) = log10(ds_mt(:,n)+eps);
                        idx = isoutlier(ds_mt_norm(:,n));
                        ds_mt(idx,n) = nan;
                        ds_mt_norm(idx,n) = nan;
                    end

                    %bosc pepisode (must be mean to calculate percentage since data is binary)
                    ds_pe = squeeze(mean(DS_pe(:,:,:,m))); %segment x band (theta, alpha, beta)
                    ds_pe_norm = ds_pe;
               
                    ds_pt = repmat(k,ds_nseg,1);
                    ds_ch = repmat(m,ds_nseg,1);

                    ds_rgnum = repmat(ds_regiontable.Region(m),ds_nseg,1);
                    ds_rglbl = repmat(ds_regiontable.RegionLabel(m),ds_nseg,1);
                    ds_chlbl = repmat(ds_regiontable.ChanLabel(m),ds_nseg,1);

                    %Only include normalized data for now. Amb is the only
                    %predictor with normalization. Power should be log10
                    %normalized by default to make more gaussian.
                    segtable = table(...
                        ds_pt,ds_ch,DS_wk,DS_stopwk,DS_gowk,DS_in,ds_rgnum,ds_rglbl,ds_chlbl,...
                        ds_xs_norm,ds_gz_norm,ds_am_norm,ds_kd_norm,ds_ed_norm,ds_gy_norm,...
                        ds_wv_norm(:,1),ds_wv_norm(:,2),ds_wv_norm(:,3),ds_wv_norm(:,4),ds_wv_norm(:,5),...
                        ds_mt_norm(:,1),ds_mt_norm(:,2),ds_mt_norm(:,3),ds_mt_norm(:,4),ds_mt_norm(:,5),...
                        ds_pe_norm(:,1),ds_pe_norm(:,2),ds_pe_norm(:,3),...
                        ds_FR(:,1),ds_FR(:,2),ds_FR(:,3),ds_FR(:,4),ds_FR(:,5),ds_FR(:,6),ds_FR(:,7),...
                        'VariableNames',...
                        {'Pt','Chan','Walk','StopWalk','GoWalk','InFlag','RegionNum','RegionLabel','ChanLabel',...
                        'Vel','Fix','Amb','KDE','EDA','HeadTurn',...
                        'wvTheta','wvAlpha','wvBeta','wvGamma','wvHG',...
                        'mtTheta','mtAlpha','mtBeta','mtGamma','mtHG',...
                        'peTheta','peAlpha','peBeta',...
                        ds_flabels{1},ds_flabels{2},ds_flabels{3},ds_flabels{4},ds_flabels{5},ds_flabels{6},ds_flabels{7},...
                        });

                    SegTable = vertcat(SegTable,segtable);
                end %chan
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                MTable{k} = obj.DTable;
            end %patient

            %%%%%%%%%%%%%%%%% Transitions %%%%%%%%%%%%%%%%
            DT_np = reshape(DT_np,size(DT_np,1),[]); %data (unwrapped by time, then trial, then chan for all patients - need this for OL2 calculation)
            DT_np_idx = reshape(DT_np_idx,1,[]);
            DT_xs_idx = reshape(DT_xs_idx,1,[]);
            DT_xs_vchg = reshape(DT_xs_vchg,1,[]);
            DT_gz_idx = reshape(DT_gz_idx,1,[]);
            DT_am_idx = reshape(DT_am_idx,1,[]);
            DT_kd_idx = reshape(DT_kd_idx,1,[]);
            DT_wv_idx = reshape(DT_wv_idx,1,[]);
            DT_bi_idx = reshape(DT_bi_idx,1,[]);
            DT_pe_idx = reshape(DT_pe_idx,1,[]);
            DT_im_idx = reshape(DT_im_idx,1,[]);
            DT_Walk = reshape(DT_Walk,1,[]);
            DT_NSamp = reshape(DT_NSamp,1,[]); %number of data samples in walk
            DT_Evnt = reshape(DT_Evnt,1,[]);
            DT_Desc = reshape(DT_Desc,1,[]);
            DT_StopWalk = reshape(DT_StopWalk,1,[]);
            DT_GoWalk = reshape(DT_GoWalk,1,[]);
            DT_OL = reshape(DT_OL,1,[]);
            DT_Region = reshape(DT_Region,1,[]);
            DT_RegionLabel = reshape(DT_RegionLabel,1,[]);
            DT_Chan = reshape(DT_Chan,1,[]);
            DT_ChanLabel = reshape(DT_ChanLabel,1,[]);
            DT_Patient = reshape(DT_Patient,1,[]);
            
            obj.MultTrans.DT_np_idx = DT_np_idx;
            obj.MultTrans.DT_xs_idx = DT_xs_idx;
            obj.MultTrans.DT_xs_vchg = DT_xs_vchg;
            obj.MultTrans.DT_gz_idx = DT_gz_idx;
            obj.MultTrans.DT_am_idx = DT_am_idx;
            obj.MultTrans.DT_kd_idx = DT_kd_idx;
            obj.MultTrans.DT_wv_idx = DT_wv_idx;
            obj.MultTrans.DT_bi_idx = DT_bi_idx;
            obj.MultTrans.DT_pe_idx = DT_pe_idx;
            obj.MultTrans.DT_im_idx = DT_im_idx;
            obj.MultTrans.Walk = DT_Walk;
            obj.MultTrans.NSamp = DT_NSamp;
            obj.MultTrans.Evnt = DT_Evnt;
            obj.MultTrans.Desc = DT_Desc;
            obj.MultTrans.StopWalk = DT_StopWalk;
            obj.MultTrans.GoWalk = DT_GoWalk;
            obj.MultTrans.OL = DT_OL; %only includes <-500 and NaNs (different than any(isnan(DT_np),1) below since all 4 chans are included in the detection in parseTransData)
            obj.MultTrans.OL2 = any(isnan(DT_np),1)|isoutlier(max(abs(zscore(DT_np)))); %NaNs plus simple IED detection (combined across all pts/chans -> should do this by pt/chan not combined)
            obj.MultTrans.Region = DT_Region;
            obj.MultTrans.RegionLabel = DT_RegionLabel;
            obj.MultTrans.Chan = DT_Chan;
            obj.MultTrans.ChanLabel = DT_ChanLabel;
            obj.MultTrans.Patient = DT_Patient;
            obj.MultTrans.TimeSamp_np = obj.DParsed.Trans.WT_np_samp;
            obj.MultTrans.TimeSec_np = obj.DParsed.Trans.WT_np_sec;
            obj.MultTrans.TimeSamp_xs = obj.DParsed.Trans.WT_xs_samp;
            obj.MultTrans.TimeSec_xs = obj.DParsed.Trans.WT_xs_sec;
            obj.MultTrans.TimeSamp_gz = obj.DParsed.Trans.WT_gz_samp;
            obj.MultTrans.TimeSec_gz = obj.DParsed.Trans.WT_gz_sec;
            obj.MultTrans.TimeSamp_am = obj.DParsed.Trans.WT_am_samp;
            obj.MultTrans.TimeSec_am = obj.DParsed.Trans.WT_am_sec;
            obj.MultTrans.TimeSamp_kd = obj.DParsed.Trans.WT_kd_samp;
            obj.MultTrans.TimeSec_kd = obj.DParsed.Trans.WT_kd_sec;
            obj.MultTrans.TimeSamp_wv = obj.DParsed.Trans.WT_wv_samp;
            obj.MultTrans.TimeSec_wv = obj.DParsed.Trans.WT_wv_sec;
            obj.MultTrans.TimeSamp_bi = obj.DParsed.Trans.WT_bi_samp;
            obj.MultTrans.TimeSec_bi = obj.DParsed.Trans.WT_bi_sec;
            obj.MultTrans.TimeSamp_pe = obj.DParsed.Trans.WT_pe_samp;
            obj.MultTrans.TimeSec_pe = obj.DParsed.Trans.WT_pe_sec;
            obj.MultTrans.TimeSamp_im = obj.DParsed.Trans.WT_im_samp;
            obj.MultTrans.TimeSec_im = obj.DParsed.Trans.WT_im_sec;
            obj.MultTrans.FS_np = obj.DParsed.Trans.FS_np;
            obj.MultTrans.FS_xs = obj.DParsed.Trans.FS_xs;
            obj.MultTrans.FS_gz = obj.DParsed.Trans.FS_gz;
            obj.MultTrans.FS_am = obj.DParsed.Trans.FS_am;
            obj.MultTrans.FS_kd = obj.DParsed.Trans.FS_kd;
            obj.MultTrans.FS_wv = obj.DParsed.Trans.FS_wv;
            obj.MultTrans.FS_bi = obj.DParsed.Trans.FS_bi;
            obj.MultTrans.FS_pe = obj.DParsed.Trans.FS_pe;
            obj.MultTrans.FS_im = obj.DParsed.Trans.FS_im;
            obj.MultTrans.F_bi = DT_bi_flds;
            obj.MultTrans.F_im = DT_im_flds;
            obj.MultTrans.F_pe = obj.DParsed.Trans.F_pe; %freq bands for pepisode
            obj.MultTrans.Freq = obj.DParsed.Trans.F_wv; %freq vector (same for all datasets)
            obj.MultTrans.WinSegSec = obj.WinSegSec; %Also smoothing kernel for specgrams in trans analysis
            obj.MultTrans.WinTransSec = obj.WinTransSec;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%% Segments %%%%%%%%%%%%%%%%%%%
            %Table of individual 2sec segments for all data streams
            %(xs/fix/amb), NaNs exist for xs/fix/amp. Data normalized for
            %each stream. However, according to chatgpt, only the response
            %variable in a glme needs to be normal (or distribution type
            %needs to be specified). The fixed effects predictors can be
            %non-normal.
            [uPtChan,~,uPtChanIdx] = unique(SegTable(:,{'Pt','Chan'}),'rows');
            [uPtWalk,~,uPtWalkIdx] = unique(SegTable(:,{'Pt','Walk'}),'rows');
            SegTable = addvars(SegTable,uPtChanIdx,'NewVariableNames','uPtChan','Before','Pt');
            SegTable = addvars(SegTable,uPtWalkIdx,'NewVariableNames','uPtWalk','Before','Pt');

            %%%%%%%%%%%%%% Checking data %%%%%%%%%%%%%%%
            % flbls = {'peTheta','peAlpha','peBeta'};
            % for m=1:length(flbls)
            %     fH = figure;
            %     fH.Name = flbls{m};
            %     for k=1:size(uPtChan,1)
            %         aH = subplot(4,5,k,'parent',fH);
            %         idx = find(uPtChanIdx==k);
            %         d = SegTable.(flbls{m})(idx);
            %         n = sum(isnan(d));
            %         histogram(d,'Parent',aH);
            %         title(aH,sprintf('P%0.0f,C%0.0f (%0.0f/%0.0f nan)',SegTable.Pt(idx(1)),SegTable.Chan(idx(1)),n,numel(d)))
            %     end
            % end
            % 
            % fH = figure;
            % fH.Name = 'HeadTurn';
            % for k=1:size(uPtChan,1)
            %     aH = subplot(4,5,k,'parent',fH);
            %     idx = find(uPtChanIdx==k);
            %     d = SegTable.(fH.Name)(idx);
            %     n = sum(isnan(d));
            %     histogram(d,'Parent',aH);
            %     title(aH,sprintf('P%0.0f,C%0.0f (%0.0f/%0.0f nan)',SegTable.Pt(idx(1)),SegTable.Chan(idx(1)),n,numel(d)))
            % end

            obj.MultSeg.SegTable = SegTable;
            obj.MultSeg.uPtChan = uPtChan;
            obj.MultSeg.uPtChanIdx = uPtChanIdx;
            obj.MultSeg.uPtWalk = uPtWalk;
            obj.MultSeg.uPtWalkIdx = uPtWalkIdx;
            obj.MultSeg.f = ds_f_spec;
            obj.MultSeg.f_mt = ds_f_spec_mt;
            obj.MultSeg.f_band = ds_f_band;
            obj.MultSeg.f_band_pe = obj.DParsed.Seg.F_pe;
            obj.MultSeg.f_lbls = {'Theta','Alpha','Beta','Gamma','HG'}; %freq band labels
            obj.MultSeg.f_lbls_pe = {'Theta','Alpha','Beta'}; %freq band labels
            obj.MultSeg.p_lbls = {'Vel','Fix','Amb','KDE','EDA'}; %predictor labels
            obj.MultSeg.WinSegSec = obj.WinSegSec; %Also smoothing kernel for specgrams in trans analysis
            obj.MultSeg.WinTransSec = obj.WinTransSec;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            obj.MultTable = MTable; %this can also be filled in getMultSegData
        end %getMultData
        
        function filterMultTransData(obj,varargin)
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
            rmoverlap = p.rmoverlap; %remove any overlapping trials

            if isempty(obj.MultTrans)
                obj.getMultData;
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
            transidx = transidx & ~obj.MultTrans.OL2; %calculated in getMultTransData

            np_idx = obj.MultTrans.DT_np_idx(transidx);
            xs_idx = obj.MultTrans.DT_xs_idx(transidx);
            xs_vchg = obj.MultTrans.DT_xs_vchg(transidx); %percent velocity change
            gz_idx = obj.MultTrans.DT_gz_idx(transidx);
            am_idx = obj.MultTrans.DT_am_idx(transidx);
            kd_idx = obj.MultTrans.DT_kd_idx(transidx);
            wv_idx = obj.MultTrans.DT_wv_idx(transidx);
            bi_idx = obj.MultTrans.DT_bi_idx(transidx);
            pe_idx = obj.MultTrans.DT_pe_idx(transidx);
            im_idx = obj.MultTrans.DT_im_idx(transidx);
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
            MT = table(pt',wk',swk',gwk',rg',lb',ch',np_idx',ns',ev',ds',ds_idx',...
                xs_idx',xs_vchg',gz_idx',am_idx',kd_idx',wv_idx',bi_idx',pe_idx',im_idx',...
                'VariableNames',{'patient','walk','stopwalk','gowalk','regionnum',...
                'regionlabel','chan','npidx','nsamp','evnt','desc','descidx',...
                'xsidx','xsvchg','gzidx','amidx','kdidx','wvidx','biidx','peidx','imidx'});

            %filter by region (amygdala/anterior hipp)
            if all(regionidx<5)
                idx = ~ismember(MT.regionnum,regionidx); %1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
                MT(idx,:) = [];
            elseif any(regionidx==6) %custom regions (nx2 matrix where cols are patient,chan)
                idx = false(size(MT,1),1);
                for k=1:size(customregion,1)
                    idx = idx | all(MT.patient==customregion(k,1) & MT.chan==customregion(k,2),2);
                end
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
                MT(idx,:) = [];
            end

            %filter by patient
            if patientidx==-1
                idx = ~ismember(MT.patient,str2num(patienttype)); %filter out all but specified walk number
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
                MT(idx,:) = [];
            end

            %filter by description
            idx = ~MT.descidx;
            MT(idx,:) = [];

            %filter by velocity
            if ~isempty(veltype)
                idx = isnan(MT.xsvchg);
                MT(idx,:) = [];

                [uPt,~,uPtIdx] = unique(MT.patient);
                mtidx = [];
                for m=1:length(uPt)
                    ptidx = find(uPtIdx==m);
                    mt = MT(ptidx,:);

                    [~,sidx] = sort(abs(mt.xsvchg),'descend'); %larger values represent largest change from median (sorted 1 to 0)

                    mid_len = floor(length(sidx)/2); %middle
                    ter_len = floor(length(sidx)/3); %tercile
                    ten_len = floor(length(sidx)*0.1); %10th percentile
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
                        case 'VelHigh10'
                            ten_idx = 1:ten_len;
                            idx = sidx(setdiff(full_idx,ten_idx));
                        case 'VelLow10'
                            ten_idx = full_len-ten_len+1:full_len;
                            idx = sidx(setdiff(full_idx,ten_idx));
                        otherwise
                            error('veltype is incorrect!')
                    end
                    mtidx = cat(1,mtidx,ptidx(idx));
                end %by pt
                
                MT(mtidx,:) = [];
            end %vel

            %filter out overlap
            if rmoverlap
                [uPWX,~,uPWXIdx] = unique(MT(:,{'patient','walk','npidx'})); %need to remove trials that have the same dataidx for a given patient/walk
                [uPW,~,uPWIdx] = unique(uPWX(:,{'patient','walk'}));

                ovrlp = diff(p.transrng); %overlap window in sec
                RMIdx = false(length(uPWXIdx),1);
                for k=1:size(uPW,1)
                    idx = find(k==uPWIdx); %index into uPWX
                    [snp,sidx] = sort(uPWX.npidx(idx));

                    %finding overlap trials and removing
                    m=1;
                    while m<length(snp)
                        dnp = (snp(m+1)-snp(m))/250; %difference in sec
                        if dnp<ovrlp
                            snp(m+1) = [];
                            sidx(m+1) = [];
                        else
                            m=m+1;
                        end
                    end
                    rmidx = setdiff(idx,idx(sidx)); %indices to remove from uPWX
                    RMIdx = RMIdx | ismember(uPWXIdx,rmidx);                  
                end
                MT(RMIdx,:) = [];
                % [PercOverlap,AvgOverlapSec] = obj.calcOverlap(MT(~RMIdx,:),p.transrng); %did it work?
            end

            obj.MultTrans.MT = MT;

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

        function getFilteredMultTransData(obj,transrng,varargin)
            %MT is generated by filterMultTransData. transrng is the
            %start/stop time in sec for the window with zero at center.
            %tsamp is a vector of samples centered at zero specifying the
            %window size.

            tsamp_np = round(transrng(1)*obj.MultTrans.FS_np):round(transrng(2)*obj.MultTrans.FS_np);
            tsec_np = tsamp_np./obj.MultTrans.FS_np;
            ntime_np = length(tsamp_np);

            tsamp_xs = round(transrng(1)*obj.MultTrans.FS_xs):round(transrng(2)*obj.MultTrans.FS_xs);
            tsec_xs = tsamp_xs./obj.MultTrans.FS_xs;
            ntime_xs = length(tsamp_xs);

            tsamp_gz = round(transrng(1)*obj.MultTrans.FS_gz):round(transrng(2)*obj.MultTrans.FS_gz);
            tsec_gz = tsamp_gz./obj.MultTrans.FS_gz;
            ntime_gz = length(tsamp_gz);

            tsamp_am = round(transrng(1)*obj.MultTrans.FS_am):round(transrng(2)*obj.MultTrans.FS_am);
            tsec_am = tsamp_am./obj.MultTrans.FS_am;
            ntime_am = length(tsamp_am);

            tsamp_kd = round(transrng(1)*obj.MultTrans.FS_kd):round(transrng(2)*obj.MultTrans.FS_kd);
            tsec_kd = tsamp_kd./obj.MultTrans.FS_kd;
            ntime_kd = length(tsamp_kd);

            tsamp_wv = round(transrng(1)*obj.MultTrans.FS_wv):round(transrng(2)*obj.MultTrans.FS_wv); 
            tsec_wv = tsamp_wv./obj.MultTrans.FS_wv;
            freq_wv = obj.MultTrans.Freq;
            ntime_wv = length(tsamp_wv);
            nfreq_wv = length(freq_wv);

            tsamp_bi = round(transrng(1)*obj.MultTrans.FS_bi):round(transrng(2)*obj.MultTrans.FS_bi);
            tsec_bi = tsamp_bi./obj.MultTrans.FS_bi;
            ntime_bi = length(tsamp_bi);
            flds_bi = obj.MultTrans.F_bi;
            nflds_bi = size(flds_bi,2); %always 8 fields of data (different order for pt5)

            tsamp_pe = round(transrng(1)*obj.MultTrans.FS_pe):round(transrng(2)*obj.MultTrans.FS_pe); 
            tsec_pe = tsamp_pe./obj.MultTrans.FS_pe;
            freq_pe = obj.MultTrans.F_pe;
            ntime_pe = length(tsamp_pe);
            nfreq_pe = size(freq_pe,1);

            tsamp_im = round(transrng(1)*obj.MultTrans.FS_im):round(transrng(2)*obj.MultTrans.FS_im);
            tsec_im = tsamp_im./obj.MultTrans.FS_im;
            ntime_im = length(tsamp_im);
            flds_im = obj.MultTrans.F_im;
            nflds_im = size(flds_im,2); %8 fields of data
            
            [uPtWkCh,~,uPtWkChIdx] = unique(obj.MultTrans.MT(:,{'patient','walk','chan'}));
            
            d_np = nan(size(uPtWkChIdx,1),ntime_np); %raw np data for running fooof, etc.
            d_xs = nan(size(uPtWkChIdx,1),ntime_xs);
            d_gz = nan(size(uPtWkChIdx,1),ntime_gz);
            d_am = nan(size(uPtWkChIdx,1),ntime_am);
            d_kd = nan(size(uPtWkChIdx,1),ntime_kd);
            d_wv = nan(size(uPtWkChIdx,1),ntime_wv,nfreq_wv);
            d_ed = nan(size(uPtWkChIdx,1),ntime_bi);
            d_pe = nan(size(uPtWkChIdx,1),ntime_pe,nfreq_pe);
            d_gy = nan(size(uPtWkChIdx,1),ntime_im);
            for k=1:size(uPtWkCh,1)
                idx = find(k==uPtWkChIdx);
                mt = obj.MultTrans.MT(idx,:);
                pt = uPtWkCh.patient(k);
                wk = uPtWkCh.walk(k);
                ch = uPtWkCh.chan(k);

                npidx = round(mt.npidx)+tsamp_np; 
                xsidx = round(mt.xsidx)+tsamp_xs; 
                gzidx = round(mt.gzidx)+tsamp_gz; 
                amidx = round(mt.amidx)+tsamp_am; 
                kdidx = round(mt.kdidx)+tsamp_kd; 
                wvidx = round(mt.wvidx)+tsamp_wv;
                biidx = round(mt.biidx)+tsamp_bi;
                peidx = round(mt.peidx)+tsamp_pe;
                imidx = round(mt.imidx)+tsamp_im;

                %NP
                if ~isempty(obj.MultTable{pt}.d_np{wk})
                    d = obj.MultTable{pt}.d_np{wk}(:,ch); %time x chan
                    ridx = any(npidx<1,2)|any(npidx>length(d),2)|any(isnan(mt.npidx),2); %remove index
                    npidx(ridx,:) = [];
                    d = reshape(d(npidx,:),size(npidx,1),size(npidx,2));
                    d_np(idx(~ridx),:) = d; %trial x time
                end

                %XS
                if ~isempty(obj.MultTable{pt}.d_xs{wk})
                    d = obj.MultTable{pt}.d_xs{wk}; 
                    ridx = any(xsidx<1,2)|any(xsidx>length(d),2)|any(isnan(mt.xsidx),2); %remove index
                    xsidx(ridx,:) = [];
                    d = reshape(d(xsidx,:),size(xsidx,1),size(xsidx,2));
                    d_xs(idx(~ridx),:) = d; %trial x time
                end

                %Gaze
                if ~isempty(obj.MultTable{pt}.d_gaze_fix{wk})
                    d = obj.MultTable{pt}.d_gaze_fix{wk}; 
                    ridx = any(gzidx<1,2)|any(gzidx>length(d),2)|any(isnan(mt.gzidx),2); %remove index
                    gzidx(ridx,:) = [];
                    d = reshape(d(gzidx,:),size(gzidx,1),size(gzidx,2));
                    d_gz(idx(~ridx),:) = d; %trial x time
                end

                %Amb
                if ~isempty(obj.MultTable{pt}.d_amb{wk})
                    d = obj.MultTable{pt}.d_amb{wk}; 
                    ridx = any(amidx<1,2)|any(amidx>length(d),2)|any(isnan(mt.amidx),2); %remove index
                    amidx(ridx,:) = [];
                    d = reshape(d(amidx,:),size(amidx,1),size(amidx,2));
                    d_am(idx(~ridx),:) = d; %trial x time
                end

                %KDE
                if ~isempty(obj.MultTable{pt}.d_kde{wk})
                    d = obj.MultTable{pt}.d_kde{wk}; %time x chan
                    ridx = any(kdidx<1,2)|any(kdidx>length(d),2)|any(isnan(mt.kdidx),2); %remove index
                    kdidx(ridx,:) = [];
                    d = reshape(d(kdidx,:),size(kdidx,1),size(kdidx,2));
                    d_kd(idx(~ridx),:) = d; %trial x time
                end

                %Wav
                if ~isempty(obj.MultTable{pt}.d_wav{wk})
                    d = obj.MultTable{pt}.d_wav{wk}(:,:,ch); %time x freq x chan
                    ridx = any(wvidx<1,2)|any(wvidx>size(d,1),2)|any(isnan(mt.wvidx),2); %remove index
                    wvidx(ridx,:) = [];
                    d = reshape(d(wvidx,:),size(wvidx,1),size(wvidx,2),nfreq_wv);
                    d_wv(idx(~ridx),:,:) = d; %trial x time x freq
                end

                %Biopac
                if ~isempty(obj.MultTable{pt}.d_bio{wk})
                    d = table2array(obj.MultTable{pt}.d_bio{wk}(:,contains(obj.MultTrans.F_bi(pt,:),'PhasicEDA'))); %time x field (data type)
                    ridx = any(biidx<1,2)|any(biidx>size(d,1),2)|any(isnan(mt.biidx),2); %remove index
                    biidx(ridx,:) = [];
                    d = reshape(d(biidx,:),size(biidx,1),size(biidx,2));
                    d_ed(idx(~ridx),:) = d; %trial x time
                end

                %Bosc Pepisode
                if ~isempty(obj.MultTable{pt}.d_pep{wk})
                    d = obj.MultTable{pt}.d_pep{wk}(:,:,ch); %time x freq x chan
                    ridx = any(peidx<1,2)|any(peidx>size(d,1),2)|any(isnan(mt.peidx),2); %remove index
                    peidx(ridx,:) = [];
                    d = reshape(d(peidx,:),size(peidx,1),size(peidx,2),nfreq_pe);
                    d_pe(idx(~ridx),:,:) = d; %trial x time x freq
                end

                %IMU
                if ~isempty(obj.MultTable{pt}.d_imu(wk))
                    d = obj.MultTable{pt}.d_imu(wk).gyroY; 
                    ridx = any(imidx<1,2)|any(imidx>size(d,1),2)|any(isnan(mt.imidx),2); %remove index
                    imidx(ridx,:) = [];
                    d = reshape(d(imidx,:),size(imidx,1),size(imidx,2));
                    d_gy(idx(~ridx),:) = d; %trial x time
                end

            end
            d_np = permute(d_np,[2,1]);
            d_xs = permute(d_xs,[2,1]);
            d_gz = permute(d_gz,[2,1]);
            d_am = permute(d_am,[2,1]);
            d_kd = permute(d_kd,[2,1]);
            d_wv = permute(d_wv,[2,3,1]); %time x freq x trial
            d_ed = permute(d_ed,[2,1]); %time x trial
            d_pe = permute(d_pe,[2,3,1]); %time x freq x trial
            d_gy = permute(d_gy,[2,1]); %time x trial

            obj.MultTrans.MTd = [];

            obj.MultTrans.MTd.d_np = d_np;
            obj.MultTrans.MTd.d_xs = d_xs;
            obj.MultTrans.MTd.d_gz = d_gz;
            obj.MultTrans.MTd.d_am = d_am;
            obj.MultTrans.MTd.d_kd = d_kd;
            obj.MultTrans.MTd.d_wv = d_wv;
            obj.MultTrans.MTd.d_ed = d_ed;
            obj.MultTrans.MTd.d_gy = d_gy;
            obj.MultTrans.MTd.d_pe1 = squeeze(d_pe(:,1,:)); %theta
            obj.MultTrans.MTd.d_pe2 = squeeze(d_pe(:,2,:)); %alpha
            obj.MultTrans.MTd.d_pe3 = squeeze(d_pe(:,3,:)); %gamma

            obj.MultTrans.MTd.tsec_np = tsec_np;
            obj.MultTrans.MTd.tsec_xs = tsec_xs;
            obj.MultTrans.MTd.tsec_gz = tsec_gz;
            obj.MultTrans.MTd.tsec_am = tsec_am;
            obj.MultTrans.MTd.tsec_kd = tsec_kd;
            obj.MultTrans.MTd.tsec_wv = tsec_wv;
            obj.MultTrans.MTd.tsec_ed = tsec_bi;
            obj.MultTrans.MTd.tsec_gy = tsec_im;
            obj.MultTrans.MTd.tsec_pe1 = tsec_pe;
            obj.MultTrans.MTd.tsec_pe2 = tsec_pe;
            obj.MultTrans.MTd.tsec_pe3 = tsec_pe;

            obj.MultTrans.MTd.tsamp_np = tsamp_np;
            obj.MultTrans.MTd.tsamp_xs = tsamp_xs;
            obj.MultTrans.MTd.tsamp_gz = tsamp_gz;
            obj.MultTrans.MTd.tsamp_am = tsamp_am;
            obj.MultTrans.MTd.tsamp_kd = tsamp_kd;
            obj.MultTrans.MTd.tsamp_wv = tsamp_wv;
            obj.MultTrans.MTd.tsamp_ed = tsamp_bi;
            obj.MultTrans.MTd.tsamp_gy = tsamp_im;
            obj.MultTrans.MTd.tsamp_pe1 = tsamp_pe;
            obj.MultTrans.MTd.tsamp_pe2 = tsamp_pe;
            obj.MultTrans.MTd.tsamp_pe3 = tsamp_pe;

            obj.MultTrans.MTd.freq_wv = freq_wv;
            obj.MultTrans.MTd.freq_pe = freq_pe;
        end

        function parseSegData(obj,varargin)
            %Parse data into 2sec segments. Run loadData
            %first. This parses based on in/out boundaries. Make sense to
            %chunk up data sequentially from the start and tag segments as
            %in or out. Need to try this at some point...
            %parseSegData;
            DS_np = []; %2 sec segments marked as inside/outside
            DS_xs = []; %velocity data for xsens in 2sec segments to match DS_np
            DS_gz = [];
            DS_am = [];
            DS_kd = [];
            DS_wv = [];
            DS_bi = [];
            DS_pe = [];
            DS_im = [];
            InFlag = logical([]); %flag for inside segments
            WalkNumSeg = []; %walk number for segments

            TT = {}; %in/out start/end times table

            for m=1:size(obj.DTable,1) %by walk
                walknum = obj.WalkNums(m);

                dbb = obj.DB(obj.DB.Walk==walknum,:);

                %NP
                d_np = obj.DTable.d_np{m};
                ntp_np = obj.DTable.ntp_np{m}; %ntp in sec
                fs_np = obj.DTable.fs_np(m);

                %XS
                d_xs = obj.DTable.d_xs{m};
                ntp_xs = obj.DTable.ntp_xs{m}; %ntp in sec
                fs_xs = obj.DTable.fs_xs(m);

                %Gaze
                d_gz = obj.DTable.d_gaze_fix{m};
                ntp_gz = obj.DTable.ntp_gaze{m}; %ntp in sec
                fs_gz = obj.DTable.fs_gaze(m);

                %Ambient
                d_am = obj.DTable.d_amb{m};
                ntp_am = obj.DTable.ntp_amb{m}; %ntp in sec
                fs_am = obj.DTable.fs_amb(m);

                %KDE
                d_kd = obj.DTable.d_kde{m};
                ntp_kd = obj.DTable.ntp_kde{m}; %ntp in sec
                fs_kd = obj.DTable.fs_kde(m);

                %Wavelet
                d_wv = obj.DTable.d_wav{m};
                ntp_wv = obj.DTable.ntp_wav{m}; %ntp in sec
                fs_wv = obj.DTable.fs_wav(m);
                f_wv = obj.DTable.f_wav{m};

                %Biopac
                if isempty(obj.DTable.d_bio{m})
                    d_bi = [];
                else
                    d_bi = table2array(obj.DTable.d_bio{m});
                end
                ntp_bi = obj.DTable.ntp_bio{m}; %ntp in sec
                fs_bi = obj.DTable.fs_bio(m);
                f_bi = obj.DTable.d_bio{find(~cellfun(@isempty,obj.DTable.d_bio),1)}.Properties.VariableNames;

                %Bosc pepisode
                d_pe = obj.DTable.d_pep{m};
                ntp_pe = obj.DTable.ntp_pep{m}; %ntp in sec
                fs_pe = obj.DTable.fs_pep(m);
                f_pe = obj.DTable.f_pep{m};

                %IMU
                d_im = struct2array(obj.DTable.d_imu(m));
                ntp_im = obj.DTable.ntp_imu{m}; %ntp in sec
                fs_im = obj.DTable.fs_imu(m);
                f_im = fieldnames(obj.DTable.d_imu(m));

                %Segmenting data into 2sec segments and classifying as
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

                % Chunk into 2sec (500sample) segments for inside/outside
                win_seg_sec = obj.WinSegSec; %default 2sec
                win_seg_np = win_seg_sec*fs_np;
                win_seg_xs = win_seg_sec*fs_xs;
                win_seg_gz = win_seg_sec*fs_gz;
                win_seg_am = win_seg_sec*fs_am;
                win_seg_kd = win_seg_sec*fs_kd;
                win_seg_wv = win_seg_sec*fs_wv;
                win_seg_bi = win_seg_sec*fs_bi;
                win_seg_pe = win_seg_sec*fs_pe;
                win_seg_im = win_seg_sec*fs_im;
                for k=1:size(T,1) %by indoor/outoor segment
                    if (T.Start_NTP(k)-ntp_np(1))<0
                        error('NP start is after start of walk! This is a problem and needs to be checked.');
                    end

                    if ~isempty(d_xs)
                        if (T.Start_NTP(k)-ntp_xs(1))<0
                            error('XS start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_gz)
                        if (T.Start_NTP(k)-ntp_gz(1))<0
                            error('Gaze start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_am)
                        if (T.Start_NTP(k)-ntp_am(1))<0
                            error('Amb start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_kd)
                        if (T.Start_NTP(k)-ntp_kd(1))<0
                            error('KDE start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_wv)
                        if (T.Start_NTP(k)-ntp_wv(1))<0
                            error('Wavelet start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_bi)
                        if (T.Start_NTP(k)-ntp_bi(1))<0
                            error('BIO start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_pe)
                        if (T.Start_NTP(k)-ntp_pe(1))<0
                            error('Bosc pepisode start is after start of walk! This is a problem and needs to be checked.');
                        end
                    end

                    if ~isempty(d_im)
                        if (T.Start_NTP(k)-ntp_im(1))<0
                            error('IMU start is after start of walk! This is a problem and needs to be checked.');
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
                    d_gz_seg = nan(win_seg_gz,size(d_np_seg,2));
                    if ~isempty(d_gz)
                        [~,start_idx] = min(abs(ntp_gz - T.Start_NTP(k))); %start samples for gaze
                        [~,end_idx] = min(abs(ntp_gz - T.End_NTP(k)));
                        dd_gz = d_gz(start_idx:end_idx);
                        dd_gz(floor(size(dd_gz,1)/win_seg_gz)*win_seg_gz+1:end,:) = []; %200Hz (4*200=800 -> 4sec)
                        dd_gz = reshape(dd_gz,win_seg_gz,[]);
                        if size(dd_gz,2)>size(d_np_seg,2)
                            d_gz_seg = dd_gz(:,1:size(d_np_seg,2));
                            fprintf('Gaze data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_gz,2)<size(d_np_seg,2)
                            d_gz_seg(:,1:size(dd_gz,2)) = dd_gz;
                            fprintf('Gaze data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_gz_seg = dd_gz;
                        end
                    end
                    DS_gz = cat(2,DS_gz,d_gz_seg);

                    %ambient light data (make the same size as d_np_seg)
                    d_am_seg = nan(win_seg_am,size(d_np_seg,2));
                    if ~isempty(d_am)
                        [~,start_idx] = min(abs(ntp_am - T.Start_NTP(k))); %start samples for amb
                        [~,end_idx] = min(abs(ntp_am - T.End_NTP(k)));
                        dd_am = d_am(start_idx:end_idx);
                        dd_am(floor(size(dd_am,1)/win_seg_am)*win_seg_am+1:end,:) = []; %15Hz (4*15=60 -> 4sec)
                        dd_am = reshape(dd_am,win_seg_am,[]);
                        if size(dd_am,2)>size(d_np_seg,2)
                            d_am_seg = dd_am(:,1:size(d_np_seg,2));
                            fprintf('Ambient data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_am,2)<size(d_np_seg,2)
                            d_am_seg(:,1:size(dd_am,2)) = dd_am;
                            fprintf('Ambient data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_am_seg = dd_am;
                        end
                    end
                    DS_am = cat(2,DS_am,d_am_seg);

                    %kde data (make the same size as d_np_seg)
                    d_kd_seg = nan(win_seg_kd,size(d_np_seg,2));
                    if ~isempty(d_kd)
                        [~,start_idx] = min(abs(ntp_kd - T.Start_NTP(k))); %start samples for kde
                        [~,end_idx] = min(abs(ntp_kd - T.End_NTP(k)));
                        dd_kd = d_kd(start_idx:end_idx);
                        dd_kd(floor(size(dd_kd,1)/win_seg_kd)*win_seg_kd+1:end,:) = []; %60Hz (4*60=240 -> 4sec)
                        dd_kd = reshape(dd_kd,win_seg_kd,[]);
                        if size(dd_kd,2)>size(d_np_seg,2)
                            d_kd_seg = dd_kd(:,1:size(d_np_seg,2));
                            fprintf('KDE data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_kd,2)<size(d_np_seg,2)
                            d_kd_seg(:,1:size(dd_kd,2)) = dd_kd;
                            fprintf('KDE data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_kd_seg = dd_kd;
                        end
                    end
                    DS_kd = cat(2,DS_kd,d_kd_seg);

                    %wavelet data (make the same size as d_np_seg)
                    d_wv_seg = nan(win_seg_wv,size(d_np_seg,2),length(f_wv),4);
                    if ~isempty(d_wv)
                        [~,start_idx] = min(abs(ntp_wv - T.Start_NTP(k))); %start samples for wav
                        [~,end_idx] = min(abs(ntp_wv - T.End_NTP(k)));
                        dd_wv = d_wv(start_idx:end_idx,:,:); %time x freq x chan
                        dd_wv(floor(size(dd_wv,1)/win_seg_wv)*win_seg_wv+1:end,:,:) = []; %60Hz (4*60=240 -> 4sec)
                        dd_wv = reshape(dd_wv,win_seg_wv,[],length(f_wv),4);
                        if size(dd_wv,2)>size(d_np_seg,2)
                            d_wv_seg = dd_wv(:,1:size(d_np_seg,2),:,:);
                            fprintf('WAV data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_wv,2)<size(d_np_seg,2)
                            d_wv_seg(:,1:size(dd_wv,2),:,:) = dd_wv;
                            fprintf('WAV data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_wv_seg = dd_wv;
                        end
                    end
                    DS_wv = cat(2,DS_wv,d_wv_seg);

                    %Biopac data (make the same size as d_np_seg)
                    d_bi_seg = nan(win_seg_bi,size(d_np_seg,2),length(f_bi));
                    if ~isempty(d_bi)
                        [~,start_idx] = min(abs(ntp_bi - T.Start_NTP(k))); %start samples for bio
                        [~,end_idx] = min(abs(ntp_bi - T.End_NTP(k)));
                        dd_bi = d_bi(start_idx:end_idx,:); %time x datatype
                        dd_bi(floor(size(dd_bi,1)/win_seg_bi)*win_seg_bi+1:end,:) = []; %trimming to nearest 2sec chunck
                        dd_bi = reshape(dd_bi,win_seg_bi,[],length(f_bi));
                        if size(dd_bi,2)>size(d_np_seg,2)
                            d_bi_seg = dd_bi(:,1:size(d_np_seg,2),:);
                            fprintf('BIO data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_bi,2)<size(d_np_seg,2)
                            d_bi_seg(:,1:size(dd_bi,2),:) = dd_bi;
                            fprintf('BIO data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_bi_seg = dd_bi;
                        end
                    end
                    DS_bi = cat(2,DS_bi,d_bi_seg);

                    %bosc pepisode data (make the same size as d_np_seg)
                    d_pe_seg = nan(win_seg_pe,size(d_np_seg,2),size(f_pe,1),4);
                    if ~isempty(d_pe)
                        [~,start_idx] = min(abs(ntp_pe - T.Start_NTP(k))); %start samples for pe
                        [~,end_idx] = min(abs(ntp_pe - T.End_NTP(k)));
                        dd_pe = d_pe(start_idx:end_idx,:,:); %time x freq x chan
                        dd_pe(floor(size(dd_pe,1)/win_seg_pe)*win_seg_pe+1:end,:,:) = []; 
                        dd_pe = reshape(dd_pe,win_seg_pe,[],size(f_pe,1),4); %2sec time x # segments x freq x chan
                        if size(dd_pe,2)>size(d_np_seg,2)
                            d_pe_seg = dd_pe(:,1:size(d_np_seg,2),:,:);
                            fprintf('PEP data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_pe,2)<size(d_np_seg,2)
                            d_pe_seg(:,1:size(dd_pe,2),:,:) = dd_pe;
                            fprintf('PEP data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_pe_seg = dd_pe;
                        end
                    end
                    DS_pe = cat(2,DS_pe,d_pe_seg);

                    %IMU data (make the same size as d_np_seg)
                    d_im_seg = nan(win_seg_im,size(d_np_seg,2),length(f_im));
                    if ~isempty(d_im)
                        [~,start_idx] = min(abs(ntp_im - T.Start_NTP(k))); %start samples for bio
                        [~,end_idx] = min(abs(ntp_im - T.End_NTP(k)));
                        dd_im = d_im(start_idx:end_idx,:); %time x datatype
                        dd_im(floor(size(dd_im,1)/win_seg_im)*win_seg_im+1:end,:) = []; %trimming to nearest 2sec chunck
                        dd_im = reshape(dd_im,win_seg_im,[],length(f_im));
                        if size(dd_im,2)>size(d_np_seg,2)
                            d_im_seg = dd_im(:,1:size(d_np_seg,2),:);
                            fprintf('IMU data for walk %0.0f (in/out segment %0.0f) is longer than NP. Truncating to fit NP.\n',walknum,k);
                        elseif size(dd_im,2)<size(d_np_seg,2)
                            d_im_seg(:,1:size(dd_im,2),:) = dd_im;
                            fprintf('IMU data for walk %0.0f (in/out segment %0.0f) is shorter than NP. Missing time is filled with NaNs.\n',walknum,k);
                        else
                            d_im_seg = dd_im;
                        end
                    end
                    DS_im = cat(2,DS_im,d_im_seg);

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

            obj.DParsed.Seg.DS_np = DS_np;
            obj.DParsed.Seg.DS_xs = DS_xs;
            obj.DParsed.Seg.DS_gz = DS_gz;
            obj.DParsed.Seg.DS_am = DS_am;
            obj.DParsed.Seg.DS_kd = DS_kd;
            obj.DParsed.Seg.DS_wv = DS_wv;
            obj.DParsed.Seg.DS_bi = DS_bi;
            obj.DParsed.Seg.DS_pe = DS_pe;
            obj.DParsed.Seg.DS_im = DS_im;
            obj.DParsed.Seg.FS_np = fs_np;
            obj.DParsed.Seg.FS_xs = fs_xs;
            obj.DParsed.Seg.FS_gz = fs_gz;
            obj.DParsed.Seg.FS_am = fs_am;
            obj.DParsed.Seg.FS_kd = fs_kd;
            obj.DParsed.Seg.FS_wv = fs_wv;
            obj.DParsed.Seg.FS_bi = fs_bi;
            obj.DParsed.Seg.FS_pe = fs_pe;
            obj.DParsed.Seg.FS_im = fs_im;
            obj.DParsed.Seg.F_wv = f_wv; %freq vector for wavelet
            obj.DParsed.Seg.F_bi = f_bi; %field names for biopac
            obj.DParsed.Seg.F_pe = f_pe; %freq bands for pepisode
            obj.DParsed.Seg.F_im = f_im; %field names for imu
            obj.DParsed.Seg.InFlag = InFlag;
            obj.DParsed.Seg.WalkNumSeg = WalkNumSeg;
            obj.DParsed.Seg.Outliers = any(any(isnan(DS_np),1),3); %<-500 outliers handled in loadData (set to nan)
            obj.DParsed.Trans.WinTransSec = obj.WinTransSec;
            obj.DParsed.Seg.WinSegSec = obj.WinSegSec;
        end %parseSegData

        function filterMultSegData(obj,varargin)
            %Filters seg data based on:
            %
            %regiontype -> 'AntHipp', 'LatTemp', 'Ent+Peri', 'PostHipp+Para', 'All Chans'
            %patienttype -> 'All Patients', '1,2,5'
            %walktype -> 'First Walks', 'Last Walks', 'Stop Walks', 'Go Walks', 'All Walks', '2,3'
            %
            %filterMultSegData('regiontype','PostHipp+Para','walktype','All Walks','patienttype','All Patients');
            p = obj.parseInputs(varargin{:});

            if isempty(p.regiontype)
                error('regiontype must be specified!');
            end
            if isempty(p.walktype)
                error('walktype must be specified!');
            end
            if isempty(p.patienttype)
                error('patienttype must be specified!');
            end

            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            customregion = p.customregion; %nx2 matrix where cols are patient,chan
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', [individual walks]
            patienttype = p.patienttype; %'All Patients', '2,3'

            if isempty(obj.MultSeg)
                obj.getMultSegData;
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
                    if ~all(ismember(regioncell,unique(obj.MultSeg.SegTable.RegionLabel)))
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

            MS = obj.MultSeg.SegTable;
          
            %filter by region (amygdala/anterior hipp)
            if all(regionidx<5)
                idx = ~ismember(MS.RegionNum,regionidx); %1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
                MS(idx,:) = [];
            elseif any(regionidx==6) %custom regions (nx2 matrix where cols are patient,chan)
                idx = false(size(MS,1),1);
                for k=1:size(customregion,1)
                    idx = idx | all(MS.Pt==customregion(k,1) & MS.Chan==customregion(k,2),2);
                end
                MS(~idx,:) = [];
            end

            %filter by walk
            if walkidx<5
                if walkidx==-1 %walktype is a number
                    idx = ~ismember(MS.Walk,str2num(walktype)); %filter out all but specified walk number
                elseif walkidx<3
                    idx = (MS.Walk>3); %keep 1st 3 walks
                    if walkidx==2 %keep last walks
                        idx = ~idx;
                    end
                elseif walkidx<4
                    idx = (MS.StopWalk~=1); %stop walks
                else
                    idx = (MS.GoWalk~=1); %go walks
                end
                MS(idx,:) = [];
            end

            %filter by patient
            if patientidx==-1
                idx = ~ismember(MS.Pt,str2num(patienttype)); %filter out all but specified walk number
                MS(idx,:) = [];
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
                idx = ~ismember(MS.Pt,ptlist); %filter out all but specified walk number
                MS(idx,:) = [];
            end
        
            obj.MultSeg.MS = MS;
          
            disp('Selected regions:')
            disp(unique(MS.RegionLabel));
            disp('Selected patients:')
            disp(unique(MS.Pt)'); %pt
            disp('Selected walks:')
            disp(unique(MS.Walk)'); %walk
            disp('Total trials:')
            disp(size(MS,1));

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

        function saveMultData(obj,varargin)
            disp('Saving...')
            multseg = obj.MultSeg;
            multtrans = obj.MultTrans;
            multtable = obj.MultTable;
            save(obj.AnalysisFile,'multseg','multtrans','multtable','-v7.3');
            disp('Finished saving.')
        end

        function loadMultData(obj,varargin)
            [~,name] = fileparts(obj.AnalysisFile);
            rootname = [name,'.mat'];
            if isfile(obj.AnalysisFile) || isfile(rootname)
                if isfile(rootname)
                    analysisfile = rootname;
                else
                    analysisfile = obj.AnalysisFile;
                end
                disp('Loading data for multiple patients from file...')
                MD = load(analysisfile,'multseg','multtrans','multtable');
                if ~isempty(MD.multseg)
                    obj.MultSeg = MD.multseg;
                end
                if ~isempty(MD.multtrans)
                    obj.MultTrans = MD.multtrans;
                end
                if ~isempty(MD.multtable)
                    obj.MultTable = MD.multtable;
                end
            else
                disp('Analysis file not found!')
            end
        end

        function calcfBOSCWalk(obj,d,fs,varargin)
            %calcfBOSCWalk(d,250); %walk1
            %This runs fooof in python. If libiomp5md.dll error, add
            %KMP_DUPLICATE_LIB_OK=1 as a system variable in windows.

            obj.FB.D = d; %data for full walk from  1 patient (4 chans)
            obj.FB.fs = fs;

            obj.FB.D(isnan(obj.FB.D)) = 0;

            obj.FB.nt = size(obj.FB.D,1); %time
            obj.FB.nc = size(obj.FB.D,2); %channels
            obj.FB.t = (0:obj.FB.nt-1)/obj.FB.fs;
            obj.FB.f_bands = [4,8;8,12;12,30]; %Theta,Alpha,Beta
            % obj.FB.f_bands = [4,8;8,12;12,30;30,70;70,120]; %this takes about 10x longer!!
            obj.FB.f_band_lbls = {'theta','alpha','beta','gamma','hg'};
           
            obj.FB.data.label = {'chan1','chan2','chan3','chan4'};
            obj.FB.data.time = {obj.FB.t};
            obj.FB.data.trial = {obj.FB.D'};
            obj.FB.data.fsample = obj.FB.fs;

            % general fBOSC setup
            obj.FB.cfg.fBOSC.F                 = min(obj.FB.f_bands(:)):0.5:max(obj.FB.f_bands(:));
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
            obj.FB.cfg.fBOSC.threshold.duration	= repmat(10, 1, numel(obj.FB.cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc, typically 3 cycles, 10 seems excessive?)
            obj.FB.cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold

            % episode post-processing
            obj.FB.cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now

            % general processing settings
            obj.FB.cfg.fBOSC.channel           = []; % select posterior channels (default: all)
            obj.FB.cfg.fBOSC.trial             = []; % select trials (default: all)
            obj.FB.cfg.fBOSC.trial_background  = []; % select trials for background (default: all)

            [obj.FB.fBOSC, obj.FB.cfg] = fBOSC_wrapper(obj.FB.cfg, obj.FB.data);

            %freq bands of interest
            % ds_f_spec = obj.DParsed.Seg.F_wv;
            % ds_f_band = [4,8;8,12;12,30;30,70;70,120]; %Theta,Alpha,Beta,Gamma,HG
            % ds_nband = size(ds_f_band,1);
            % ds_idx60 = (ds_f_spec>=59 & ds_f_spec<=61);
            % ds_fidx = false(length(ds_f_spec),5);
            % for n=1:ds_nband
            %     ds_fidx(:,n) = (ds_f_spec>=ds_f_band(n,1) & ds_f_spec<ds_f_band(n,2)) & ~ds_idx60;
            % end

            % collapsing detected oscillations to bands of interest
            t = obj.FB.t;
            etable = obj.FB.fBOSC.episodes;
            f_band = obj.FB.f_bands;
            nband = size(f_band,1);
            nchan = obj.FB.nc;
            ntime = length(t);
            min_dur = 10/fs*3; %minimum duration in sec, PEP is downsampled later by x10, so this will ensure that detections will contain at least 3 samples after the downsample process
            PEP = zeros(ntime,nband,nchan); %detected pepisodes
            for n=1:nband
                fidx = (etable.FrequencyMean>=f_band(n,1) & etable.FrequencyMean<f_band(n,2));
                ett = etable(fidx,:);
                for m=1:nchan
                    cidx = (ett.Channel==m);
                    et = ett(cidx,:); et(et.DurationS<min_dur,:) = [];
                    dt = zeros(ntime,1);
                    for r=1:size(et,1)
                        dt(t>=et.Onset(r) & t<et.Offset(r)) = 1;
                    end
                    PEP(:,n,m) = dt;
                end
            end
            obj.FB.PEP = PEP;

        end %calcfBOSCWalk

        function [PercOverlap,AvgOverlapSec] = calcOverlap(~,MT,transrng)
            %Calculates the amount of overlap between transition events
            uPWX = unique(MT(:,{'patient','walk','npidx'})); %need to remove chans that have the same dataidx for a given patient/walk
            [uPW,~,uPWIdx] = unique(uPWX(:,{'patient','walk'}));
            DiffSec = []; %difference between events in same patient/walk converted to sec
            for k=1:size(uPW,1)
                idx = (k==uPWIdx);
                diffsec = diff(sort(uPWX.npidx(idx)))./250; %difference between events in same patient/walk converted to sec
                DiffSec = cat(1,DiffSec,diffsec);
            end
            OverlapWinSec = diff(transrng); %size of specgram window in sec
            OverlapIdx = DiffSec<OverlapWinSec;
            PercOverlap = sum(OverlapIdx)./numel(DiffSec)*100; %?? Need to think this through carefully
            AvgOverlapSec = mean(OverlapWinSec-DiffSec(OverlapIdx));
        end

    end %methods

    %%%%%%%%%%%%%%%%%%%%%
    methods %general plotting

        function openMultTransSpecGramGUI(obj,varargin)
            MultTransSpecGramGUI(obj);
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
            plotpowercomp = p.plotpowercomp;
            plotfooofcomp = p.plotfooofcomp;
            boxcomprng = p.boxcomprng; %2x4 matrix where rows are each box and cols are low/high ranges for time/freq
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

            %Getting trials and specgram data
            obj.filterMultTransData(varargin{:}); %table of all trials and corresponding info (this creates MT and must be run before getFilteredMultTransData)
            obj.getFilteredMultTransData(transrng); %raw power for all trials (time x freq x trial)

            %Init some params
            tsamp = obj.MultTrans.MTd.tsamp_wv; %+-10sec at 25Hz (downsampled by 10 from 250)
            tsec = obj.MultTrans.MTd.tsec_wv;
            baseidx = tsec>=normrng(1) & tsec<=normrng(2); %baseline (normalization) indices
            nperm = 1000; %permutations

            pwr = obj.MultTrans.MTd.d_wv;
            tsec = obj.MultTrans.MTd.tsec_wv;
            freq = obj.MultTrans.MTd.freq_wv;
            np = obj.MultTrans.MTd.d_np;
            tsec_np = obj.MultTrans.MTd.tsec_np;

            predlist = [{'xs','gz','kd','am','ed','pe1','pe2','pe3','gy'};{'Vel','Fix','KDE','Amb','Eda','PeT','PeA','PeB','HeadTurn'}];

            nfreq = length(freq);
            ntime = length(tsamp);
            ntrials = size(obj.MultTrans.MT,1);

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
            pH = nan(size(predlist,2),3);
            yyaxis(aH,'right');
            for k=1:size(predlist,2) %xs,gz,kd,am,ed,pe1,pe2,pe3,gy
                dx = obj.MultTrans.MTd.(['d_',predlist{1,k}]);
                tx = obj.MultTrans.MTd.(['tsec_',predlist{1,k}]);
                if contains(predlist{1,k},'pe') || contains(predlist{1,k},'gy')
                    dxx = dx;
                else
                    dxx = dx(:,~isoutlier(mean(dx)));
                end
                dxx_me = mean(dxx,2,"omitnan");
                dxx_se = 1.96*std(dxx,0,2,'omitnan')./sqrt(size(dxx,2)); %95ci
                pH(k,1) = plot(aH,tx,dxx_me,'-k','LineWidth',2);
                pH(k,2) = plot(aH,tx,dxx_me-dxx_se,':k','LineWidth',0.5);
                pH(k,3) = plot(aH,tx,dxx_me+dxx_se,':k','LineWidth',0.5);
            end  
            predidx = strcmp(predlist(2,:),p.predtype);
            set(pH(~predidx,:),'visible','off');
            ylabel(aH,p.predtype);
            switch p.predtype %vel:0-1.5, gaze:0-1, kde:1-7, amb: 0-500, eda:-0.2-0.5, pe:0-1
                case 'Vel'
                    ylim(aH,[0,1.5])
                case 'Fix'
                    ylim(aH,[0,1])
                case 'KDE'
                    ylim(aH,[1,7])
                case 'Amb'
                    ylim(aH,[0,500])
                case 'Eda'
                    ylim(aH,[-0.2,0.5])
                case {'PeT','PeA','PeB'}
                    ylim(aH,[0,0.3])
                case 'HeadTurn' 
                    ylim(aH,[-40,40])
                otherwise
                    aH.YAxis(2).Visible = 'off';
            end
            yyaxis(aH,'left');
            if plottrials
                plot(aH,[tsec(1),tsec(end)],repmat(trialsfreqrng,2,1),':k')
            end
            if plotpowercomp || plotfooofcomp
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
            ptstr = regexprep(num2str(unique([obj.MultTrans.MT.patient])'),'\s+','');
            if isempty(str2num(walktype))
                wktypestr = walktype;
            else
                wktypestr = 'Sel Walks';
            end
            wkstr = regexprep(num2str(unique([obj.MultTrans.MT.walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr = cell2mat(rgcell(unique(obj.MultTrans.MT.regionnum)));
            if p.fullwalknorm
                normstr = 'full';
            else
                normstr = regexprep(num2str(normrng),'\s+','t');
            end
            transstr = regexprep(num2str(transrng),'\s+','t');
            [PercOverlap,AvgOverlapSec] = obj.calcOverlap(obj.MultTrans.MT,transrng);
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

                [mt,sidx] = sortrows(obj.MultTrans.MT,{'patient','chan','walk'});
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

            [uPtCh,~,uPtChIdx] = unique(obj.MultTrans.MT(:,{'patient','regionnum'}));
            mrk_style = {'o','+','*','x'}; %by region -> 1=AntHipp, 2=LatTemp, 3=Ent+Peri, 4=PostHipp+Para
            mrk_size = [6,10,10,10];
            mrk_color = {'r','g','b','c','m'}; %by patient -> 1:5
            ptnum = uPtCh.patient;
            rgnum = uPtCh.regionnum;

            fH3 = [];
            if plotpowercomp
                tidx_box1 = tsec>boxcomprng(1,1) & tsec<boxcomprng(1,2);
                tidx_box2 = tsec>boxcomprng(2,1) & tsec<boxcomprng(2,2);
                fidx_box1 = freq>boxcomprng(1,3) & freq<boxcomprng(1,4);
                fidx_box2 = freq>boxcomprng(2,3) & freq<boxcomprng(2,4);
                
                % pwr_box1 = squeeze(mean(mean(10*log10(pwr(tidx_box1,fidx_box1,:)),1),2));
                % pwr_box2 = squeeze(mean(mean(10*log10(pwr(tidx_box2,fidx_box2,:)),1),2));

                pwr_box1 = squeeze(mean(mean(pwr(tidx_box1,fidx_box1,:),1),2)); %average before converting to dB
                pwr_box2 = squeeze(mean(mean(pwr(tidx_box2,fidx_box2,:),1),2));

                % tidx_xs_box1 = tsec_xs>boxcomprng(1,1) & tsec_xs<boxcomprng(1,2);
                % tidx_xs_box2 = tsec_xs>boxcomprng(2,1) & tsec_xs<boxcomprng(2,2);

                % dxx = dx;
                % dxx(:,isoutlier(mean(dx))) = nan;
                % dx_box1 = mean(dxx(tidx_xs_box1,:))';
                % dx_box2 = mean(dxx(tidx_xs_box2,:))';

                mpwr_box1 = nan(size(uPtCh,1),1);
                mpwr_box2 = nan(size(uPtCh,1),1);
                % ccdx_box1 = nan(size(uPtCh,1),1);
                % ccdx_box2 = nan(size(uPtCh,1),1);
                for k=1:size(uPtCh,1)
                    idx = (uPtChIdx==k);
                    pwr_b1 = pwr_box1(idx);
                    pwr_b2 = pwr_box2(idx);
                    % dx_b1 = dx_box1(idx);
                    % dx_b2 = dx_box2(idx);
                    mpwr_box1(k) = mean(pwr_b1);
                    mpwr_box2(k) = mean(pwr_b2);
                    % nan_idx = isnan(dx_b1)|isnan(dx_b2);
                    % cc_b1 = corrcoef(dx_b1(~nan_idx),10*log10(pwr_b1(~nan_idx))); %convert pwr to dB before correlation
                    % cc_b2 = corrcoef(dx_b2(~nan_idx),10*log10(pwr_b2(~nan_idx)));
                    % ccdx_box1(k) = cc_b1(1,2);
                    % ccdx_box2(k) = cc_b2(1,2);
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
            end %power comparison

            fH4 = [];
            if plotfooofcomp
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
            end %fooof comparison

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
           
            %Getting trials and specgram data
            obj.filterMultTransData('transtype',transtype{1},'regiontype',regiontype{1},'walktype',walktype{1},'desctype',desctype{1},'patienttype',patienttype{1},'veltype',veltype{1},'customregion',customregion{1}); %table of all trials and corresponding info
            obj.getFilteredMultTransData(transrng); %raw power for all trials (time x freq x trial)

            %Init some params
            MTd1 = obj.MultTrans.MTd;
            MT1 = obj.MultTrans.MT;
            pwr1 = MTd1.d_wv;
            tsamp = MTd1.tsamp_wv; %+-10sec at 25Hz (downsampled by 10 from 250)
            tsec = MTd1.tsec_wv;
            freq = MTd1.freq_wv;
            ntrials1 = size(MT1,1);
            nfreq = length(freq);
            ntime = length(tsamp);

            obj.filterMultTransData('transtype',transtype{2},'regiontype',regiontype{2},'walktype',walktype{2},'desctype',desctype{2},'patienttype',patienttype{2},'veltype',veltype{2},'customregion',customregion{2}); %table of all trials and corresponding info
            obj.getFilteredMultTransData(transrng); %raw power for all trials (time x freq x trial)

            MTd2 = obj.MultTrans.MTd;
            MT2 = obj.MultTrans.MT;
            pwr2 = MTd2.d_wv;
            ntrials2 = size(MT2,1);
            % baseidx = tsec>=normrng(1) & tsec<=normrng(2); %baseline (normalization) indices
            nperm = 1000; %permutations

            predlist = [{'xs','gz','kd','am','ed'};{'Vel','Fix','KDE','Amb','Eda'}];

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
            pH = nan(size(predlist,2),3);
            yyaxis(aH,'right');
            for k=1:size(predlist,2) %xs,gz,kd,am,ed
                dx1 = MTd1.(['d_',predlist{1,k}]);
                dx2 = MTd2.(['d_',predlist{1,k}]);
                tx = MTd1.(['tsec_',predlist{1,k}]);
                dxx_me = mean(dx1(:,~isoutlier(mean(dx1))),2,"omitnan")-mean(dx2(:,~isoutlier(mean(dx2))),2,"omitnan");
                dxx_se = nan(size(dxx_me));
                pH(k,1) = plot(aH,tx,dxx_me,'-k','LineWidth',2);
                pH(k,2) = plot(aH,tx,dxx_me-dxx_se,':k','LineWidth',0.5);
                pH(k,3) = plot(aH,tx,dxx_me+dxx_se,':k','LineWidth',0.5);
            end  
            predidx = strcmp(predlist(2,:),p.predtype);
            set(pH(~predidx,:),'visible','off');
            ylabel(aH,p.predtype);
            switch p.predtype
                case 'Vel'
                    ylim(aH,[-2,2])
                case 'Fix'
                    ylim(aH,[-1,1])
                case 'KDE'
                    ylim(aH,[-2,2])
                case 'Amb'
                    ylim(aH,[-250,250])
                case 'Eda'
                    ylim(aH,[-1,1])
                otherwise
                    aH.YAxis(2).Visible = 'off';
            end
            yyaxis(aH,'left');
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
            %Correlation coefficient of features
            p = parseInputs(obj,varargin{:});

            if isempty(p.regiontype)
                error('regiontype must be specified!');
            end
            if isempty(p.walktype)
                error('walktype must be specified!');
            end
            if isempty(p.patienttype)
                error('patienttype must be specified!');
            end
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,5'
            patienttype = p.patienttype; %'All Patients', '2,4'

            %Getting trials and specgram data
            obj.filterMultSegData(varargin{:});

            ms = obj.MultSeg.MS;
            uptwk = obj.MultSeg.uPtWalk;
            lbls = {'wvTheta','wvAlpha','wvBeta','wvGamma','wvHG','Vel','Fix','KDE','Amb','EDA'};
            % lbls = {'mtTheta','mtAlpha','mtBeta','mtGamma','mtHG','Vel','Fix','KDE','Amb','EDA'};
            for k=1:size(uptwk,1)
                idx = (ms.uPtWalk==k);
                ms(idx,lbls) = ms(idx,lbls)-mean(ms(idx,lbls),'omitnan'); %subtract mean by ptwk
            end
            ms_array = table2array(ms(:,lbls));
            idx = any(isnan(ms_array),2);
            ms(idx,:) = [];
            ms_array(idx,:) = [];
           
            [C,P] = corrcoef(ms_array);
            C(P>=0.001) = 0;
            
            cmap = parula;
            % cmap(129,:) = [0,0,0];

            fH = figure('Colormap',cmap);
            aH = axes('Parent',fH);
            imagesc(C,'Parent',aH);
            clim(aH,[-0.5,0.5])
            set(aH,'xtick',1:length(lbls),'xticklabel',lbls,'ytick',1:length(lbls),'yticklabel',lbls)

            if isempty(str2num(patienttype))
                pttypestr = patienttype;
            else
                pttypestr = 'Sel Patients';
            end
            ptstr = regexprep(num2str(unique([ms.Pt])'),'\s+','');
            if isempty(str2num(walktype))
                wktypestr = walktype;
            else
                wktypestr = 'Sel Walks';
            end
            wkstr = regexprep(num2str(unique([ms.Walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr = cell2mat(rgcell(unique(ms.RegionNum)));
            ntrials = size(ms,1);
            
            ttlstr = sprintf('%s-%s, %s-%s, %s-%s (n=%s)',...
                regiontype,rgstr,pttypestr,ptstr,...
                wktypestr,wkstr,num2str(ntrials));
            title(aH,ttlstr); 
            cb = colorbar(aH);
            ylabel(cb,'correlation')
           
            if p.saveflag
                fstr = [regexprep(ttlstr,'\s+',''),'.png'];
                print(fH,fullfile('D:\RWN\ConfMatrix',fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function fH = plotMultInOutByRegion(obj,varargin)
            %Generates a boxplot of power comparing in/out segments across
            %all patients and channels.
            %
            %plotMultInOutByRegion('predtype','wvTheta','saveflag',true); %Theta, Alpha, Beta, Gamma, HG (wv or mt prefix for wavelet or multitaper)
            p = obj.parseInputs(varargin{:});
            if isempty(obj.MultSeg)
                obj.getMultData;
            end
            if isempty(p.predtype)
                error('predtype is missing!')
            end

            SegTable = obj.MultSeg.SegTable; %TD: added Fix and Amb 20240418 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)

            if ~any(contains(SegTable.Properties.VariableNames,p.predtype))
                error('predtype was incorrect!');
            end

            flist = obj.MultSeg.f_band;
            flbls = obj.MultSeg.f_lbls;
            fidx = find(contains(flbls,regexprep(p.predtype,'^(wv|mt)','')),1);
            frng = flist(fidx,:);

            %Creating table that is mean across uPtChan for glme
            SS = []; %side-by-side comparison of in/out
            SSmdl = []; %in/out is inline for glme
            for k=1:size(obj.MultSeg.uPtChan,1)
                idx = (obj.MultSeg.uPtChanIdx==k);

                pt = obj.MultSeg.uPtChan.Pt(k);
                chan = obj.MultSeg.uPtChan.Chan(k);
                pred = SegTable.(p.predtype)(idx); 
                inflag = SegTable.InFlag(idx);

                if p.rmoutliers
                    idx = isoutlier(pred);
                    pred(idx) = [];
                    inflag(idx) = [];
                end

                pred_in = pred(inflag);
                pred_out = pred(~inflag);

                in_n = sum(~isnan(pred_in));
                out_n = sum(~isnan(pred_out));

                if contains(p.predtype,{'Fix','peTheta','peAlpha','peBeta'})
                    mpred_in = mean(pred_in,'omitnan');
                    mpred_out = mean(pred_out,'omitnan');
                else
                    mpred_in = median(pred_in,'omitnan');
                    mpred_out = median(pred_out,'omitnan');
                end

                SS = vertcat(SS,table(...
                        k,pt,chan,...
                        in_n,out_n,...
                        mpred_in,mpred_out,...
                        'VariableNames',{...
                        'uPtChan','Pt','Chan',...
                        'In_N','Out_N',...
                        ['In_',p.predtype],['Out_',p.predtype],...
                        }));

                SSmdl = vertcat(SSmdl,table(...
                        [k;k],[pt;pt],[chan;chan],...
                        [true;false],[in_n;out_n],...
                        [mpred_in;mpred_out],...
                        'VariableNames',{...
                        'uPtChan','Pt','Chan',...
                        'InFlag','InFlag_N',...
                        p.predtype,...
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
            if contains(p.predtype,flbls,"IgnoreCase",true)
                ss = [SS.(['In_',p.predtype]),SS.(['Out_',p.predtype])].*10; %log10 already applied to median normalized (across walk) power
            else
                ss = [SS.(['In_',p.predtype]),SS.(['Out_',p.predtype])];
            end
            ss = ss-mean(ss,2); %since log10, need to subtract instead of divide by mean (i.e. negative numbers mess up the division -> if mean is negative it can actually flip the direction of the response!!) 
            [h,pval,ci,stats] = ttest2(ss(:,1),ss(:,2));

            InOutTable = removevars(SS,{['In_',p.predtype],['Out_',p.predtype]});
            InOutTable = addvars(InOutTable,ss(:,1),ss(:,2),'NewVariableNames',{['In_',p.predtype],['Out_',p.predtype]});

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
            if contains(p.predtype,flbls,"IgnoreCase",true)
                fstr2 = [regexprep(num2str(frng),'\s+','-'),' Hz'];
                ylabel(aH,['Normalized (dB) Power (',fstr2,')'])
            else
                ylabel(aH,['Normalized ',p.predtype])
            end
            set(aH,'FontSize',14);
            legend([pH_pt,pH_rg],{'Patient1','Patient2','Patient3','Patient4','Patient5','AntHipp','LatTemp','Ent+Peri','PostHipp'},'location','northwest')

%             [p,tbl,stats] = kruskalwallis(ss,labels,'on'); %use raw responses here since these are stats by subject
%             figure;
%             [c,m] = multcompare(stats,'alpha',0.05,'display','on');
            
            title(aH,sprintf('Inside vs Outside By Region (%s, p=%0.1e)',p.predtype,pval))
            % ylim(aH,ylimit)

            if p.saveflag
                fstr = sprintf('InOut_Comparison_%d-%dHz_wAllLabels.png',frng(1),frng(2));
                fdir = fullfile('D:\RWN\InOutFreq2');
                if ~isfolder(fdir)
                    mkdir(fdir);
                end
                print(fH,fullfile(fdir,fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function fH = plotMultInOutByWalk(obj,varargin)
            %Generates a boxplot of power comparing in/out segments across
            %all patients and channels for specified predictor.
            %
            %plotMultInOutPred('predtype','Vel','saveflag',true); 
            p = obj.parseInputs(varargin{:});
            if isempty(obj.MultSeg)
                obj.getMultData;
            end
            if isempty(p.predtype)
                error('predtype is missing!')
            end

            SegTable = obj.MultSeg.SegTable; %TD: added Fix and Amb 20240418 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)

            if ~any(contains(SegTable.Properties.VariableNames,p.predtype))
                error('predtype was incorrect!');
            end

            flist = obj.MultSeg.f_band;
            flbls = obj.MultSeg.f_lbls;
            fidx = find(contains(flbls,regexprep(p.predtype,'^(wv|mt)','')),1);
            frng = flist(fidx,:);

            %Creating table that is mean across uPtWalk
            SS = []; %side-by-side comparison of in/out
            SSmdl = []; %in/out is inline for glme
            for k=1:size(obj.MultSeg.uPtWalk,1)
                idx = (obj.MultSeg.uPtWalkIdx==k);

                sstrials = SegTable(idx,:);
                % sstrials = sstrials(sstrials.Chan==1,:);

                pt = obj.MultSeg.uPtWalk.Pt(k);
                walk = obj.MultSeg.uPtWalk.Walk(k);

                pred = sstrials.(p.predtype); 
                inflag = sstrials.InFlag;

                if p.rmoutliers
                    idx = isoutlier(pred);
                    pred(idx) = [];
                    inflag(idx) = [];
                end

                pred_in = pred(inflag);
                pred_out = pred(~inflag);

                in_n = sum(~isnan(pred_in));
                out_n = sum(~isnan(pred_out));

                if contains(p.predtype,{'Fix','peTheta','peAlpha','peBeta'})
                    in_pred = mean(pred_in,'omitnan');
                    out_pred = mean(pred_out,'omitnan');
                else
                    in_pred = median(pred_in,'omitnan');
                    out_pred = median(pred_out,'omitnan');
                end

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
            if contains(p.predtype,flbls,"IgnoreCase",true)
                ss = [SS.In_Pred,SS.Out_Pred].*10; %convert to dB for power
            else
                ss = [SS.In_Pred,SS.Out_Pred];
            end
            ss = ss-mean(ss,2);
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
            rg_style = {'o','x','*','+','square','diamond','^','>'}; %walk 1,2,3,4,5,6,7,8
            rg_size = [7,10,10,10,10,10,10,10];
            for m=1:size(ss,1)
                rg_idx = SS.Walk(m);
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
            if contains(p.predtype,flbls,"IgnoreCase",true)
                fstr2 = [regexprep(num2str(frng),'\s+','-'),' Hz'];
                ylabel(aH,['Normalized (dB) Power (',fstr2,')'])
            else
                ylabel(aH,['Normalized ',p.predtype])
            end
            set(aH,'FontSize',14);
            legend([pH_pt,pH_rg],{'Patient1','Patient2','Patient3','Patient4','Patient5','Walk1','Walk2','Walk3','Walk4','Walk5','Walk6','Walk7','Walk8'},'location','northwest')

%             [p,tbl,stats] = kruskalwallis(ss,labels,'on'); %use raw responses here since these are stats by subject
%             figure;
%             [c,m] = multcompare(stats,'alpha',0.05,'display','on');
            
            title(aH,sprintf('Inside vs Outside By Walk (%s, p=%0.1e)',p.predtype,pval))

            if p.saveflag
                fstr = sprintf('InOut_Comparison_%s_wAllLabels.png',p.predtype);
                fdir = fullfile('D:\RWN\InOutPred2');
                if ~isfolder(fdir)
                    mkdir(fdir);
                end
                print(fH,fullfile(fdir,fstr),'-dpng','-r300');
                close(fH);
            end
        end

        function openMultSegGLMEGUI(obj,varargin)
            MultSegGLMEGUI(obj);
        end

        function varargout = plotMultGLME(obj,varargin)
            %Runs glme for specified model
            %
            %plotMultGLME('regiontype','All Chans','walktype','All Walks','patienttype','All Patients','glmemdl','nTheta ~ InFlag + (1|uPtChan)');
            p = obj.parseInputs(varargin{:});

            if isempty(p.glmemdl)
                error('glmemdl must be specified!');
            end
            regiontype = p.regiontype; %'AntHipp','LatTemp','Ent+Peri','PostHipp+Para','All Chans','Custom'
            walktype = p.walktype; %'First Walks','Last Walks','Stop Walks','Go Walks','All Walks', '1,2,5'
            patienttype = p.patienttype; %'All Patients', '2,4'
            glmemdl = p.glmemdl; %nTheta ~ InFlag + (1|uPtChan)

            obj.filterMultSegData(varargin{:});

            ms = obj.MultSeg.MS;
            ms.HeadTurn = abs(ms.HeadTurn); %left turn is negative, right is positive, need to take abs or these will cancel out
            % uptch = obj.MultSeg.uPtChan;

            glme = fitglme(ms,glmemdl);

            [~,fENames] = fixedEffects(glme);
            fENames = fENames.Name;
            fENames(strcmp(fENames,'(Intercept)')) = [];
            fENames = regexprep(fENames,'_1','');

            % [~,rENames] = randomEffects(glme);

            % nperm = 1000;
            % npred = length(fENames)+1; %+1 to include intercept
            % nobs = size(ms,1); %includes nans (ntrials below excludes nans)
            % respvar = glme.ResponseName; 
            % T = glme.Coefficients.tStat';
            % T_lbls = regexprep(glme.Coefficients.Name','_1','');
            %
            % T_perm = nan(nperm,npred);
            % parfor k=1:nperm
            %     ms_rand = ms;
            %     for m=1:size(uptch,1)
            %         idx = find(ms_rand.uPtChan==m);
            %         ms_rand.(respvar)(idx) = ms.(respvar)(idx(randperm(length(idx))));
            %     end
            %     glme_perm = fitglme(ms_rand,glmemdl);
            %     T_perm_lbls = regexprep(glme_perm.Coefficients.Name','_1','');
            %     if all(ismember(T_lbls,T_perm_lbls)) && (length(T_lbls)==length(T_perm_lbls))
            %         T_perm(k,:) = glme_perm.Coefficients.tStat';
            %     end
            % end
            % 
            % if any(isnan(T_perm(:)))
            %     disp('Labels did not match in permutation testing!');
            % end
            % 
            % pvals = sum(T>=T_perm,1)./nperm;
            % pvals = min(pvals,1-pvals);
            % pvals(pvals<=0) = 1/nperm;
            % pvals = 2*pvals; %correct for 2-tail test
            % P = pvals; %pvalues for the predictors based on permutation testing
            
            % LL = glme.LogLikelihood;
            % AC = glme.ModelCriterion.AIC;
            % TS = {[glme.Coefficients.Name,num2cell(glme.Coefficients.tStat),num2cell(glme.Coefficients.pValue)]};

            %SSR = sum of squares for regression
            %SSE = sum of squares for error (also sum of squared residuals)
            %SST = sum of squares total
            %SST = SSR+SSE;
            %Ordinary R2 = 1 - (SSE/SST);
            %Adjusted R2 = 1 - (SSE/(n-p-1))/(SST/(n-1)); where n is number of observations and p is number of predictors -> penalizes more predictors
            %Partial R2 = (SSE_reduced - SSE_full)/SSE_reduced; SSE_reduced is caculcated by leaving out the predictor of interest
            %
            % The partial R^2 values do not sum up to equal the full R^2.
            % Here's why:
            %
            % ### Understanding Partial R^2
            % Partial R^2 measures the unique contribution of each
            % predictor to the variance explained by the model, holding all
            % other predictors constant. It highlights the specific
            % explanatory power of each individual predictor.
            % 
            % ### Understanding Full R^2
            % Full R^2 represents the total proportion of variance
            % explained by the entire model with all predictors included.
            % 
            % ### Why They Don’t Sum Up
            % Partial R^2 values, when summed, often exceed the full R^2
            % because each partial R^2 reflects the incremental explanatory
            % power of each predictor over and above all others. The total
            % R^2 captures the combined effect of all predictors, which
            % includes overlaps and interactions between predictors.
            % 
            % Imagine the model’s explanatory power as a Venn diagram. Each
            % predictor's contribution overlaps with the others, making
            % their unique contributions (partial R^2 values) important,
            % but not directly additive to the overall model’s R^2.

            pnames = regexp(glmemdl,'\([1|a-zA-Z\+]+\)','match'); %get random effects
            mdlstr = regexprep(glmemdl,'\([1|a-zA-Z\+]+\)',''); %remove random effects
            mdlstr = strtrim(regexp(mdlstr,'~','split')); mdlstr(1) = [];
            mdlstr = strtrim(regexp(cell2mat(mdlstr),'\+','split')); 
            mdlstr(cellfun(@isempty,mdlstr)) = [];
            pnames = [mdlstr,pnames];
            if contains(pnames(1),'-1')
                pnames(1) = [];
            end
            R2_partial = nan(length(pnames),1);
            if length(pnames)>1
                for k=1:length(pnames)
                    mdl = pnames{k};
                    mdl = regexprep(mdl,'\(','\\(');
                    mdl = regexprep(mdl,'\)','\\)');
                    mdl = regexprep(mdl,'\|','\\|');
                    mdl = regexprep(mdl,'\+','\\+');
                    mdl = regexprep(glmemdl,[mdl,'\s*\+*\s*'],'');
                    mdl = regexprep(mdl,'(\s*\+*\s*)$','');
                    mdl = regexprep(mdl,'~\s+(','~ 1 + (');
                    mdl = regexprep(mdl,'\+\|','|');
                    disp(mdl);
                    glm = fitglme(ms,mdl);
                    R2_partial(k) = (glm.SSE-glme.SSE)/glm.SSE;
                end
            end

            if isempty(str2num(patienttype))
                pttypestr = patienttype;
            else
                pttypestr = 'Sel Patients';
            end
            ptstr = regexprep(num2str(unique([ms.Pt])'),'\s+','');
            if isempty(str2num(walktype))
                wktypestr = walktype;
            else
                wktypestr = 'Sel Walks';
            end
            wkstr = regexprep(num2str(unique([ms.Walk])'),'\s+','');
            rgcell = {'a','l','e','p'}; %AntHipp, LatTemp, Ent+Peri, PostHipp+Para
            rgstr = cell2mat(rgcell(unique(ms.RegionNum)));
            ntrials = glme.NumObservations;
            
            ttlstr = sprintf('%s-%s, %s-%s, %s-%s (n=%s)',...
                regiontype,rgstr,pttypestr,ptstr,...
                wktypestr,wkstr,num2str(ntrials));

            out_str = ttlstr;
            out_str = sprintf('%s\n\n%s',out_str,regexprep(evalc('disp(glme)'),'\%|<strong>|</strong>|',''));
            out_str = sprintf('%s\nR2 = %0.4f\n',out_str,glme.Rsquared.Ordinary);
            for k=1:length(pnames)
                out_str = sprintf('%sR2_%s = %0.4f\n',out_str,pnames{k},R2_partial(k));
            end
            % out_str = sprintf('%s\nPermutation Results:\n',out_str);
            % for k=1:length(T_lbls)
            %     out_str = sprintf('%s%s p = %0.4f\n',out_str,T_lbls{k},P(k));
            % end
         
            fH = figure('position',[50,50,1800,1000],'Visible','on');
            fH.Name = [glmemdl,', ',ttlstr];
            tl = tiledlayout(2,10,'Parent',fH,'TileSpacing','compact','Padding','compact');
            aH = nexttile(tl,1,[2,3]);
            text(aH,0,1,out_str,'VerticalAlignment','top','Interpreter','none'); 
            axis(aH,'off');

            for k=1:length(fENames)
                Adjusted = calcAdjRespGLMEFit(ms,glme,fENames{k});
                fobj = Adjusted.PolyFit;
                x = Adjusted.Predictor;
                yadj = Adjusted.Response;
                idx = contains(coeffnames(fobj),'p1');
                val = coeffvalues(fobj);
                % v95 = confint(fobj);
                r2 = Adjusted.R2;

                if k>3
                    aH2 = nexttile(tl,(k-1)*2+9,[1,2]);
                else
                    aH2 = nexttile(tl,(k-1)*2+5,[1,2]);
                end
               
                plot(aH2,fobj,x,yadj);
                t = title(aH2,sprintf('est=%0.5f, r2=%0.5f',val(idx),r2));
                p = glme.Coefficients.pValue(contains(glme.Coefficients.Name,fENames{k}));
                if p<0.05
                    t.FontWeight = 'bold';
                else
                    t.FontWeight = 'normal';
                end
                legend(aH2,'off');
                xlabel(aH2,fENames{k})
                ylabel(aH2,glme.ResponseName)
            end

            if nargout
                varargout{1} = fH;
            end
           
        end

        function plotMultGLMEAllCombos(obj,varargin)
            %Runs glme models for all predictor combos
            if isempty(obj.MultSeg)
                obj.getMultSegData;
            end

            p = obj.parseInputs(varargin{:});

            SegTable = obj.MultSeg.SegTable; %TD: added Fix and Amb 20240418 and KDE 20240513 (contains raw pwr for all segments - not mean - and includes xs/fix/amb which may contain nans)
            % f = obj.MultSeg.f;

            %Model list
            pred_list = {'InFlag','nVel','nFix','nAmb','nKDE','nEDA'};
            mdl_template = {};
            for k=1:length(pred_list)
                idx = nchoosek(1:length(pred_list),k);
                for n=1:size(idx,1)
                    mdl = ['nPwr ~ ',cell2mat(cellfun(@(x)[x,' + '],pred_list(idx(n,:)),'UniformOutput',false)),'(1|uPtChan)'];
                    mdl_template = cat(1,mdl_template,mdl);
                end
            end
            % load('GLME_sidx.mat','sidx');
            % sidx = 1:63;

            % flist = [[2,5];[5,8];[8,12];[12,30];[30,60];[60,85];[85,120]];
            flist = obj.MultSeg.f_band;
            flbls = obj.MultSeg.f_lbls;
            LL = nan(length(mdl_template),size(flist,1));
            AC = nan(length(mdl_template),size(flist,1));
            TS = cell(length(mdl_template),size(flist,1));
            for m=1:size(flist,1)
                frng = flist(m,:);
                % fidx = f>frng(1)&f<=frng(2);

                %Mean across freq (5-8Hz) in log10 space
                % resp_var = ['nPwr',regexprep(num2str(frng),'\s+','')];
                resp_var = ['n',flbls{m}];
                % SegTable = addvars(SegTable,mean(SegTable.nPwr(:,fidx),2),'NewVariableNames',resp_var,'Before','Pwr'); %nPwr is log10

                %Model list
                mdl_list = regexprep(mdl_template,'nPwr',resp_var);

                for k=1:length(mdl_list)
                    disp(mdl_list{k});
                    glme = fitglme(SegTable,mdl_list{k});
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

                if m==1
                    [~,sidx] = sort(AC(:,m),'ascend');
                end
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

        function plotInOutSpecDiff(obj,varargin)
            %Plots spectrum comparison of inside vs outside
            %plotInOutDiff('walknum',1); only walk 1
            %plotInOutDiff('fpass',[2,85]); different fpass
            %plotInOutDiff; all walks
            p = obj.parseInputs(varargin{:});

            DS = obj.DParsed.Seg.DS_np;
            InFlag = obj.DParsed.Seg.InFlag;
            WalkNumSeg = obj.DParsed.Seg.WalkNumSeg;

            % DS_xs = obj.DParsed.Seg.DS_xs;
            % DS_xs(DS_xs>5) = nan;
            % DS_xs = mean(DS_xs,'omitnan');
            % 
            % nan_idx = isnan(DS_xs);
            % DS_xs(nan_idx) = [];
            % DS(:,nan_idx,:) = [];
            % InFlag(nan_idx) = [];
            % WalkNumSeg(nan_idx) = [];
            % 
            % InFlag = DS_xs>1; %velocity greater than 1

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

            % InFlag = InFlag(randperm(length(InFlag)));

            params.Fs = 250;
            params.trialave = 0;
            params.err = [2,0.05];
            params.tapers = [3,5];
            params.fpass = [fpass(1)-1,fpass(2)+1];

            fH = figure('Position',[50,50,1700,500],'Visible','on');
            for k=1:4 %by channel
                disp(k);
                ds = DS(:,:,k);
                inflag = InFlag;
                nanidx = any(isnan(ds));
                ds(:,nanidx) = [];
                inflag(nanidx) = [];

                params.trialave = 0;
                [S,f] = mtspectrumc(ds,params);
                % [PVal,PMsk] = kstestPerm_mex(single(permute(S,[2,3,1])),InFlag(:)',200,0.05);
                [PVal,PMsk] = kstestPerm(single(permute(S,[2,3,1])),inflag(:)',200,0.05);

                SE = strel('rectangle',[3 1]);
                PMsk = imdilate(PMsk,SE);
                PMsk = imerode(PMsk,SE);
                B = regionprops(PMsk(:)','BoundingBox');
                B = cat(1,B.BoundingBox);
                % disp(B(:,3))
                % if ~isempty(B)
                %     B(B(:,3)<3,:) = [];
                % end

                params.trialave = 1;
                [S_in,~,~,Serr_in] = mtspectrumc(ds(:,inflag),params);
                [S_out,f,~,Serr_out] = mtspectrumc(ds(:,~inflag),params);

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

            % fstr = sprintf('%s_AllValidWalks_Spectrum_%0.0f-%0.0fHz_Inside-Outside_n%0.0f-n%0.0f.png',obj.PatientID,fpass(1),fpass(2),sum(InFlag),sum(~InFlag));
            % print(fH,fullfile('C:\Users\Administrator\Dropbox\Work\Code\Analysis\Inman\RealWorld\figs',obj.PatientID,fstr),'-dpng','-r300');
            % close(fH);
        end

        function plotMultInOutSpecDiff(obj,varargin)
            %Generates a power difference plot comparing in/out segments
            %across all patients and channels.
            %
            %plotMultInOutDiff;
            SpecTable = [];
            for k=1:length(obj.PatientList)
                disp(k);
                if ~any(obj.PatientIdx==k)
                    obj.loadData('PatientIdx',k); %this now generates downsampled/normalized wavelet data
                end

                obj.parseSegData;

                DS_np = obj.DParsed.Seg.DS_np;
                DS_wv = obj.DParsed.Seg.DS_wv;
                DS_wk = obj.DParsed.Seg.WalkNumSeg;
                DS_in = obj.DParsed.Seg.InFlag;

                DS_ol = obj.DParsed.Seg.Outliers; %these have nans for <-500 condition
                DS_np(:,DS_ol,:) = [];
                DS_wv(:,DS_ol,:,:) = [];
                DS_wk(DS_ol) = [];
                DS_in(DS_ol) = [];

                %multitaper params
                ds_params.Fs = obj.DTable.fs_np(1); %250
                ds_params.trialave = 0;
                ds_params.tapers = [5,9];
                ds_params.fpass = [2,120];
                ds_params.pad = 0;
                [~,ds_f_spec_mt] = mtspectrumc(DS_np(:,:,1),ds_params);

                %wavelet freq
                ds_f_spec_wv = obj.DParsed.Seg.F_wv;

                ds_nseg = size(DS_np,2); %number of segments after outlier removal
                for m=1:4 %by channel
                    ds_S = mtspectrumc(DS_np(:,:,m),ds_params); ds_S = ds_S'; %nseg x freq
                    ds_W = squeeze(median(DS_wv(:,:,:,m),1)); %median across time
                    ds_pt = repmat(k,ds_nseg,1);
                    ds_ch = repmat(m,ds_nseg,1);
                    spectable = table(...
                        ds_pt,ds_ch,DS_wk,DS_in,ds_S,ds_W,...
                        'VariableNames',...
                        {'Pt','Chan','Walk','InFlag','mtSpec','wvSpec',...
                        });
                    SpecTable = vertcat(SpecTable,spectable);
                end
            end

            SS = [];
            % f = ds_f_spec_mt;
            f = fliplr(ds_f_spec_wv');
            [uPtChan,~,uPtChanIdx] = unique(SpecTable(:,{'Pt','Chan'}),'rows');
            for k=1:size(uPtChan,1)
                idx = (uPtChanIdx==k);
                spectable = SpecTable(idx,:);
                % spec_in = mean(spectable.mtSpec(spectable.InFlag,:));
                % spec_out = mean(spectable.mtSpec(~spectable.InFlag,:));
                spec_in = fliplr(mean(spectable.wvSpec(spectable.InFlag,:)));
                spec_out = fliplr(mean(spectable.wvSpec(~spectable.InFlag,:)));
                ss = table(k,spectable.Pt(1),spectable.Chan(1),...
                    spec_in,spec_out,...
                    'VariableNames',...
                    {'uPtChan','Pt','Chan','InSpec','OutSpec'});
                SS = vertcat(SS,ss);
            end

            %Significance for average unique pt/chan response (similar to
            %boxplot ttest)
            P = []; T = []; MD = []; SE = [];
            for n=1:length(f)
                ss = [SS.InSpec(:,n),SS.OutSpec(:,n)];
                ss = ss./mean(ss,2);
                [h,p,ci,stats] = ttest2(ss(:,1),ss(:,2));
                P(n) = p;
                T(n) = stats.tstat;
                MD(n,:) = mean(ss);
                SE(n,:) = std(ss)./sqrt(size(ss,1));
            end
            T(P.*length(f)>0.05) = 0;

            se = strel('rectangle',[3 1]);
            PMsk = imdilate(reshape(T~=0,[],1),se);
            PMsk = imerode(PMsk,se);
            B = regionprops(PMsk(:)','BoundingBox');
            I = regionprops(PMsk(:)','PixelIdxList');
            B = cat(1,B.BoundingBox);
            disp(B(:,3))
            if ~isempty(B)
                idx = B(:,3)<3;
                B(idx,:) = [];
                rm_idxs = I(idx,:);
                rm_idxs = cat(1,rm_idxs.PixelIdxList);
                T(rm_idxs) = 0;
            end

            figure
            plot(f,T)
            xlabel('Hz')
            ylabel('tstat')
            title('Inside/Outside Difference Across Freq (all pts, all chans, t-test, bonferonni corr)')

            mInPwr = smoothdata(MD(:,1)'); sInPwr = smoothdata(SE(:,1))'.*2;
            mOutPwr = smoothdata(MD(:,2)'); sOutPwr = smoothdata(SE(:,2))'.*2;

            % fH = figure;
            % aH = axes('Parent',fH);
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

            fH = figure;
            aH = axes('Parent',fH);
            mOutInDiff = mOutPwr - mInPwr;
            sOutInDiff = sqrt(sInPwr.^2 + sOutPwr.^2);
            hold on;
            ylimit = [-0.02,0.07];
            pH = nan(size(B,1),1);
            for m=1:size(B,1)
                b = B(m,[1,3]);
                x = floor([b(1),sum(b),sum(b),b(1)]);
                x(x<1) = 1; x(x>length(f)) = length(f);
                y = [ylimit(1),ylimit(1),ylimit(2),ylimit(2)];
                pH(m) = patch(f(x),y,'k','FaceAlpha',0.05,'EdgeColor','none','parent',aH);
            end
            plot(f,mOutInDiff,'b');
            plot(round(f([1,end])),[0,0],'k')
            patch([f,fliplr(f)],([mOutInDiff+sOutInDiff,fliplr(mOutInDiff-sOutInDiff)]),'b','facealpha',0.3,'edgecolor','none');
            xlabel('Hz')
            ylabel('Normalized Power Difference (Outside-Inside)')
            title('Spectrum Difference (All Pts, All Chans, All Walks, 95% CI)')
            xlim(round(f([1,end])))
            ylim(ylimit)
            
        end %plotMultInOutDiff

        %%%%%%%%%%%%%% Old/Unused %%%%%%%%%%%%%%%%%%

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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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





