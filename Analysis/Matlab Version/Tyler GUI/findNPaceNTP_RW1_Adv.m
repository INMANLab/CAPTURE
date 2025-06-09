function [ntp,D,Fs] = findNPaceNTP_RW1_Adv(Files,varargin)
% For RW1, uses data from photon_export, then runs it throught the typical
% neuropace pipeline. Also works for other participants.
%
% fname_np = "\\RolstonServer\D\Data\RealWorldNavigationCory\RW2\NeuroPace_PHI\UCLA_NEA_4344798_ECoG_Catalog.csv";
% root_np = "\\RolstonServer\D\Data\RealWorldNavigationCory\RW2\NeuroPace_PHI\UCLA_NEA_4344798 Data EXTERNAL #PHI";
% fname_rp = "\\RolstonServer\D\Data\RealWorldNavigationCory\RW2\Original\Walk1\Raspberry\RP_marks_2021-08-09_10-31-32_871905.txt";

A = []; S = [];
if ~isempty(Files.rp_marks_file)
    fid = fopen(Files.rp_marks_file);
    A = fread(fid,[1,Inf],'*char');
    fclose(fid);

    A = regexp(A,'\n','split');
    S = A(contains(A,'START')); S = regexp(S,'START: ','split','once'); S = cat(1,S{:}); S(:,1) = [];
    A = A(contains(A,'MARK')); A = regexp(A,'MARK: ','split','once'); A = cat(1,A{:}); A(:,1) = []; %ntp times for current walk
end

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
if contains(Files.rp_marks_file,'RW1')
    load(Files.photon_export_file,"photon_export");

    if photon_export.raspberry.ntp.START~=0
        error('RP START time in photon_export was expected to be zero!');
    end

    D = photon_export.rns.data;
    nchans = size(D,1);
    Fs = photon_export.rns.fs;
    dt = find(diff(round(photon_export.rns.time.*Fs))~=1); %data start/end times (samples) for stitched segments in D
    dt = [[1,dt+1];[dt,size(D,2)]];

    if size(dt,2)>8 %if more than 8 data segments
        error('Too many data segments! Need to check...');
    end

    DTable = cell(size(dt,2),4);
    for k=1:size(dt,2)
        dat = D(:,dt(1,k):dt(2,k))';

        ts = photon_export.rns.time(dt(1,k)); %sec relative to RP START
        ts = ts + datenum(S{1})*60*60*24; %ntp time for start sample of np data from photon_export corrected for RP START (sec)
        dl = size(dat,1)./Fs; %data length in sec

        DTable(k,1) = {dat};
        DTable(k,2) = {ts};
        DTable(k,3) = {dl};
        DTable(k,4) = {size(dat,1)};
    end

else
    if contains(Files.rp_marks_file,'RW4')
        opts = detectImportOptions(Files.np_csv_file);
        opts = setvaropts(opts,{'Timestamp','RawUTCTimestamp','RawLocalTimestamp'},'InputFormat','MM/dd/uuuu HH:mm:ss.SSS');
        nptable = readtable(Files.np_csv_file, opts);
    else
        nptable = readtable(Files.np_csv_file);
    end
    nptable = sortrows(nptable,"Timestamp");

    if isempty(A) %RW5 walks 6-8
        world_tbl = readtable(Files.pupil_ts_file); %in units nanoseconds
        drift_tbl = readtable(Files.drift_csv_file); %in units milliseconds
        ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec
        ntp_pupil = datetime(world_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
        ntp_pupil = datenum(ntp_pupil)*60*60*24; %convert from days to seconds (use this to find start/stop time of walk when rp is missing)
        ntp_pupil = ntp_pupil + ntp_offset;

        idx = find(nptable.Timestamp>=datetime(ntp_pupil(1)/(60*60*24),'ConvertFrom','datenum') & nptable.Timestamp<=datetime(ntp_pupil(end)/(60*60*24),'ConvertFrom','datenum'));
        idx = [idx(1)-1;idx;idx(end)+1];
        nptable = nptable(idx,:); %remove all but current walk

        walk_str = regexp(Files.walk_dir,'Walk\d{1}$','match','once');
        switch walk_str
            case 'Walk6'
                nptable.Timestamp(2) = nptable.Timestamp(2)-seconds(0.5); %manually adjust timestamp for each data file (this is probably easier than detecting led reliably)
            case 'Walk7'
                nptable.Timestamp(2) = nptable.Timestamp(2)-seconds(1.16);
            case 'Walk8'
                nptable.Timestamp(2) = nptable.Timestamp(2)-seconds(0.96);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %need to populate A with times from LED onset in video
%         % Find sync led in video (need to match green color to improve sensitivity - white lights also have high green component - tried this and sensitivity got worse)
%         vid = VideoReader(Files.pupil_vid_file); 
%         fr = vid.frameRate;
%         nof = vid.NumFrames;
%         II = zeros(vid.Height,vid.Width);
%         I = zeros(nof,1);
%         wait_msg = parfor_wait(nof);
%         for k=1:nof
%             wait_msg.Send;
%             img = read(vid,k);
%             I(k) = mean(reshape(img(1:25,425:475,2),[],1));
%             %     image(img(1:25,425:475,2));
%             %     drawnow;
%         end
%         wait_msg.Destroy;
% 
%         %find threshold crossings
%         thresh = 225; %RW2 - walk 1
%         dthr = diff([false;I>thresh]);
%         [ts,tsidx] = find(dthr==1);
%         tsend = find(dthr==-1);
%         if length(ts)>length(tsend)
%             while 1
%                 idx = find((tsend-ts(1:length(tsend)))<0,1);
%                 ts(idx) = [];
%                 tsidx(idx) = [];
%                 if (length(ts)==length(tsend))
%                     break;
%                 end
%                 if isempty(idx)
%                     ts(end) = [];
%                     tsidx(end) = [];
%                 end
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        idx = find(nptable.Timestamp>=A{1} & nptable.Timestamp<=A{end});
        idx = [idx(1)-1;idx;idx(end)+1];
        nptable = nptable(idx,:); %remove all but current walk
    end

    nchans = nptable.WaveformCount(1);
    Fs = nptable.SamplingRate(1);

    DTable = cell(size(nptable,1),4);
    for k=1:size(nptable,1)
        fid = fopen(fullfile(Files.np_data_dir, nptable.Filename{k}));
        dat = reshape(fread(fid,'int16'),nchans,[])'-512;
        fclose(fid);

        ts = datenum(nptable.Timestamp(k))*60*60*24; %starting timestamp in seconds
        dl = size(dat,1)./Fs; %data length in sec

        DTable(k,1) = {dat};
        DTable(k,2) = {ts};
        DTable(k,3) = {dl};
        DTable(k,4) = {size(dat,1)};
    end
end
DTable = cell2table(DTable,'VariableNames',{'data','ntp','sec','samp'});

[ntp,D,TS] = stitchData(DTable,Fs,nchans);
ntp_diff = diff(ntp)-(1/Fs);
if min(ntp_diff)<-0.0001 || max(ntp_diff)>0.0001 %if more than 1/10ms jitter, send alert
    mesg = 'NTP time vector has too much variance! Need to check...';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end

[ts,chan] = findMarks(D); %ts is in samples

if length(ts)<10
    mesg = 'Not enough marks were detected in the neuropace data! Using neuropace spreadsheet ntp times instead of RP times.';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
    plotResult(ntp,D(:,chan));
    return;
end

%Finding offset between data marks and RP marks in sec
if isempty(A)
    t1 = ntp(ts);
else
    t1 = datenum(A)*60*60*24;  %ntp times when mark sent to np data (sec) from RP text file
end
t2 = ntp(ts);  %ntp times for actual marks in sec using timestamps for data files (these might not be accurate)

mesg = 'Removing NP and RP marks that do not match within 8sec of overlap...';
if nargin>3
    uialert(varargin{1},mesg,'');
end
disp(mesg);

%Find overlap between RP/NP marks and remove the rest.
[t1,t2,n] = findOvlp(t1,t2);

if length(t1)<10 || (length(t1)~=length(t2))
    mesg = 'Too many marks were removed! Using neuropace spreadsheet ntp times instead of RP times.';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
    plotResult(ntp,D(:,chan));
    return;
else
    mesg = sprintf('%0.0f marks were removed...',n);
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end

%Classifying marks based on data segment (last segment might not have a mark)
t2_seg = nan(length(t2),1);
for k=1:size(TS,1)
    idx = t2>=TS(k,2) & t2<=TS(k,4);
    t2_seg(idx) = k;
end

dt = t1-t2;
ol = isoutlier(abs(dt));
if any(ol)
    mesg = 'Outliers found in time difference between RP and neuropace marks! Double check this data!';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end

shift = mean(dt);
plotResult(ntp,D(:,chan),t1,t2,t2_seg);
if abs(shift)>1
    mesg = sprintf('Mean NTP difference between RP marks and neuropace timestamps is large (%0.2f sec)! Using neuropace spreadsheet ntp times instead of RP times.\n',shift);
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
    return;
end

%Shifting data segments to match marks
mesg = sprintf('Mean NTP difference between RP marks and neuropace timestamps is %0.2f sec. Shifting neuropace spreadsheet ntp times (by data segment) to match RP times.\n',shift);
if nargin>3
    uialert(varargin{1},mesg,'');
end
disp(mesg);

DTableShift = DTable;
useg = unique(t2_seg); %these correspond to rows in DTable/TS
shifts = zeros(size(DTableShift,1),1);
for k=1:length(useg)
    idx = (t2_seg==useg(k));
    shifts(useg(k)) = median(dt(idx));
    DTableShift.ntp(useg(k)) = DTable.ntp(useg(k)) + shifts(useg(k));
end
disp('Shifts by segment:')
disp(shifts');

%Process the data again with the shifts applied
[ntp,D,TS] = stitchData(DTableShift,Fs,nchans);
ntp_diff = diff(ntp)-(1/Fs);
if min(ntp_diff)<-0.0001 || max(ntp_diff)>0.0001 %if more than 1/10ms jitter, send alert
    mesg = 'NTP time vector has too much variance! Need to check...';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end
[ts,chan2] = findMarks(D); %ts is in samples
t2 = ntp(ts);
[t1,t2] = findOvlp(t1,t2);

t2_seg = nan(length(t2),1);
for k=1:size(TS,1)
    idx = t2>=TS(k,2) & t2<=TS(k,4);
    t2_seg(idx) = k;
end
plotResult(ntp,D(:,chan),t1,t2,t2_seg);

if length(t1)<10 || (length(t1)~=length(t2)) || chan~=chan2
    mesg = 'Something went wrong while shifting the data segments! Need to check this manually...';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
else
    mesg = 'Data alignment was successful!';
    if nargin>3
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ntp,D,TS] = stitchData(DTable,Fs,nchans)
TotSamp = round(((DTable.ntp(end)+DTable.sec(end)) - DTable.ntp(1))*Fs); %total number of samples based on ntp times

D = nan(TotSamp,nchans);
TS = nan(size(DTable,1),4);
for k=1:size(DTable,1)
    cs = round((DTable.ntp(k)-DTable.ntp(1))*Fs); %current sample
    dl = DTable.samp(k); %data samples

    D((1:dl)+cs,:) = DTable.data{k};
    TS(k,1) = cs+1;
    TS(k,2) = DTable.ntp(k);
    TS(k,3) = cs+dl;
    TS(k,4) = DTable.ntp(k)+DTable.sec(k);
end

ntp = interp1(TS(:,1),TS(:,2),1:size(D,1),'linear','extrap')'; %ntp times in sec for all data samples


function [ts,chan] = findMarks(D)
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


function [t1,t2,n] = findOvlp(t1,t2)
%Find RP/NP marks that overlap and remove the rest. t1 is array of RP ntp
%marks. t2 is array of NP ntp marks. n is number of marks removed.
nt1 = round(t1*1000); %convert to ms
nt2 = round(t2*1000);

offset = min(nt1(1),nt2(1));
nt1 = nt1 - offset + 1; %in ms and normalized to 1st value (RP marks)
nt2 = nt2 - offset + 1; %data marks

N = max([nt1;nt2]);

T1 = false(1,N);
T1(nt1) = true;
T1 = imdilate(T1,true(1,2000));

T2 = false(1,N);
T2(nt2) = true;
T2 = imdilate(T2,true(1,2000));

ovlp = T1 & T2;
ovlp = imdilate(ovlp,true(1,2000));

t1 = t1(ovlp(nt1)); %only keeping times that overlap within 8sec border (2000 samples)
t2 = t2(ovlp(nt2));

n = sum(~ovlp(nt1))+sum(~ovlp(nt2));


function plotResult(ntp,d,varargin)
figure; 
plot(ntp,d);
hold on;
diff_ntp = diff(ntp);
fprintf('Time difference min/max (%0.3e, %0.3e)\n',min(diff_ntp),max(diff_ntp));
if nargin>2
    t1 = varargin{1};
    pH1 = plot([t1,t1],[-600,200],'r');
end
if nargin>3
    t2 = varargin{2};
    pH2 = plot([t2,t2],[-600,200],'g');
    dt = t1-t2;
    title(sprintf('RP/NP mark time difference (mean=%0.3f,min=%0.3f,max=%0.3f)',mean(dt),min(dt),max(dt)));
    legend([pH1(1),pH2(1)],{'RP','NP'})
end
if nargin>4
    t2_seg = varargin{3};
    ut = unique(t2_seg);
    figure;
    dt = t1-t2;
    hold on;
    for k=1:length(ut)
        idx = (t2_seg==ut(k));
        plot(t2(idx),dt(idx),'o-')
    end
    title('NP/RP time differences by data segment')
    ylabel('sec')
end



