%%
rootdir = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2";
fname_np = fullfile(rootdir,"NeuroPace_PHI\UCLA_NEA_4344798_ECoG_Catalog.csv");
fname_vid = fullfile(rootdir,"Original\Walk1\Pupil\5d738720_0.0-1556.477.mp4");
fname_rp = fullfile(rootdir,"Original\Walk1\Raspberry\RP_marks_2021-08-09_10-31-32_871905.txt");
root_np = fullfile(rootdir,"NeuroPace_PHI\UCLA_NEA_4344798 Data EXTERNAL #PHI");
fname_world = fullfile(rootdir,"Original\Walk1\Pupil\world_timestamps.csv");
fname_drift = fullfile(rootdir,"Original\Walk1\PupilPhoneData\timeDriftRecord-2021-08-09_10-18-49-102.csv");


%% Match green led flashes in video to raspberry pi ntp times
[ntp_vid,vid] = findPupilVidNTP(fname_vid,fname_drift,fname_world); %1st is led, 2nd is csv timestamps
[ntp_np,D] = findNPaceNTP(fname_np,fname_rp,root_np);
t = datetime(ntp_np./(60*60*24),'convertFrom','datenum'); %convert from seconds to datetime format



%% Spectrum inside vs outside (rw2 - walk1, need to add offset=12 frames to these timestamps to match ntp_np)
datestr = '09-Aug-2021 ';
tIn = {
    '10:32:36.985','10:33:51.902'
    '10:36:44.564','10:39:40.760'
    '10:40:18.969','10:40:27.865'
    '10:41:01.682','10:41:17.678'
    '10:41:55.867','10:42:27.848'
    '10:47:52.536','10:49:38.729'
    };
tOut = {
    '10:33:58.845','10:36:34.768'
    '10:39:45.784','10:40:09.785'
    '10:40:32.661','10:40:56.662'
    '10:41:21.687','10:41:50.863'
    '10:42:31.840','10:47:47.812'
    };
tIn_sec = cellfun(@(x)datenum(datetime([datestr,x]))*60*60*24,tIn);
tOut_sec = cellfun(@(x)datenum(datetime([datestr,x]))*60*60*24,tOut);

%chunk into 4sec (1000sample) segments
DDin = [];
for k=1:size(tIn,1)
    idx = ntp_np>tIn_sec(k,1) & ntp_np<tIn_sec(k,2);
    d = D(idx,:);
    d(floor(size(d,1)/1000)*1000+1:end,:) = [];
    dd = reshape(d,1000,[],4);
    DDin = cat(2,DDin,dd);
end
b1 = any(any(DDin<-500,1),3);
b2 = any(any(isnan(DDin),1),3);
DDin(:,b1|b2,:) = [];

DDout = [];
for k=1:size(tOut,1)
    idx = ntp_np>tOut_sec(k,1) & ntp_np<tOut_sec(k,2);
    d = D(idx,:);
    d(floor(size(d,1)/1000)*1000+1:end,:) = [];
    dd = reshape(d,1000,[],4);
    DDout = cat(2,DDout,dd);
end
b1 = any(any(DDout<-500,1),3); %remove segments with mark artifact
b2 = any(any(isnan(DDout),1),3); %remove segments with nans
DDout(:,b1|b2,:) = [];

params.Fs = 250;
params.trialave = 1;
params.err = [1,0.05];
params.tapers = [3,5];
params.fpass = [0,50];

fH = figure('Position',[50,50,1700,500]);
for k=1:4
    [S_in,~,~,Serr_in] = mtspectrumc(DDin(:,:,k),params);
    [S_out,f,~,Serr_out] = mtspectrumc(DDout(:,:,k),params);

    aH = subplot(1,4,k,'parent',fH);
    patch([f,fliplr(f)],log10([Serr_in(1,:),fliplr(Serr_in(2,:))]),'b','facealpha',0.3,'edgecolor','none','parent',aH)
    hold on
    patch([f,fliplr(f)],log10([Serr_out(1,:),fliplr(Serr_out(2,:))]),'r','facealpha',0.3,'edgecolor','none','parent',aH)
%     axis(aH,[0,20,0,2.2])
    xlim([0,15])
    xlabel('Hz')
    ylabel('log10(uV2/Hz)')
    title(['Chan',num2str(k)])
    legend(aH,{'In','Out'})
end


%% Transitions
datestr = '09-Aug-2021 ';
tInOut = {
    '10:33:56.273'
    '10:39:44.957'
    '10:40:30.030'
    '10:41:19.950'
    '10:42:30.440'
    };
tOutIn = {
    '10:36:41.332'
    '10:40:17.000'
    '10:40:59.552'
    '10:41:53.840'
    '10:47:49.806'
    };
tInOut_sec = cellfun(@(x)datenum(datetime([datestr,x]))*60*60*24,tInOut);
tOutIn_sec = cellfun(@(x)datenum(datetime([datestr,x]))*60*60*24,tOutIn);

%InOut - 10sec segments centered on transition
DDio = nan(5001,size(tInOut,1),4); %time x trial x chan
for k=1:size(tInOut,1)
    idx = ntp_np>tInOut_sec(k)-10 & ntp_np<tInOut_sec(k)+10;
    d = D(idx,:);
    if length(d)<5001
        d(end+1:5001,:) = d(end,:);
    else
        d = d(1:5001,:);
    end
    DDio(:,k,:) = d;
end
b1 = any(any(DDio<-500,1),3);
b2 = any(any(isnan(DDio),1),3);
DDio(:,b1|b2,:) = [];

%OutIn - 10sec segments centered on transition
DDoi = nan(5001,size(tOutIn,1),4); %time x trial x chan
for k=1:size(tOutIn,1)
    idx = ntp_np>tOutIn_sec(k)-10 & ntp_np<tOutIn_sec(k)+10;
    d = D(idx,:);
    if length(d)<5001
        d(end+1:5001,:) = d(end,:);
    else
        d = d(1:5001,:);
    end
    DDoi(:,k,:) = d;
end
b1 = any(any(DDoi<-500,1),3);
b2 = any(any(isnan(DDoi),1),3);
DDoi(:,b1|b2,:) = [];

%%
chan = 1;
fH = figure('Position',[50,50,1000,450]);

aH = subplot(1,2,1,'parent',fH);
PSG = PermSpecGram(DDio(:,:,chan),'Fs',250,'FPass',[2,64],'TimeRng',[-10,10],'NormRng',[-5,0]);
PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','ErrPerc',0.05,'Correction','Pixel','Smoothing',[250,1]);
PSG.plotSpecGram('MaskFlag',false,'Clim',[-10,10],'aH',aH)
c = colorbar; c.Label.String = 'dB';
title(aH,sprintf('In->Out (chan%0.0f)',chan))
hold(aH,'on')
plot([0,0],[1,6],'k');

aH = subplot(1,2,2,'parent',fH);
PSG = PermSpecGram(DDoi(:,:,chan),'Fs',250,'FPass',[2,64],'TimeRng',[-10,10],'NormRng',[-5,0]);
PSG.calcSpecGram('AnalysisType','Pwr','NormType','Mean','ErrPerc',0.05,'Correction','Pixel','Smoothing',[250,1]);
PSG.plotSpecGram('MaskFlag',false,'Clim',[-10,10],'aH',aH)
c = colorbar; c.Label.String = 'dB';
title(aH,sprintf('Out->In (chan%0.0f)',chan))
hold(aH,'on')
plot([0,0],[1,6],'k');






