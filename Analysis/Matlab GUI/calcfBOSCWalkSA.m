function fB = calcfBOSCWalkSA(d,fs)
%Calculates pepisode for entire walk. This is the standalone version for
%use in parfor.
%calcfBOSCWalk(d,250); %walk1
%This runs fooof in python. If libiomp5md.dll error, add
%KMP_DUPLICATE_LIB_OK=1 as a system variable in windows.

fB.D = d; %data for full walk from  1 patient (4 chans)
fB.fs = fs;

fB.D(isnan(fB.D)) = 0;

fB.nt = size(fB.D,1); %time
fB.nc = size(fB.D,2); %channels
fB.t = (0:fB.nt-1)/fB.fs;
fB.f_bands = [4,8;8,12;12,30]; %Theta,Alpha,Beta
fB.f_band_lbls = {'theta','alpha','beta'};

fB.data.label = {'chan1','chan2','chan3','chan4'};
fB.data.time = {fB.t};
fB.data.trial = {fB.D'};
fB.data.fsample = fB.fs;

% general fBOSC setup
fB.cfg.fBOSC.F                 = 4:0.5:30;
fB.cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
fB.cfg.fBOSC.fsample           = fB.fs;         % current sampling frequency of EEG data

% padding
fB.cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
fB.cfg.fBOSC.pad.detection_s   = 0.1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
fB.cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof parameters - fit with fixed line or allow a knee
fB.cfg.fBOSC.fooof.aperiodic_mode    = 'knee'; %old = eBOSC not fooof
fB.cfg.fBOSC.fooof.version           = 'python'; %'matlab'

% threshold settings
fB.cfg.fBOSC.threshold.duration	= repmat(10, 1, numel(fB.cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc, typically 3 cycles, 10 seems excessive?)
fB.cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold

% episode post-processing
fB.cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now

% general processing settings
fB.cfg.fBOSC.channel           = []; % select posterior channels (default: all)
fB.cfg.fBOSC.trial             = []; % select trials (default: all)
fB.cfg.fBOSC.trial_background  = []; % select trials for background (default: all)

[fB.fBOSC, fB.cfg] = fBOSC_wrapper(fB.cfg, fB.data);

% collapsing detected oscillations to bands of interest
t = fB.t;
etable = fB.fBOSC.episodes;
f_band = fB.f_bands;
nband = size(f_band,1);
nchan = fB.nc;
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
fB.PEP = PEP;

