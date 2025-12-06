%% Code to apply fBOSC to RW data. 
%This runs fooof in python. If libiomp5md.dll error, add
%KMP_DUPLICATE_LIB_OK=1 as a system variable in windows.
clear; clc;

A = load("\\rolstonserver\D\Data\RealWorldNavigationCory\RW2\Original\Walk1\RWNApp_RW2_Walk1.mat");

d = A.d_np';
d(d<-500) = 0;
d(isnan(d)) = 0;
% d = d(:,6.5e4:10e4);

fs = A.fs_np;
nc = size(d,1); %channels
nt = size(d,2); %time
t = (0:nt-1)/fs;

data.label = {'chan1','chan2','chan3','chan4'};
data.time = {t};
data.trial = {d};
data.fsample = fs;


%% Set-up fBOSC parameters

% general setup
cfg.fBOSC.F                 = 2:0.5:20;
cfg.fBOSC.wavenumber        = 6;           % wavelet family parameter (time-frequency tradeoff)
cfg.fBOSC.fsample           = fs;         % current sampling frequency of EEG data

% padding
cfg.fBOSC.pad.tfr_s         = 0.1;      % padding following wavelet transform to avoid edge artifacts in seconds (bi-lateral)
cfg.fBOSC.pad.detection_s   = 0.1;       % padding following rhythm detection in seconds (bi-lateral); 'shoulder' for BOSC eBOSC.detected matrix to account for duration threshold
cfg.fBOSC.pad.background_s  = 0.1;      % padding of segments for BG (only avoiding edge artifacts)

% fooof parameters - fit with fixed line or allow a knee
cfg.fBOSC.fooof.aperiodic_mode    = 'knee';
cfg.fBOSC.fooof.version           = 'python'; %'matlab'

% threshold settings
cfg.fBOSC.threshold.duration	= repmat(6, 1, numel(cfg.fBOSC.F)); % vector of duration thresholds at each frequency (previously: ncyc)
cfg.fBOSC.threshold.percentile  = .99;                              % percentile of background fit for power threshold

% episode post-processing
cfg.fBOSC.postproc.use      = 'no';        % Post-processing turned off for now

% general processing settings
cfg.fBOSC.channel           = []; % select posterior channels (default: all)
cfg.fBOSC.trial             = []; % select trials (default: all)
cfg.fBOSC.trial_background  = []; % select trials for background (default: all)


%%
[fBOSC, cfg] = fBOSC_wrapper(cfg, data);


%% Plot the Results of the 1/f fit
cfg.log_freqs = 1;
cfg.plot_old = 0;
fBOSC_fooof_plot(cfg,fBOSC)


%% fBOSC.detected (channel x trial x freq x time)
chan = 1;
trial = 1;
fidx = cfg.fBOSC.F > 5 & cfg.fBOSC.F < 8;

tmpDetected = single(squeeze(nanmean(fBOSC.detected(chan,trial,fidx,:),3))>0);
tmpDetected(tmpDetected==0) = NaN;
    
origData = data.trial{trial}(cfg.fBOSC.channel(chan),cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);
origData_time = data.time{trial}(cfg.fBOSC.pad.total_sample+1:end-cfg.fBOSC.pad.total_sample);

figure;
plot(origData_time,squeeze(origData), 'k');
hold on;
plot(origData_time,squeeze(origData).*tmpDetected', 'r');
ylabel({'Power';'(a.u.)'});
