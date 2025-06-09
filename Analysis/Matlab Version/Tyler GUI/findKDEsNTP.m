function [ntp,D,Fs] = findKDEsNTP(Files,varargin)
%finds ntp times for student identified events (KDEs)

rwn_vars = load(Files.save_file);
kde_tbl = readtable(Files.kde_file);
frame_tbl = readtable(Files.frame_comp_file);

pt = regexp(Files.walk_dir,'RW(\d)','tokens','once'); pt = [pt{:}];
wk = regexp(Files.walk_dir,'Walk(\d)','tokens','once'); wk = [wk{:}];

frame_tbl = frame_tbl(strcmp(frame_tbl.Walk,[pt,'_',wk]),:);
ntp_offset = frame_tbl.CutFramesBeginning./frame_tbl.fpsOriginal; %offset from true ntp time in sec

ntp_orig = kde_tbl.X + ntp_offset + rwn_vars.ntp_gp(1);
D_orig = kde_tbl.Y;
% Fs = frame_tbl.fpsSynced; %there is some spiking in the sampling rate - resampling to make consistent
Fs = 60; %make fixed so its the same for all datasets

N = round((ntp_orig(end)-ntp_orig(1))*Fs); %number of samples based on ntp times and Fs
ntp = ((0:N-1)/Fs + ntp_orig(1))'; %time vector in sec that is evenly sampled at Fs and adjusted to match start of ntp
D = interp1(ntp_orig,D_orig,ntp,'linear','extrap'); 








