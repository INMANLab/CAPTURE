function [ntp,gazeX,gazeY,fixD,Fs] = findPupilGazeNTP(Files,varargin)

%get timestamps from csv file (need to adjust for drift)
drift_tbl = readtable(Files.drift_csv_file); %in units milliseconds
ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

orig_flag = true;
if orig_flag
    gaze_tbl = readtable(Files.pupil_gaze_file,'Delimiter',',');
else
    gaze_tbl = readtable(Files.gaze_csv_file,'Delimiter',','); %reprocessed
end

% a = datetime(ntp_all(1)./(60*60*24),'convertfrom','datenum','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
% b = datetime(world_tbl.timestamp_ns_(1)./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');

ntp_orig = datetime(gaze_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
ntp_orig = datenum(ntp_orig)*60*60*24; %convert from days to seconds
ntp_orig = ntp_orig + ntp_offset;

gazeX_orig = gaze_tbl.gazeX_px_;
gazeY_orig = gaze_tbl.gazeY_px_;

Fs = 200; %there is some spiking in the sampling rate but mean is close to 200 (resampling to make consistent)
N = round((ntp_orig(end)-ntp_orig(1))*Fs); %number of samples based on ntp times and Fs
ntp = ((0:N-1)/Fs + ntp_orig(1))'; %time vector in sec that is evenly sampled at Fs and adjusted to match start of ntp

gazeX = interp1(ntp_orig,gazeX_orig,ntp,'linear','extrap');
gazeY = interp1(ntp_orig,gazeY_orig,ntp,'linear','extrap');
if orig_flag
    fixD = calcRWAFixations(gazeX,gazeY,Fs,20,500); %gazeX, gazeY, sampleRate (Hz), maxDispersion (pixels), minDuration (ms)
else
    fixD = interp1(ntp_orig,gaze_tbl.fixationId,ntp,'nearest','extrap'); %these are fixation counts so nearest is needed
end







