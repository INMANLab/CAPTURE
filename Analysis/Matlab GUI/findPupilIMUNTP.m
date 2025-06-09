function [ntp,IMU,Fs] = findPupilIMUNTP(Files,varargin)

%get timestamps from csv file (need to adjust for drift)
drift_tbl = readtable(Files.drift_csv_file); %in units milliseconds
ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

imu_tbl = readtable(Files.pupil_imu_file,'Delimiter',',');

% a = datetime(ntp_all(1)./(60*60*24),'convertfrom','datenum','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
% b = datetime(world_tbl.timestamp_ns_(1)./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');

ntp_orig = datetime(imu_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
ntp_orig = datenum(ntp_orig)*60*60*24; %convert from days to seconds
ntp_orig = ntp_orig + ntp_offset;

gyroX_orig = imu_tbl.gyroX_deg_s_;
gyroY_orig = imu_tbl.gyroY_deg_s_;
gyroZ_orig = imu_tbl.gyroZ_deg_s_;
accelX_orig = imu_tbl.accelerationX_G_;
accelY_orig = imu_tbl.accelerationY_G_;
accelZ_orig = imu_tbl.accelerationZ_G_;
roll_orig = imu_tbl.roll;
pitch_orig = imu_tbl.pitch;

Fs = 200; %there is some spiking in the sampling rate but mean is close to 200 (resampling to make consistent)
N = round((ntp_orig(end)-ntp_orig(1))*Fs); %number of samples based on ntp times and Fs
ntp = ((0:N-1)/Fs + ntp_orig(1))'; %time vector in sec that is evenly sampled at Fs and adjusted to match start of ntp

IMU = [];
IMU.gyroX = interp1(ntp_orig,gyroX_orig,ntp,'linear','extrap');
IMU.gyroY = interp1(ntp_orig,gyroY_orig,ntp,'linear','extrap');
IMU.gyroZ = interp1(ntp_orig,gyroZ_orig,ntp,'linear','extrap');
IMU.accelX = interp1(ntp_orig,accelX_orig,ntp,'linear','extrap');
IMU.accelY = interp1(ntp_orig,accelY_orig,ntp,'linear','extrap');
IMU.accelZ = interp1(ntp_orig,accelZ_orig,ntp,'linear','extrap');
IMU.roll = interp1(ntp_orig,roll_orig,ntp,'linear','extrap');
IMU.pitch = interp1(ntp_orig,pitch_orig,ntp,'linear','extrap');






