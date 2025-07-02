function [ntp,D,Fs] = findChestAccelNTP(Files,varargin)
%Accel data from "ChestPhoneData" folder. 

%get timestamps from csv file (need to adjust for drift)
drift_tbl = readtable(Files.chest_drift_csv_file); %in units milliseconds
ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

accel_tbl = readtable(Files.chest_accel_csv_file,'Delimiter',',');

% a = datetime(ntp_all(1)./(60*60*24),'convertfrom','datenum','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
% b = datetime(world_tbl.timestamp_ns_(1)./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');

ntp_orig = datetime(accel_tbl.Var1,'timezone','America/Los_Angeles','Format','uuuu-MM-dd_HH-mm-ss-SSS');
ntp_orig = datenum(ntp_orig)*60*60*24; %convert from days to seconds
ntp_orig = ntp_orig + ntp_offset;

D_orig = [accel_tbl.Var2,accel_tbl.Var3,accel_tbl.Var4];

[uNTP,~,uNTPIdx] = unique(ntp_orig);

%Removing duplicate ntp times
rm_flag = false(size(ntp_orig));
for k=1:length(uNTP)
    idx = (uNTPIdx==k);
    if sum(idx)>1
        idx(find(idx,1)) = false;
        rm_flag(idx) = true;
    end
end
ntp_orig(rm_flag) = [];
D_orig(rm_flag,:) = [];

Fs = 100; 
N = round((ntp_orig(end)-ntp_orig(1))*Fs); %number of samples based on ntp times and Fs
ntp = ((0:N-1)/Fs + ntp_orig(1))'; %time vector in sec that is evenly sampled at Fs and adjusted to match start of ntp

D = interp1(ntp_orig,D_orig,ntp,'linear','extrap'); 







