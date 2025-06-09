function [ntp,vid,aud,aud_data] = findPupilVidNTP(fname_vid,fname_drift,fname_world)

% fname_vid = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2\Original\Walk1\Pupil\5d738720_0.0-1556.477.mp4";
% fname_rp = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2\Original\Walk1\Raspberry\RP_marks_2021-08-09_10-31-32_871905.txt";
% fname_drift = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2\Original\Walk1\PupilPhoneData\timeDriftRecord-2021-08-09_10-18-49-102.csv";
% fname_world = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2\Original\Walk1\Pupil\world_timestamps.csv";

% Loading video
vid = VideoReader(fname_vid); 
[aud_data.y,aud_data.Fs] = audioread(fname_vid);
aud = audioplayer(aud_data.y,aud_data.Fs);

%get timestamps from csv file (need to adjust for drift)
world_tbl = readtable(fname_world); %in units nanoseconds
drift_tbl = readtable(fname_drift); %in units milliseconds
% mdrift = median(drift_tbl.DriftCorrection)/1000; %median drift in sec?
ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

% a = datetime(ntp_all(1)./(60*60*24),'convertfrom','datenum','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
% b = datetime(world_tbl.timestamp_ns_(1)./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');

ntp = datetime(world_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
ntp = datenum(ntp)*60*60*24; %convert from days to seconds
% ntp = ntp + mdrift; %is this correct? Not needed since drift rarely occurs and is zero for most datasets
ntp = ntp + ntp_offset;



% %%%%%%%%%%%%%%%%%%%%% Find sync using led - not needed %%%%%%%%%%%%%%%%
% % Find sync led in video (need to match green color to improve sensitivity - white lights also have high green component - tried this and sensitivity got worse)
% fr = vid.frameRate;
% nof = vid.NumFrames;
% I = zeros(nof,1);
% wait_msg = parfor_wait(nof);
% for k=1:nof
%     wait_msg.Send;
%     img = read(vid,k);
%     I(k) = mean(reshape(img(1:25,425:475,2),[],1));
% %     image(img(1:25,425:475,2));
% %     drawnow;
% end
% wait_msg.Destroy; 
% 
% %find threshold crossings
% thresh = 225; %RW2 - walk 1
% dthr = diff([false;I>thresh]);
% [ts,tsidx] = find(dthr==1);
% tsend = find(dthr==-1);
% if length(ts)>length(tsend)
%     while 1
%         idx = find((tsend-ts(1:length(tsend)))<0,1);
%         ts(idx) = [];
%         tsidx(idx) = [];
%         if (length(ts)==length(tsend))
%             break;
%         end
%         if isempty(idx)
%             ts(end) = [];
%             tsidx(end) = [];
%         end
%     end
% end
% 
% fid = fopen(fname_rp);
% A = fread(fid,[1,Inf],'*char');
% fclose(fid);
% 
% A = regexp(A,'\n','split');
% A = A(contains(A,'PUPIL'));
% A = regexp(A,'PUPIL_VIDEO: ','split','once');
% A = cat(1,A{:});
% 
% %convolution method to find match
% t1 = datenum(A(:,2))*60*60*24;  %ntp times when led flashed for pupil camera in sec
% t2 = (ts./fr);  %actual led flashes in sec
% 
% nt1 = round(t1*1000); 
% nt2 = round(t2*1000);
% 
% nt1 = nt1 - nt1(1) + 1; %in ms and normalized to 1st value
% nt2 = nt2 - nt2(1) + 1;
% 
% N = max([nt1;nt2]);
% 
% T1 = false(1,N); 
% T1(nt1) = true;
% T1 = imdilate(T1,true(1,2000));
% 
% T2 = false(1,N); 
% T2(nt2) = true;
% T2 = imdilate(T2,true(1,2000));
% 
% [c,lags] = xcorr(T1,T2,'coeff'); %limit to 200ms in both directions
% [~,midx] = max(c);
% shift = lags(midx);
% 
% [~,midx] = min(abs(nt2+shift)); %t2(midx) matches t1(1)
% % midx = find(diff(t2(2:end))>1000,1,'first')+2; %better way to find this index? find 1st led blink after the rapid series of blinks...
% 
% CC = bwconncomp([T1;circshift(T2,shift)]);
% stats = regionprops(CC,'boundingbox');
% bb = cat(1,stats.BoundingBox);
% t1_match_idxs = (bb(:,3)>1000 & bb(:,4)>1); %t1 indices that match t2 starting at midx
% t2_match_idxs = (midx:midx+sum(t1_match_idxs)-1)';
% 
% ntp = t1(t1_match_idxs); %time in sec
% frms = ts(t2_match_idxs); %frame number
% 
% ntp_led = interp1(frms,ntp,1:nof,'linear','extrap')'; %ntp times for every frame in sec
% 
% % figure; 
% % plot(T1.*2,'b'); 
% % hold on; 
% % plot(circshift(T2,shift),'r')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









