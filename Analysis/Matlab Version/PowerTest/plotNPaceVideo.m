function plotNPaceVideo
% Creates video time synced with neuropace channels. Left/right arrows
% cycle through video frames. Shift + left/right arrows skips through video
% frames. A black vertical line indicates location of video with respect to
% data. Keypress "o" allows user to input an offset factor to align video
% with data. Keypress "r" starts recording of video to root folder and a
% second press of "r" stops recording.
%
% Tyler Davis 20220516

% Selecting files
rootdir = "\\rolstonserver\D\Data\Real World Navigation Cory\RW2";
fname_np = fullfile(rootdir,"NeuroPace_PHI\UCLA_NEA_4344798_ECoG_Catalog.csv");
fname_vid = fullfile(rootdir,"Original\Walk1\Pupil\5d738720_0.0-1556.477.mp4");
fname_rp = fullfile(rootdir,"Original\Walk1\Raspberry\RP_marks_2021-08-09_10-31-32_871905.txt");
root_np = fullfile(rootdir,"NeuroPace_PHI\UCLA_NEA_4344798 Data EXTERNAL #PHI");
fname_world = fullfile(rootdir,"Original\Walk1\Pupil\world_timestamps.csv");
fname_drift = fullfile(rootdir,"Original\Walk1\PupilPhoneData\timeDriftRecord-2021-08-09_10-18-49-102.csv"); %used to adjust/correct times from fname_world
fname_gopro = fullfile(rootdir,"Original\Walk1\Gopro\Gopro_RW2_Walk1_Combined.MP4");


% [file,path] = uigetfile("\\rolstonserver\d\data\Real World Navigation Cory\*.csv","Select the NeuroPace Catalog csv file...");
% fname_np = fullfile(path,file);
% 
% root_np = dir(path);
% idx = [root_np.isdir] & ~(strcmp({root_np.name},'.')|strcmp({root_np.name},'..'));
% root_np = fullfile(root_np(idx).folder,root_np(idx).name);
% 
% [file,path] = uigetfile(path + "*.txt","Select the RaspberryPi RP_Marks txt file...");
% fname_rp = fullfile(path,file);
% 
% [file,path] = uigetfile(path + "*.mp4","Select the PupilLabs mp4 file...");
% fname_vid = fullfile(path,file);
% 
% fname_world = dir(path);
% idx = contains({fname_world.name},'world');
% fname_world = fullfile(fname_world(idx).folder,fname_world(idx).name);
% 
% [file,path] = uigetfile(path + "*.csv","Select the PupilPhone drift file...");
% fname_drift = fullfile(path,file);

%Loading data/timestamps
[ntp_vid,vid] = findPupilVidNTP(fname_vid,fname_drift,fname_world); 
[ntp_gp,vid_gp] = findGoProVidNTP(fname_vid,fname_gopro,ntp_vid);
[ntp_np,D] = findNPaceNTP(fname_np,fname_rp,root_np);
t = datetime(ntp_np./(60*60*24),'convertFrom','datenum'); %convert from seconds to datetime format

%Plotting raster
fH = figure;
aH_data = subplot(1,2,2,'parent',fH);
plot(aH_data,t,D-(0:3)*1000)
ylim(aH_data,[-3500,500])
set(aH_data,'yticklabels',[])
set(fH,'units','normalized','renderer','opengl');

% Plotting video
fr = vid.frameRate;
nof = vid.NumFrames;
offset = 0; %offset in frames
% cf = 1; %current frame
cf = round((ntp_np(1)-ntp_vid(1))*fr); %video starts 1st so fast forward to start of data
[~,cf_gp] = min(abs(ntp_vid(cf-offset)-ntp_gp)); %current gopro frame

warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded')

[~,cs] = min(abs(ntp_vid(cf-offset)-ntp_np)); %current data sample
hold(aH_data,'on');
pH = plot(aH_data,[t(cs),t(cs)],[-3500,500],'k');

aH_vid = subplot(1,2,1,'parent',fH);
img = read(vid,cf);
iH = image(img,'parent',aH_vid);
axis(aH_vid,'image')
tt = t(cs); tt.Format = 'dd-MMM-uuuu HH:mm:ss.SSS';
title(aH_vid,string(tt))

% Saving data to figure to access in updatePlot callback
SS = [];
SS.pH = pH;
SS.vid = vid;
SS.vid_gp = vid_gp;
SS.fr = fr;
SS.nof = nof;
SS.cf = cf;
SS.cf_gp = cf_gp;
SS.iH = iH;
SS.offset = offset;
SS.rf = false; %recording flag
SS.writerObj = [];
SS.aH_data = aH_data;
SS.aH_vid = aH_vid;
SS.ntp_vid = ntp_vid;
SS.ntp_np = ntp_np;
SS.ntp_gp = ntp_gp;
SS.t = t;
SS.cs = cs;
SS.gpflag = false;

setappdata(fH,'SS',SS)

set(fH,'windowkeypressfcn',@keypressCallback)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Keypress callback function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypressCallback(varargin)

keystr = varargin{1, 2}.Key;
if isempty(varargin{1, 2}.Modifier)
    modstr = [];
else
    modstr = varargin{1, 2}.Modifier{1,1};
end
fH = varargin{1};
SS = getappdata(fH,'SS');

switch keystr
    case 'rightarrow'
        if strcmp(modstr,'shift')
            SS.cf = SS.cf+120;
        else
            SS.cf = SS.cf+1;
        end
        SS = updateFcn(SS);
    case 'leftarrow'
        if strcmp(modstr,'shift')
            SS.cf = SS.cf-120;
        else
            SS.cf = SS.cf-1;
        end
        SS = updateFcn(SS);
    case 'o'
        SS.offset = inputdlg(sprintf('Enter #Frames Offset (%0.3f)',SS.offset));
        SS.offset = str2double(SS.offset{:});
        SS = updateFcn(SS);
    case 'g' %go to location
        %10:32:36.985
        tt = SS.t(SS.cs); tt.Format = 'HH:mm:ss.SSS';
        ts = inputdlg(sprintf('Go to timestamp (%s)',string(tt)));
        tt.Format = 'dd-MMM-uuuu';
        ts_sec = datenum(datetime([char(tt),' ',ts{1}]))*60*60*24;
        [~,SS.cf] = min(abs(SS.ntp_vid-ts_sec));
        SS = updateFcn(SS);
    case 's' %switch to gopro
        SS.gpflag = ~SS.gpflag;
        SS = updateFcn(SS);
    case 'r' %record (not working)
%         if SS.rf
%             SS.rf = false;
%             set(pH,'color','k')
%             if isobject(SS.writerObj)
%                 close(SS.writerObj);
%                 delete(SS.writerObj);
%             end
%         else
%             SS.rf = true;
%             set(pH,'color','r')
%             SS.writerObj = VideoWriter(datestr(clock,'yyyymmdd-HHMMSS'),'MPEG-4');
%             SS.writerObj.FrameRate = SS.fr;
%             SS.writerObj.Quality = 100;
%             open(SS.writerObj);
%         end
%         setappdata(fH,'SS',SS)
%         while 1
%             try
%                 SS.cf = SS.cf+1;
%                 clc; fprintf('recording %0.0f\n',SS.cf);
%                 SS = updateFcn(SS);
%                 writeVideo(SS.writerObj,getframe(fH));
%                 xwin = get(SS.aH_data,'xlim');
%                 if SS.t(SS.cs)>xwin(2)
%                     os = SS.t(SS.cs) - xwin(1);
%                     set(SS.aH_data,'xlim',xwin + os)
%                 end
%             catch
%                 break;
%             end
%         end
end %switch
setappdata(fH,'SS',SS)




%%%%%%%%%%%%%%%%%%%%%%%%% Update frame rate variables %%%%%%%%%%%%%%%%%%%%%
function SS = updateFcn(SS)

SS.cf(SS.cf<1) = 1;
SS.cf(SS.cf>SS.nof) = SS.nof;

[~,SS.cs] = min(abs(SS.ntp_vid(SS.cf-SS.offset)-SS.ntp_np)); %current data sample
[~,SS.cf_gp] = min(abs(SS.ntp_vid(SS.cf-SS.offset)-SS.ntp_gp)); %current gopro frame

xwin = get(SS.aH_data,'xlim');
if SS.t(SS.cs)>xwin(2)
    os = SS.t(SS.cs) - xwin(1);
    set(SS.aH_data,'xlim',xwin + os)
end

if SS.gpflag
    img = read(SS.vid_gp,SS.cf_gp);
else
    img = read(SS.vid,SS.cf);
end
set(SS.iH,'cdata',img)
set(SS.pH,'xdata',[SS.t(SS.cs),SS.t(SS.cs)])
tt = SS.t(SS.cs); tt.Format = 'dd-MMM-uuuu HH:mm:ss.SSS';
title(SS.aH_vid,string(tt))





