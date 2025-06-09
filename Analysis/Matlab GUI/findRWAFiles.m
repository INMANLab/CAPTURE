function Files = findRWAFiles(walk_dir, varargin)

% walk_dir = '\\155.100.91.44\D\Data\RealWorldNavigationCory\RW2\Original\Walk1'

pt_dir = walk_dir(1:regexp(walk_dir,'RW\d+','end'));

%file for saving gui data
pt = regexp(walk_dir,'RW\d','match','once');
wk = regexp(walk_dir,'Walk\d','match','once');
save_file = fullfile(walk_dir,sprintf('RWNApp_%s_%s.mat',pt,wk));

if contains(pt,'RW1')
    np_data_dir = fullfile(walk_dir,"RNS");
    np_csv_file = [];
    photon_export_file = dir(fullfile(walk_dir,'photon_export.mat'));
    if ~isempty(photon_export_file)
        photon_export_file = fullfile(photon_export_file(1).folder,photon_export_file(1).name);
    end
else
    photon_export_file = [];

    %neuropace
    np_dir = dir(fullfile(pt_dir,"NeuroPace_PHI"));
    np_data_dir = np_dir(find(contains({np_dir.name},'data','IgnoreCase',true) & [np_dir.isdir],1));
    np_data_dir = fullfile(np_data_dir.folder,np_data_dir.name);

    np_csv_file = np_dir(find(contains({np_dir.name},'.csv'),1));
    np_csv_file = fullfile(np_csv_file.folder,np_csv_file.name);
end

%pupil
pupil_dir = dir(fullfile(walk_dir,"Pupil/**/*"));
pupil_ts_file = pupil_dir(find(contains({pupil_dir.name},'world_timestamps.csv'),1));
pupil_ts_file = fullfile(pupil_ts_file.folder,pupil_ts_file.name);

pupil_vid_file = pupil_dir(find(contains({pupil_dir.name},'.mp4'),1));
pupil_vid_file = fullfile(pupil_vid_file.folder,pupil_vid_file.name);

%raspberry
rp_dir = dir(fullfile(walk_dir,"Raspberry"));
rp_marks_file = rp_dir(find(contains({rp_dir.name},'marks') & contains({rp_dir.name},'.txt'),1));
if isempty(rp_marks_file)
    rp_marks_file = '';
else
    rp_marks_file = fullfile(rp_marks_file.folder,rp_marks_file.name);
end

%drift (Applies only to data in "Pupil" and "PupilPhoneData" folders?)
drift_dir = dir(fullfile(walk_dir,"PupilPhoneData"));
drift_csv_file = drift_dir(find(contains({drift_dir.name},'drift','IgnoreCase',true) & contains({drift_dir.name},'.csv'),1));
drift_csv_file = fullfile(drift_csv_file.folder,drift_csv_file.name);

%chest drift (Applies only to data in "ChestPhoneData" folder?)
chest_drift_dir = dir(fullfile(walk_dir,"ChestPhoneData"));
chest_drift_csv_file = chest_drift_dir(find(contains({chest_drift_dir.name},'drift','IgnoreCase',true) & contains({chest_drift_dir.name},'.csv'),1));
if isempty(chest_drift_csv_file)
    chest_drift_csv_file = '';
else
    chest_drift_csv_file = fullfile(chest_drift_csv_file.folder,chest_drift_csv_file.name);
end

%ambient light from pupil phone (bad data)
ambient_csv_file = drift_dir(find(contains({drift_dir.name},'ambient','IgnoreCase',true) & contains({drift_dir.name},'.csv'),1));
if isempty(ambient_csv_file)
    ambient_csv_file = '';
else
    ambient_csv_file = fullfile(ambient_csv_file.folder,ambient_csv_file.name);
end

%ambient light from chest phone (data looks better)
chest_ambient_csv_file = chest_drift_dir(find(contains({chest_drift_dir.name},'ambient','IgnoreCase',true) & contains({chest_drift_dir.name},'.csv'),1));
if isempty(chest_ambient_csv_file)
    chest_ambient_csv_file = '';
else
    chest_ambient_csv_file = fullfile(chest_ambient_csv_file.folder,chest_ambient_csv_file.name);
end

%accel data from pupil phone
accel_csv_file = drift_dir(find(contains({drift_dir.name},'accel','IgnoreCase',true) & contains({drift_dir.name},'.csv'),1));
if isempty(accel_csv_file)
    accel_csv_file = '';
else
    accel_csv_file = fullfile(accel_csv_file.folder,accel_csv_file.name);
end

%accel data from chest phone (used for drop sync with biopac)
chest_accel_csv_file = chest_drift_dir(find(contains({chest_drift_dir.name},'accel','IgnoreCase',true) & contains({chest_drift_dir.name},'.csv'),1));
if isempty(chest_accel_csv_file)
    chest_accel_csv_file = '';
else
    chest_accel_csv_file = fullfile(chest_accel_csv_file.folder,chest_accel_csv_file.name);
end

%gaze data from pupil phone
pupil_gaze_file = pupil_dir(find(contains({pupil_dir.name},'gaze','IgnoreCase',true),1));
if isempty(pupil_gaze_file)
    pupil_gaze_file = '';
else
    pupil_gaze_file = fullfile(pupil_gaze_file.folder,pupil_gaze_file.name);
end

%imu data from pupil phone
pupil_imu_file = pupil_dir(find(contains({pupil_dir.name},'imu','IgnoreCase',true),1));
if isempty(pupil_imu_file)
    pupil_imu_file = '';
else
    pupil_imu_file = fullfile(pupil_imu_file.folder,pupil_imu_file.name);
end

%gaze (reprocessed to include fixations)
gaze_dir = dir(fullfile(walk_dir,"GLM_Gaze"));
gaze_csv_file = gaze_dir(find(contains({gaze_dir.name},'gaze','IgnoreCase',true) & contains({gaze_dir.name},'.csv'),1));
if isempty(gaze_csv_file)
    gaze_csv_file = '';
else
    gaze_csv_file = fullfile(gaze_csv_file.folder,gaze_csv_file.name);
end

%gopro
gopro_dir = dir(fullfile(walk_dir,"Gopro"));
gopro_vid_file = gopro_dir(find(contains({gopro_dir.name},'combined','IgnoreCase',true) & contains({gopro_dir.name},'.mp4','IgnoreCase',true),1));
if isempty(gopro_vid_file)
    disp('GoPro combined video is missing! Combining individual files now...')
    if nargin>1
        varargin{1}.Message = 'GoPro combined video is missing! Combining individual files now...';
    end
    combineGoProVids(fullfile(walk_dir,"Gopro"));
    gopro_dir = dir(fullfile(walk_dir,"Gopro"));
    gopro_vid_file = gopro_dir(find(contains({gopro_dir.name},'combined','IgnoreCase',true) & contains({gopro_dir.name},'.mp4','IgnoreCase',true),1));
end
gopro_vid_file = fullfile(gopro_vid_file.folder,gopro_vid_file.name);

%gopro synced
gopro_synced_dir = dir(fullfile(regexprep(walk_dir,'Original','Synced'),"Gopro"));
gopro_synced_vid_file = gopro_synced_dir(find(contains({gopro_synced_dir.name},'video','IgnoreCase',true) & contains({gopro_synced_dir.name},'.mp4','IgnoreCase',true),1));
if isempty(gopro_synced_vid_file)
    disp('Synced GoPro combined video is missing!')
    gopro_synced_vid_file = '';
else
    gopro_synced_vid_file = fullfile(gopro_synced_vid_file.folder,gopro_synced_vid_file.name);
end

%xsens
xsens_dir = fullfile(walk_dir,'Xsens');

%frame comparison table for synced versus original gopro
frame_comp_file = fullfile(fileparts(pt_dir),'FrameComparison.xlsx');

%KDE directory (undergrad event times for synced gopro)
ptnum = regexp(pt,'\d$','match','once');
wknum = regexp(wk,'\d$','match','once');
kde_dir = dir(fullfile(fileparts(pt_dir),'KDEs',['RWN',ptnum,'_',wknum,'_KDE.csv']));
if isempty(kde_dir)
    disp('kde_file is missing!');
    kde_file = '';
else
    kde_file = fullfile(kde_dir(1).folder,kde_dir(1).name);
end

%biopac
biopac_dir = dir(fullfile(walk_dir,'Biopac'));
biopac_header_file = biopac_dir(contains({biopac_dir.name},'header','IgnoreCase',true));
if isempty(biopac_header_file)
    biopac_header_file = '';
else
    biopac_header_file = fullfile(biopac_header_file.folder,biopac_header_file.name);
end
biopac_data_file = biopac_dir(contains({biopac_dir.name},'cleanscr','IgnoreCase',true));
if isempty(biopac_data_file)
    biopac_data_file = '';
else
    biopac_data_file = fullfile(biopac_data_file.folder,biopac_data_file.name);
end


%save to structure
Files = [];

Files.walk_dir = walk_dir;
Files.pt_dir = pt_dir;
Files.np_data_dir = np_data_dir;
Files.np_csv_file = np_csv_file;
Files.pupil_ts_file = pupil_ts_file;
Files.pupil_vid_file = pupil_vid_file;
Files.rp_marks_file = rp_marks_file;
Files.drift_csv_file = drift_csv_file;
Files.gopro_vid_file = gopro_vid_file;
Files.xsens_dir = xsens_dir;
Files.save_file = save_file;
Files.photon_export_file = photon_export_file;
Files.ambient_csv_file = ambient_csv_file; %pupil phone
Files.chest_ambient_csv_file = chest_ambient_csv_file; %chest phone
Files.chest_accel_csv_file = chest_accel_csv_file; %chest phone
Files.accel_csv_file = accel_csv_file; %pupil phone
Files.gaze_csv_file = gaze_csv_file; %reprocessed
Files.pupil_gaze_file = pupil_gaze_file; %pupil phone unprocessed (no fixations)
Files.pupil_imu_file = pupil_imu_file; %pupil phone imu data
Files.chest_drift_csv_file = chest_drift_csv_file; %chest phone
Files.frame_comp_file = frame_comp_file;
Files.kde_file = kde_file;
Files.gopro_synced_vid_file = gopro_synced_vid_file;
Files.biopac_header_file = biopac_header_file;
Files.biopac_data_file = biopac_data_file;









