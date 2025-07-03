function updateRWNAppFiles(RootDir)
%Adds new data streams to RWNApp*.mat files without using the GUI.

RWNAppFiles = dir([RootDir,'\RWNApp*.mat']);
OrigDir = '\\155.100.91.44\d\Data\RealWorldNavigationCory';

for k=1:length(RWNAppFiles)
    disp(k/length(RWNAppFiles));

    rwnappfile = fullfile(RWNAppFiles(k).folder,RWNAppFiles(k).name);
    ptname = regexp(RWNAppFiles(k).name,'RW\d','match','once');
    walkname = regexp(RWNAppFiles(k).name,'Walk\d','match','once');
    origdir = fullfile(OrigDir,ptname,'Original',walkname);
    % vnames = who('-file',rwnappfile);
    Files = findRWAFiles(origdir);

    disp(rwnappfile);
    disp(origdir);

    %%%%%%%%% Updating evnts_tbl in destination using original %%%%%%%%%% 
    %Kirsten updated the events in the original files. Copying these events
    %to the files in D:\Tyler\RealWorld on the server. Then, new datasets
    %like biopac, gaze, etc. do not need to be added to these files.
    orig_data = load(Files.save_file);
    dest_data = load(rwnappfile);

    orig_fields = fieldnames(orig_data);
    same_flag = false(length(orig_fields),1);
    for m=1:length(orig_fields)
        orig_var = orig_data.(orig_fields{m});
        dest_var = dest_data.(orig_fields{m});
        if all(size(orig_var)==size(dest_var))
            same_flag(m) = true;
        else
            fprintf('%s does not match for %s %s!\n',orig_fields{m},ptname,walkname)
        end
    end

    dest_data.evnts_tbl = orig_data.evnts_tbl;
    fprintf('Overwriting evnts_tbl in destination for %s %s!\n',ptname,walkname)
    save(rwnappfile,'-struct','dest_data','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ntp_gaze = []; d_gaze_x = []; d_gaze_y = []; d_gaze_fix = []; fs_gaze = [];
    % if ~(isempty(Files.drift_csv_file)||isempty(Files.pupil_gaze_file))
    %     [ntp_gaze,d_gaze_x,d_gaze_y,d_gaze_fix,fs_gaze] = findPupilGazeNTP(Files);
    % end
    % save(rwnappfile,'ntp_gaze','d_gaze_x','d_gaze_y','d_gaze_fix','fs_gaze','-append');
    % 
    % ntp_amb = []; d_amb = []; fs_amb = [];
    % if ~(isempty(Files.chest_drift_csv_file)||isempty(Files.chest_ambient_csv_file))
    %     [ntp_amb,d_amb,fs_amb] = findChestAmbNTP(Files);
    % end
    % save(rwnappfile,'ntp_amb','d_amb','fs_amb','-append');
    % 
    % ntp_kde = []; d_kde = []; fs_kde = [];
    % if ~(isempty(Files.kde_file)||isempty(Files.frame_comp_file))
    %     [ntp_kde,d_kde,fs_kde] = findKDEsNTP(Files);
    % end
    % save(rwnappfile,'ntp_kde','d_kde','fs_kde','-append');
    % 
    % ntp_bio = []; d_bio = []; fs_bio = [];
    % if ~(isempty(Files.biopac_data_file)||isempty(Files.chest_drift_csv_file))
    %     [ntp_bio,d_bio,fs_bio] = findBiopacNTP(Files);
    % end
    % save(rwnappfile,'ntp_bio','d_bio','fs_bio','-append');

end