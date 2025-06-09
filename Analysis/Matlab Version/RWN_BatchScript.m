%%
ptlist = 2:5;
RWA = RWAnalysis;
for k=2:length(ptlist)
    fprintf('Patient %0.0f\n',ptlist(k));
    RWA.loadData('PatientIdx',ptlist(k));
    RWA.plotInOutDiff;
    RWA.plotInOutDiff('fpass',[2,85]);
    RWA.plotTransSpecGram;
    RWA.batchWalkSpecGram;
end

%% Still need to run this using updated findNPaceNTP with segment shifts for RW4, and RW5 (20230313)
root_list = {...
    'D:\Data\RealWorldNavigationCory\RW2\Original',...
    'D:\Data\RealWorldNavigationCory\RW3\Original',...
    'D:\Data\RealWorldNavigationCory\RW4\Original',...
    'D:\Data\RealWorldNavigationCory\RW5\Original',...
    };

warning('off','MATLAB:table:ModifiedAndSavedVarnames');

for k=1:length(root_list)
    rt = root_list{k};
    fname_np = dir(fullfile(fileparts(rt),'NeuroPace_PHI\*.csv'));
    fname_np = fullfile(fname_np(1).folder,fname_np(1).name);
    root_np = dir(fullfile(fileparts(rt),'NeuroPace_PHI\*PHI'));
    root_np = fullfile(root_np(1).folder,root_np(1).name);
    wk = dir([rt,'\**\RWNApp_RW?_Walk?.mat']);
    for m=1:length(wk)
        disp([k,m]);

        Files = findRWAFiles(wk(m).folder);
        fname_mat = fullfile(wk(m).folder,wk(m).name);
        fname_rp = dir([wk(m).folder,'\Raspberry\RP*.txt']);
        fname_rp = fullfile(fname_rp(1).folder,fname_rp(1).name);
        fname_xs = fullfile(wk(m).folder,'Xsens');

        app_data = load(fname_mat);

        %get timestamps from csv file (need to adjust for drift)
        world_tbl = readtable(Files.pupil_ts_file); %in units nanoseconds
        drift_tbl = readtable(Files.drift_csv_file); %in units milliseconds
        ntp_offset = drift_tbl.CurrNTPOffset(1)/1000; %offset from true ntp time in sec

        ntp = datetime(world_tbl.timestamp_ns_./1e9,'convertfrom','posixtime','timezone','America/Los_Angeles','Format','dd-MMM-uuuu HH:mm:ss.SSS');
        ntp = datenum(ntp)*60*60*24; %convert from days to seconds
        ntp = ntp + ntp_offset;

        if any(size(app_data.ntp_pupil)~=size(ntp))
            error('ntp length does not match ntp_pupil');
        end
        app_data.ntp_pupil = ntp;
        app_data.ntp_gp = app_data.ntp_gp + ntp_offset;
        app_data.ntp_gp_orig = app_data.ntp_gp_orig + ntp_offset;

        app_data.ntp_np = app_data.ntp_np_orig;
        app_data.offset_np = 0;

        app_data.ntp_xs = app_data.ntp_xs_orig;
        app_data.offset_xs = 0;

        for n=1:size(app_data.evnts_tbl,1)
            app_data.evnts_tbl.NTP(n) = app_data.ntp_pupil(app_data.evnts_tbl.PupilFrame(n));
            [~,gp_frame] = min(abs(app_data.evnts_tbl.NTP(n)-app_data.ntp_gp));
            if gp_frame~=app_data.evnts_tbl.GoProFrame(n)
                disp('GP frames do not match!');
            end
            app_data.evnts_tbl.GoProFrame(n) = gp_frame;
            [~,app_data.evnts_tbl.NPSample(n)] = min(abs(app_data.evnts_tbl.NTP(n)-app_data.ntp_np));
        end

        save(fname_mat,'-struct','app_data','-append');


%         [ntp_np_orig,d_np,Fs] = findNPaceNTP(char(fname_np),fname_rp,root_np);
%         ntp_np = ntp_np_orig;
%         save(fname_mat,'ntp_np','ntp_np_orig','d_np','-append');

%         [ntp_xs_orig,d_xs,fs_xs] = findXsensNTP(fname_xs);
%         ntp_xs = ntp_xs_orig;
%         offset_xs = 0;
%         save(fname_mat,'ntp_xs','ntp_xs_orig','d_xs','offset_xs','fs_xs','-append');

        close all;
    end
end



