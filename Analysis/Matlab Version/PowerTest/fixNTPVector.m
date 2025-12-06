function fixNTPVector
%Not going to worry about this right now. The pupil ntp vector is not
%monotonic either, so maybe it doesn't matter for the np vector as long as
%the variance from 1/fs is small.

[path,name,ext] = lastPath('RWNApp*.mat','Choose RWN file...');

SS = load(fullfile(path,[name,ext]));

SS_names = fieldnames(SS);
for k=1:length(SS_names)
    ss_name = SS_names{k};
    if regexp(ss_name,'^ntp_','once')
        ss_parts = regexp(ss_name,'_','split','once');
        fs_name = SS_names{~cellfun(@isempty,regexp(SS_names,['(fs|fr)_',ss_parts{2}],'match','once'))};
        t = SS.(ss_name);
        fs = SS.(fs_name);
        tsamp = floor(t.*fs); tsamp = tsamp - tsamp(1) + 1; %these should be monotonically increasing samples
        tsec = ((0:length(t)-1)./fs+t(1))'; %time in sec comparable to t but using monotonically increasing samples
        diff_tsamp = diff(tsamp);
        if any(abs(t-tsec)>(1/fs))
            fprintf(['Neuropace NTP time vector for %s has some jitter.\nMin/max jitter from expected value of 1 is %0.0f/%0.0f samples.' ...
                '\nTotal # jitter values %0.0f/%0.0f.\nMax deviation from monotonic is %0.2e sec.\n'],obj.PatientID,min(diff_tsamp),...
                max(diff_tsamp),sum(diff_tsamp~=1),length(diff_tsamp),max(abs(t-tsec)));
        end
    end
end