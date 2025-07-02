function [ntp,D,Fs] = findXsensNTP(xs_folder,varargin)

% xs_folder = "D:\Data\RealWorldNavigationCory\RW2\Original\Walk1\Xsens";

xs_xlsx = sortrows(struct2table(dir(fullfile(xs_folder,'*.xlsx'))),'name');
xs_stamps = sortrows(struct2table(dir(fullfile(xs_folder,'xsense_stamps*.txt'))),'name');

ntp = []; D = []; Fs = [];
if (size(xs_xlsx,1)~=size(xs_stamps,1)) || isempty(xs_xlsx)
    mesg = 'findXsensNTP: Number of excel files does not match number of timestamp files OR no excel files exist!';
    if nargin>1
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
    return;
end

for k=1:size(xs_xlsx,1)
    if size(xs_xlsx,1)==1
        fname_xsens = fullfile(xs_xlsx.folder,xs_xlsx.name);
    else
        fname_xsens = fullfile(xs_xlsx.folder{k},xs_xlsx.name{k});
    end
    fname_xsens_header = regexprep(fname_xsens,'\.xlsx$','_header.csv');
    fname_xsens_velocity = regexprep(fname_xsens,'\.xlsx$','_velocity.csv');
    if ~isfile(fname_xsens_header)
        fname_excel2csv = fullfile(fileparts(mfilename('fullpath')),'Excel2Csv.vbs');
        cmdstr = [
            '"',char(fname_excel2csv),'" ',...
            '"',char(fname_xsens),'" ',...
            '"General Information" ',...
            '"',char(fname_xsens_header),'" ',...
            '"Segment Velocity" ',...
            '"',char(fname_xsens_velocity),'" ',...
            ];
        mesg = 'Converting xsens excel file to csv...';
        if nargin>3
            varargin{2}.Message = mesg;
        end
        disp(mesg);
        [status,cmdout] = system(cmdstr,'-echo');
        if status
            mesg = 'Conversion from xlsx to csv failed!';
            if nargin>2
                uialert(varargin{1},mesg,'');
            end
            error(mesg);
        end
    end
end

warning('off','MATLAB:table:ModifiedAndSavedVarnames')

DTable = cell(size(xs_stamps,1),7);
for k=1:size(xs_stamps,1)
    if size(xs_stamps,1)==1
        fname_stamps = fullfile(xs_stamps.folder,xs_stamps.name);
        fname_xsens = fullfile(xs_xlsx.folder,xs_xlsx.name);
    else
        fname_stamps = fullfile(xs_stamps.folder{k},xs_stamps.name{k});
        fname_xsens = fullfile(xs_xlsx.folder{k},xs_xlsx.name{k});
    end
    fname_xsens_header = regexprep(fname_xsens,'\.xlsx$','_header.csv');
    fname_xsens_velocity = regexprep(fname_xsens,'\.xlsx$','_velocity.csv');

    fid = fopen(fname_stamps);
    A = fread(fid,[1,Inf],'*char');
    fclose(fid);

    A = regexp(A,'\n','split');

    if ~isempty(A{1})
        ts_start = A{contains(A,'START')};
        ts_start = regexp(ts_start,'START: ','split','once');
        ts_start = datenum(datetime(ts_start{2},'TimeZone','America/Los_Angeles'))*60*60*24;

        ts_stop = A{contains(A,'STOP')};
        ts_stop = regexp(ts_stop,'STOP: ','split','once');
        ts_stop = datenum(datetime(ts_stop{2},'TimeZone','America/Los_Angeles'))*60*60*24;

        fid = fopen(fname_xsens_header,'r');
        hdr_table = fread(fid,'*char');
        fclose(fid);

        hdr_table = regexp(regexp(hdr_table','\n','split'),',','split');
        hdr_table = [hdr_table{:}];
        hdr_table = reshape(hdr_table(1:end-1),2,[]);
        hdr_table = cell2table(hdr_table(2,:),'VariableNames',hdr_table(1,:));

        vel_table = readtable(fname_xsens_velocity);

        dat = sqrt(sum(vel_table.PelvisX.^2 + vel_table.PelvisY.^2 + vel_table.PelvisZ.^2,2));

        ts_start_hdr = datetime(hdr_table.("Recorded Date UTC"),'TimeZone','UTC');
        ts_start_hdr.TimeZone = 'America/Los_Angeles';
        ts_start_hdr = datenum(ts_start_hdr)*60*60*24;
        Fs = str2double(hdr_table.("Frame Rate"));
        dl = size(dat,1)./Fs; %data length in sec

        DTable(k,1) = {dat};
        DTable(k,2) = {ts_start};
        DTable(k,3) = {ts_stop};
        DTable(k,4) = {dl};
        DTable(k,5) = {size(dat,1)};
        DTable(k,6) = {Fs};
        DTable(k,7) = {ts_start_hdr};
    end
end
DTable = cell2table(DTable,'VariableNames',{'data','ntp','ntp_stop','sec','samp','Fs','ntp_hdr'});

% DTable.ntp_hdr = DTable.ntp_hdr + (DTable.ntp(1)-DTable.ntp_hdr(1)); %adjust to match 1st ntp time (this does not help with the sync problem for multiple segments - aborting)

fs = DTable.Fs;
if ~iscell(fs)
    fs = num2cell(fs);
end
if any(cellfun(@isempty,fs))
    mesg = 'Something is wrong with the xs data! Need to manually check...';
    if nargin>1
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
    return;
else
    if any(DTable.Fs(1)~=DTable.Fs) || any(diff(DTable.ntp)<0)
        mesg = 'Something is wrong with the xs data! Need to manually check...';
        if nargin>1
            uialert(varargin{1},mesg,'');
        end
        disp(mesg);
        return;
    end
end

TotSamp = round(((DTable.ntp(end)+DTable.sec(end)) - DTable.ntp(1))*Fs); %total number of samples based on ntp times
% TotSamp = round(((DTable.ntp_hdr(end)+DTable.sec(end)) - DTable.ntp_hdr(1))*Fs); %total number of samples based on ntp times 

%Stitch data together
D = nan(TotSamp,1);
TS = nan(size(DTable,1),4);
for k=1:size(DTable,1)
    cs = round((DTable.ntp(k)-DTable.ntp(1))*Fs); %current sample
%     cs = round((DTable.ntp_hdr(k)-DTable.ntp_hdr(1))*Fs); %current sample
    dl = DTable.samp(k); %data samples

    D((1:dl)+cs,:) = DTable.data{k};
    TS(k,1) = cs+1;
    TS(k,2) = DTable.ntp(k);
%     TS(k,2) = DTable.ntp_hdr(k);
    TS(k,3) = cs+dl;
    TS(k,4) = DTable.ntp(k)+DTable.sec(k);
%     TS(k,4) = DTable.ntp_hdr(k)+DTable.sec(k);
end

if size(TS,1)==1
    ntp = interp1([TS(1,1);size(D,1)],[TS(1,2);TS(1,2)+(1/Fs)*size(D,1)],1:size(D,1),'linear','extrap')'; 
else
    ntp = interp1(TS(:,1),TS(:,2),1:size(D,1),'linear','extrap')'; %ntp times in sec for all data samples
end

ntp_diff = diff(ntp)-(1/Fs);
if min(ntp_diff)<-0.0001 || max(ntp_diff)>0.0001 %if more than 1/10ms jitter, send alert
    mesg = 'NTP time vector has too much variance! Need to check...';
    if nargin>1
        uialert(varargin{1},mesg,'');
    end
    disp(mesg);
end


