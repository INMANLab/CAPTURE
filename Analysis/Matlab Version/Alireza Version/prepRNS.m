function [RNS] = prepRNS(patient_ID)
% 
% prepRNS Imports the .dat files for a patient of interest, 
% marks the timing of TTL pulses and exports a .mat file (RNS.mat)
% with information about the recording(s).
%
% ArgIn: 
%    - patient_ID, unique patient identifier (e.g., RW1, RW2, RW3...) [str]
% ArgOut:
%    - RNS, workspace variable with data separated by walk and .dat file [.mat]
%
% Author:    Justin Campbell
% Contact:   justin.campbell@hsc.utah.edu 
% Version:   05-04-2022

%% Setup environment & RNS params

% get patient ID
RNS.patient_ID = patient_ID;
disp(['Loading data for ', patient_ID]);

% set paths to RWN Box folder and cd into dir
data_path = 'C:\Users\Justin\Box\RealWorldNavigationCory\';
addpath 'C:\Users\Justin\Box\INMANLab\Code\Justin\RealWorldNavigation'
cd(data_path);
cd(patient_ID);

% pull relevant files
RNS.csv = dir('*\*.csv');
RNS.dat_files = dir('*\*\*.dat');
RNS.table = readtable(fullfile(RNS.csv.folder, RNS.csv.name));

% RNS data params
RNS.params = struct;
RNS.params.Fs = 250;
RNS.params.waveform_count = 4;
RNS.params.ttl_size = -512; % fixed TTL size across all recordings for synchronization

% define TTL signal shape for template matching
RNS.params.ttlSignal = [0, -512, -512, 0, 0, -512, -512, 0, 0, 0, 0, 0, 0, -512, -512, 0];
RNS.params.ttlMaxDist = 10000; % sensitivity param for template matching
% findsignal(data, RNS.params.ttlSignal) to preview

% channel map
RNS.chanmap(1) = RNS.table.Ch1Name(1);
RNS.chanmap(2) = RNS.table.Ch2Name(1);
RNS.chanmap(3) = RNS.table.Ch3Name(1);
RNS.chanmap(4) = RNS.table.Ch4Name(1);

% define filters
[b60,a60] = iirnotch(60/(RNS.params.Fs/2), (60/(RNS.params.Fs/2))/25); % 60Hz IIR filter
[b120,a120] = iirnotch(120/(RNS.params.Fs/2), (120/(RNS.params.Fs/2))/25); % 120Hz IIR filter

%% Import & preprocess data

% calculate number of walks
num_walks = length(unique(RNS.table.Walk(~isnan(RNS.table.Walk))));
disp(['Number of walks: ', num2str(num_walks)])

% initialize for storage
RNS.raw_data = cell(num_walks,1);
RNS.processed_data = cell(num_walks,1);
RNS.timevecs = cell(num_walks,1);
RNS.ttls = cell(num_walks,1);

% loop through walks
for i = 1:num_walks
    dat_files = RNS.table(RNS.table.Walk == i, :).Filename; % pull file names
    % initialize for storage
    raw_data = {}; 
    processed_data = {}; 
    ttls = {}; 
    tvecs = {}; 


    % loop through .dat files within walk
    for ii = 1:size(dat_files,1)
        fid = fopen(fullfile(RNS.dat_files(ii).folder, dat_files{ii}));
        dat = fread(fid,'int16');
        fclose(fid);

        % loop through channels
        dat_data=[];
        for iii = 1:4
            dat_data(iii,:) = dat(iii:RNS.params.waveform_count:end)'-512;
        end % chan loop

        % get raw data
        raw_data{ii} = dat_data;

        % find TTLs
        [ttlStart, ttlStop, ttlDist] = findsignal(raw_data{ii}(1,:), RNS.params.ttlSignal, 'MaxDistance', RNS.params.ttlMaxDist); % find TTLs
        if length(ttlDist) == 1;  % if only one TTL was found
            ttls{ii} = ttlStart:ttlStop;
        else
            ttls{ii} = []; % save all the TTLs separately
            for ttl = 1:length(ttlDist)
                ttls{ii} = [ttls{ii}; ttlStart(ttl): ttlStop(ttl)];
            end
        end
        if ttlDist > RNS.params.ttlMaxDist; % eliminate likely false positive matches
            ttls{ii} = [];
        end

        % get time vectors
        tvecs{ii} = 1:length(raw_data{ii});

        % preprocess
        processed_data{ii} = double(raw_data{ii} - median(raw_data{ii}, 1))'; % common median re-reference
        processed_data{ii} = filtfilt(b60, a60, processed_data{ii}); % 60Hz notch
        processed_data{ii} = filtfilt(b120, a120, processed_data{ii}); % 120Hz notch
        processed_data{ii} = processed_data{ii}';

    end % .dat loop
        
    RNS.raw_data{i} = raw_data;
    RNS.processed_data{i} = processed_data;
    RNS.ttls{i} = ttls;
    RNS.timevecs{i} = tvecs;

end % walk loop
