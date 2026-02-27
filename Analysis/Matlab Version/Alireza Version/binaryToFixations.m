function fixTable = binaryToFixations(fixMask, time, Fs, minDur_s)
% Convert binary fixation mask to fixation table
%
% Inputs:
%   fixMask  : Nx1 binary/logical array (1 = fixation)
%   time     : Nx1 time vector (seconds)
%   Fs       : sampling rate (Hz)
%   minDur_s : optional minimum fixation duration (sec), e.g. 0.06–0.1
%
% Output:
%   fixTable : table of fixation intervals

if nargin < 4
    minDur_s = 0;   % no duration filtering by default
end

fixMask = fixMask(:) ~= 0;
time = time(:);
N = length(fixMask);

% --- find contiguous runs of ones ---
d = diff([false; fixMask; false]);
idx0 = find(d == 1);
idx1 = find(d == -1) - 1;

% --- compute times and durations ---
start_t = time(idx0);
end_t   = time(idx1);
dur_s   = (idx1 - idx0 + 1) / Fs;

% --- apply minimum duration filter ---
keep = dur_s >= minDur_s;

idx0 = idx0(keep);
idx1 = idx1(keep);
start_t = start_t(keep);
end_t   = end_t(keep);
dur_s   = dur_s(keep);

% --- build table ---
fixTable = table(start_t, end_t, dur_s, idx0, idx1);

end
