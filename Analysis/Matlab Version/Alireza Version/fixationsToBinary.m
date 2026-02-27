function fixMask = fixationsToBinary(time, Fs, fixTable)
% Convert fixation intervals (in time) to binary mask
%
% Inputs:
%   time     : Nx1 time vector (seconds)
%   Fs       : sampling rate (Hz) — used only for safety checks
%   fixTable : table with start_t and end_t (seconds)
%
% Output:
%   fixMask  : Nx1 double array (1 = fixation, 0 = not)

N = length(time);
fixMask = false(N,1);

if isempty(fixTable)
    fixMask = double(fixMask);
    return
end

for k = 1:height(fixTable)

    t0 = fixTable.start_t(k);
    t1 = fixTable.end_t(k);

    % find nearest sample indices
    i0 = max(1, floor((t0 - time(1))*Fs) + 1);
    i1 = min(N, ceil((t1 - time(1))*Fs) + 1);

    fixMask(i0:i1) = true;
end

fixMask = double(fixMask);
end
