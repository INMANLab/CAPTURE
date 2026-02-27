function events = FixationCalculation(t, xPix, yPix, fs, pxPerDeg, mdur)
% Method used IVT
% Standard I-VT (Velocity-Threshold) eye movement event detection.
%
% Inputs:
%   t         : Nx1 time in seconds (or empty -> uses sample index/fs)
%   xPix,yPix : Nx1 gaze position in pixels (NaNs allowed)
%   fs        : sampling rate (Hz)
%   pxPerDeg  : pixels-per-degree conversion (from your setup; e.g., ~30-60)
%
% Output (struct):
%   events.fixations : table [start_t, end_t, dur_s, x_mean, y_mean, idx0, idx1]
%   events.saccades  : table [start_t, end_t, dur_s, amp_deg, peakVel_dps, idx0, idx1]

% --------------------------
% Parameters (typical defaults)
% --------------------------
params.vth_dps      = 30;     % velocity threshold for saccade (deg/s); common 20-40
params.minFix_s     = mdur;   % min fixation duration (s); common 60-100 ms
params.minSac_s     = 0.010;  % min saccade duration (s); common 10 ms
params.smooth_win_s = 0.02;  % smoothing window (s); common 5-15 ms

N = numel(xPix);
if isempty(t)
    t = (0:N-1)'/fs;
else
    t = t(:);
end
xPix = xPix(:); yPix = yPix(:);

% --------------------------
% Handle missing samples
% --------------------------
valid = isfinite(xPix) & isfinite(yPix) & isfinite(t);
x = xPix; y = yPix;
% light interpolation across short gaps (optional but common)
x(~valid) = interp1(t(valid), xPix(valid), t(~valid), 'linear', 'extrap');
y(~valid) = interp1(t(valid), yPix(valid), t(~valid), 'linear', 'extrap');

% --------------------------
% Smooth position (moving average)
% --------------------------
% w = max(1, round(params.smooth_win_s * fs));
% xS = movmean(x, w, 'Endpoints','shrink');
% yS = movmean(y, w, 'Endpoints','shrink');


fltord = 60;
lowpasfrq = 30;
nyqfrq = fs ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]); %30 Hz low pass filter
xS = filtfilt(flt,1,x);
yS = filtfilt(flt,1,y);

% --------------------------
% Velocity (deg/s)
% --------------------------
dx = [0; diff(xS)];
dy = [0; diff(yS)];
v_pxps = sqrt(dx.^2 + dy.^2) * fs;     % px/s
v_dps  = v_pxps / pxPerDeg;           % deg/s

% --------------------------
% I-VT labeling: saccade if v > threshold
% --------------------------
isSac = v_dps > params.vth_dps;

% Break at invalid original samples (common: don’t bridge blinks)
isSac(~valid) = false;

% Convert boolean runs to events
sacRuns = runsFromMask(isSac);
fixRuns = runsFromMask(~isSac & valid);

% --------------------------
% Filter by minimum duration
% --------------------------
minFixN = round(params.minFix_s * fs);
minSacN = round(params.minSac_s * fs);

fixRuns = fixRuns((fixRuns(:,2)-fixRuns(:,1)+1) >= minFixN, :);
sacRuns = sacRuns((sacRuns(:,2)-sacRuns(:,1)+1) >= minSacN, :);

% --------------------------
% Build fixation table
% --------------------------
nF = size(fixRuns,1);
fixStart = zeros(nF,1); fixEnd = zeros(nF,1); fixDur = zeros(nF,1);
fixX = zeros(nF,1); fixY = zeros(nF,1);
for i = 1:nF
    i0 = fixRuns(i,1); i1 = fixRuns(i,2);
    fixStart(i) = t(i0);
    fixEnd(i)   = t(i1);
    fixDur(i)   = fixEnd(i) - fixStart(i);
    fixX(i)     = mean(xS(i0:i1), 'omitnan');
    fixY(i)     = mean(yS(i0:i1), 'omitnan');
end
events.fixations = table(fixStart, fixEnd, fixDur, fixX, fixY, ...
    fixRuns(:,1), fixRuns(:,2), ...
    'VariableNames', {'start_t','end_t','dur_s','x_mean','y_mean','idx0','idx1'});

% --------------------------
% Build saccade table
% --------------------------
nS = size(sacRuns,1);
sacStart = zeros(nS,1); sacEnd = zeros(nS,1); sacDur = zeros(nS,1);
sacAmp = zeros(nS,1); sacPeakV = zeros(nS,1);
for i = 1:nS
    i0 = sacRuns(i,1); i1 = sacRuns(i,2);
    sacStart(i) = t(i0);
    sacEnd(i)   = t(i1);
    sacDur(i)   = sacEnd(i) - sacStart(i);

    % amplitude in degrees (start->end)
    amp_px = hypot(xS(i1)-xS(i0), yS(i1)-yS(i0));
    sacAmp(i) = amp_px / pxPerDeg;

    sacPeakV(i) = max(v_dps(i0:i1), [], 'omitnan');
end
events.saccades = table(sacStart, sacEnd, sacDur, sacAmp, sacPeakV, ...
    sacRuns(:,1), sacRuns(:,2), ...
    'VariableNames', {'start_t','end_t','dur_s','amp_deg','peakVel_dps','idx0','idx1'});

% Also return samplewise signals if useful
events.v_dps = v_dps;
events.isSac = isSac;
events.xS = xS; events.yS = yS;
events.valid = valid;

end

% -------- helper: contiguous runs from a logical mask --------
function runs = runsFromMask(mask)
mask = mask(:)~=0;
d = diff([false; mask; false]);
starts = find(d==1);
ends   = find(d==-1)-1;
runs = [starts ends];
end
