clear;
clc;

GUIhandle = RWAnalysis2;
% GUIhandle.AnalysisFile = 'D:\CAPTURE Project\GUI\Data\RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;


% h = GUIhandle.openMultFRankingGUI;
% 
% GUIhandle.openMultSegGLMEGUI;

MultTransSpecGramGUI_AK(GUIhandle);
GUIhandle.openMultTransSpecGramGUI;

%% load data
T = readtable("D:\CAPTURE Project\GUI\Data\MultSegData.csv");
chanMap = unique(T(:,["Pt","Chan","RegionLabel","ChanLabel"]));
region = "AntHipp";


files = dir('*.mat');
fileNames = {files.name};
tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
if isempty(tokens) % Check if patient files were found
    disp('No patient file has been found');
    folderPath = uigetdir(pwd, 'Select the folder that contains individual data');
    if folderPath ~= 0 % If the user selects a folder
        files = dir(fullfile(folderPath,'*.mat'));
        fileNames = {files.name};
        tokens = TextProcess.ExtractPatientWalkFromFileNames(fileNames, 'RW(\d+)_Walk(\d+)\.mat');
        addpath(genpath(folderPath));
    else
        error('Patient files not found');
    end
end
Patients = tokens(:,1);
Walks = tokens(:,2);


dat = cell(1,length(Patients));
for idx = 1:length(Patients)
    d = load("RWNApp_RW"+Patients(idx)+"_Walk"+Walks(idx)+".mat");
    dat{idx} = d;
end

%%
pID = 1;
wID = 5;
chanIdx = strcmp(chanMap.RegionLabel,region);
chanIdx = chanIdx(chanMap.Pt==pID);
idx = find(Patients==pID & Walks == wID);

d = dat{idx};
evTable = d.evnts_tbl;
startNTP = evTable.NTP(evTable.Event=="Walk Beg");
endNTP = evTable.NTP(evTable.Event=="Walk End");
% startNTP = evTable.NTP(evTable.Event=="Lost Beg");
% endNTP = evTable.NTP(evTable.Event=="Lost End");

% startNTP = startNTP(1)-10;
% endNTP = startNTP+20;

npDat = d.d_np(:,chanIdx);
fsNP = d.fs_np;
tNTP = d.ntp_np;
[~,startIdx] = min(abs(tNTP-startNTP(1)));
[~,endIdx] = min(abs(tNTP-endNTP(1)));
npDat = npDat(startIdx:endIdx,:);
npDat(isnan(npDat))=0;
npDat(npDat<-400)=NaN;

fxDat = d.d_gaze_fix;
fsGaze = d.fs_gaze;
tNTP = d.ntp_gaze;
[~,startIdx] = min(abs(tNTP-startNTP(1)));
[~,endIdx] = min(abs(tNTP-endNTP(1)));
fxDat = fxDat(startIdx:endIdx,:);
fxDat = diff([0;fxDat]);
fxDat(fxDat<0)=0;

% -------------------- USER SETTINGS --------------------
bandName = 'theta';  % 'delta','theta','alpha','beta','gamma','custom'
customBand = [1 40]; % used only if bandName='custom' (Hz)

doZscore = true;     % zscore fixation and EEG before correlation
corrType = 'Pearson';% 'Pearson' or 'Spearman'
maxLagSec = 0;       % 0 -> no lag (simple corr). If >0, compute xcorr within +/- maxLagSec.

% -------------------- SELECT BAND ----------------------
switch lower(bandName)
    case 'delta', bandHz = [1 4];
    case 'theta', bandHz = [4 8];
    case 'alpha', bandHz = [8 12];
    case 'beta',  bandHz = [13 30];
    case 'gamma', bandHz = [30 80];   % adjust as needed
    case 'custom', bandHz = customBand;
    otherwise
        error('Unknown bandName.');
end

% -------------------- INPUTS ---------------------------
fixGaze = fxDat;  % force column
eegRaw  = npDat;                     % [Nnp x nCh]

Nnp = size(eegRaw,1);
nCh = size(eegRaw,2);

% -------------------- RESAMPLE FIXATION TO fsNP --------
% Build time vectors and resample fixation onto EEG sampling grid.
tFix = (0:numel(fixGaze)-1)'/fsGaze;
tEEG = (0:Nnp-1)'/fsNP;

% Convert to double (interp1 needs numeric); keep as 0/1
fixEEG = interp1(tFix, double(fixGaze), tEEG, 'previous', 'extrap');

fixEEG = diff([0;fixEEG]);
fixEEG(fixEEG<0)=0;

% Optional: smooth a little if you want less “step-like” regressor
% fixEEG = movmean(fixEEG, round(0.050*fsNP)); % 50 ms window   


%
idx = find(fixEEG);
idx = idx(idx>round((6*fsNP)));
% idx = idx(idx<idx(end)+round((6*fsNP)));
EEGMat = zeros(length(idx(1)-round((2*fsNP)):idx(1)+round((2*fsNP))),length(idx));

for i=1:(length(idx)-1)
    EEGMat(:,i) = eegRaw(idx(i)-round((2*fsNP)):idx(i)+round((2*fsNP)),1);
end
size(EEGMat,1)
figure
plot((1:1001)/fsNP,mean(EEGMat,2,"omitnan"))
xline(2)
%%
% -------------------- BANDPASS FILTER EEG --------------

% Zero-phase IIR bandpass (stable, easy). For more control use FIR.
d = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', bandHz(1), ...
    'HalfPowerFrequency2', bandHz(2), ...
    'SampleRate', fsNP);

[N, nCh] = size(eegRaw);

eegFilt = nan(N, nCh);                  % keep NaNs where data are missing (optional)
badMaskAll = ~isfinite(eegRaw);         % NaN/Inf locations

for ch = 1:nCh
    x = eegRaw(:, ch);
    bad = ~isfinite(x);

    if all(bad)
        warning('Channel %d is entirely NaN/Inf. Skipping.', ch);
        continue;
    end

    % Fill NaN/Inf so filtfilt can run
    % - linear interpolation inside
    % - nearest at the edges (if NaNs at start/end)
    xFilled = fillmissing(x, 'linear', 'EndValues', 'nearest');

    % Safety: if any Inf still present (rare), treat as missing too
    xFilled(~isfinite(xFilled)) = 0;

    % Zero-phase filtering
    y = filtfilt(d, xFilled);

    % Option A (recommended): restore missing samples to NaN so downstream analyses know
    y(bad) = NaN;

    eegFilt(:, ch) = y;
end



% Inputs:
%   fixEEG  : [N x 1] fixation samples on EEG grid (0/1 or numeric)
%   eegFilt : [N x nCh] filtered EEG (may include NaNs)
%   corrType: 'Pearson' or 'Spearman'

xAll = double(fixEEG(:));
Y    = double(eegFilt);

[N, nCh] = size(Y);

r = nan(nCh,1);
p = nan(nCh,1);
nUsed = zeros(nCh,1);

for ch = 1:nCh
    yAll = Y(:,ch);

    % keep only finite pairs
    good = isfinite(xAll) & isfinite(yAll);
    x = xAll(good);
    y = yAll(good);
    nUsed(ch) = numel(x);

    % Need enough samples
    if nUsed(ch) < 10
        continue
    end

    % If either vector has ~zero variance, corr is undefined -> NaN
    if std(x) < 1e-12 || std(y) < 1e-12
        continue
    end

    % Correlation
    if strcmpi(corrType,'Spearman')
        [r(ch), p(ch)] = corr(x, y, 'Type','Spearman');
    else
        [r(ch), p(ch)] = corr(x, y, 'Type','Pearson');
    end
end

results = table((1:nCh)', nUsed, r, p, ...
    'VariableNames', {'Channel','NSamplesUsed','r','pValue'});
disp(results);


%
% Spectrogram + phase plots per channel (NaN/Inf-safe),
% with fixation dots overlaid.
%
% Inputs expected in workspace:
%   npDat   : [N x nCh] raw EEG
%   fsNP    : EEG sampling rate (Hz)
%   fixEEG  : [N x 1] fixation samples on EEG time base (0/1)
%   bandHz  : [f1 f2] EEG band (Hz), e.g., [8 12]
%
% Notes on NaNs:
%   - We FILL NaNs/Inf just for the spectrogram computation (required).
%   - We also compute a "badness fraction" per TF time-bin and mask those
%     time bins (set to NaN) if too much missing data was present.

eegRaw = double(npDat);
fixEEG = double(fixEEG);
[N, nCh] = size(eegRaw);
tEEG = (0:N-1)'/fsNP;

% ----- Spectrogram parameters -----
winSec   = 1.0;         % seconds
ovlpFrac = 0.80;        % overlap fraction
winSamp  = max(16, round(winSec*fsNP));
win      = hamming(winSamp, 'periodic');
noverlap = round(ovlpFrac*winSamp);
nfft     = 2^nextpow2(winSamp);

% ----- Band + center frequency -----
f0 = mean(bandHz);

% ----- Masking rule: if a TF time-bin has too many missing samples, hide it -----
maxMissingFracPerBin = 0.20;  % e.g., if >20% of samples in window were NaN/Inf, mask that time-bin

for ch = 1:nCh
    x = eegRaw(:, ch);

    % Identify invalid samples
    bad = ~isfinite(x);

    if all(bad)
        warning('Channel %d is entirely NaN/Inf. Skipping.', ch);
        continue;
    end

    % Fill for spectrogram computation (spectrogram requires finite values)
    xFilled = fillmissing(x, 'linear', 'EndValues', 'nearest');
    xFilled(~isfinite(xFilled)) = 0; % extra safety

    % ---- Compute complex spectrogram ----
    [S, F, T] = spectrogram(xFilled, win, noverlap, nfft, fsNP);

    % ---- Compute how "bad" each TF time bin was (based on original bad mask) ----
    % spectrogram time centers correspond to window centers:
    % window start indices are 1 + k*(winSamp-noverlap), for k = 0..nFrames-1
    hop = winSamp - noverlap;
    nFrames = numel(T);
    startIdx = 1 + (0:nFrames-1)*hop;
    endIdx   = startIdx + winSamp - 1;

    % clamp to signal bounds (should already be in bounds, but be safe)
    endIdx(endIdx > N) = N;

    missFrac = zeros(1, nFrames);
    for k = 1:nFrames
        segBad = bad(startIdx(k):endIdx(k));
        missFrac(k) = mean(segBad);
    end

    timeBinBad = missFrac > maxMissingFracPerBin;

    % ---- Power + band selection ----
    P   = abs(S).^2;
    PdB = 10*log10(P + eps);

    fMask = (F >= bandHz(1)) & (F <= bandHz(2));
    Fband = F(fMask);
    PdBband = PdB(fMask, :);

    % Mask bad time bins in the spectrogram display
    PdBband(:, timeBinBad) = NaN;

    % ---- Phase at center frequency (from closest bin), mask bad bins ----
    [~, iF0] = min(abs(F - f0));
    phi = angle(S(iF0, :));
    phi(timeBinBad) = NaN;

    % Optionally unwrap (unwrap ignores NaNs poorly, so do it carefully)
    % If you want unwrapped phase:
    phiUnwrap = phi;
    goodPhi = isfinite(phiUnwrap);
    if any(goodPhi)
        phiUnwrap(goodPhi) = unwrap(phiUnwrap(goodPhi));
    end

    % ---- Fixation dots on spectrogram time grid ----
    % Map fixation to T; also avoid placing dots where the TF bin is masked
    fixAtT = interp1(tEEG, fixEEG, T, 'previous', 'extrap') > 0.5;
    fixAtT = fixAtT & ~timeBinBad;

    dotT = T(fixAtT);

    % ---- Plot ----
    figure('Name', sprintf('Channel %d', ch), 'Color', 'w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % (1) Power spectrogram in band
    nexttile;
    imagesc(T, Fband, PdBband);
    axis xy;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Ch %d: Power Spectrogram (%g-%g Hz) [NaN-safe]', ch, bandHz(1), bandHz(2)));
    colorbar;

    hold on;
    yDots1 = bandHz(2) * ones(size(dotT));
    plot(dotT, yDots1, '.', 'MarkerSize', 12);
    hold off;

    % (2) Phase at f0
    nexttile;
    plot(T, phiUnwrap, 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    title(sprintf('Ch %d: Phase at f0 = %.2f Hz (masked bins removed)', ch, f0));
    grid on;

    hold on;
    yl = ylim;
    yDots2 = (yl(2) - 0.05*(yl(2)-yl(1))) * ones(size(dotT));
    plot(dotT, yDots2, '.', 'MarkerSize', 12);
    plot(T,fixAtT,'*')
    hold off;
end

% -------------------- (OPTIONAL) QUICK PLOT ------------
% Plot fixation regressor and one example channel
exampleCh = 1;
figure;
plot(tEEG, x); hold on;
plot(tEEG, Y(:,exampleCh));
legend('Fixation (EEG grid)','Filtered EEG (example ch)');
xlabel('Time (s)');
title(sprintf('Band: %s [%g-%g] Hz, Ch %d', bandName, bandHz(1), bandHz(2), exampleCh));

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

pID = 5;
wID = 5;
chanIdx = strcmp(chanMap.RegionLabel,region);
chanIdx = chanIdx(chanMap.Pt==pID);
idx = find(Patients==pID & Walks == wID);

d = dat{idx};
evTable = d.evnts_tbl;
startNTP = evTable.NTP(evTable.Event=="Walk Beg");
endNTP = evTable.NTP(evTable.Event=="Walk End");

npDat = d.d_np(:,chanIdx);
fsNP = d.fs_np;
tNTP = d.ntp_np;
[~,startIdx] = min(abs(tNTP-startNTP));
[~,endIdx] = min(abs(tNTP-endNTP));
npDat = npDat(startIdx:endIdx,:);
npDat(isnan(npDat))=0;
npDat(npDat<-400)=NaN;

fxDat = d.d_gaze_fix;
fsGaze = d.fs_gaze;
tNTP = d.ntp_gaze;
[~,startIdx] = min(abs(tNTP-startNTP));
[~,endIdx] = min(abs(tNTP-endNTP));
fxDat = fxDat(startIdx:endIdx,:);

% Fixation–Phase Association (NaN-safe)
%
% Inputs expected:
%   fxDat.logicalFixation : [Nfix x 1] logical (0/1) at fsGaze
%   fsGaze                : fixation sampling rate (Hz)
%   npDat                 : [Neeg x nCh] raw EEG at fsNP
%   fsNP                  : EEG sampling rate (Hz)
%
% What this does (per channel):
%   1) Resample fixation vector onto EEG time base.
%   2) Bandpass filter EEG in a selected band (NaN-safe) and compute analytic phase via Hilbert.
%   3) Extract phase samples during fixation vs non-fixation.
%   4) Quantify fixation–phase association using:
%        - Mean Resultant Length (MRL) of phase during fixations (phase concentration)
%        - Phase-Locking Value (PLV) between a binary fixation regressor and unit phasor
%        - Rayleigh test p-value for non-uniformity of fixation-phase angles (approximate)
%      Plus a permutation p-value (recommended) that is robust and avoids assumptions.
%
% Notes:
%   - If you care about fixation ONSETS specifically (more interpretable), set useOnsets=true.

% -------------------- SETTINGS --------------------
bandName = 'theta';        % 'delta','theta','alpha','beta','gamma','custom'
customBand = [8 12];       % only if bandName='custom'
useOnsets = true;          % true = analyze fixation onset samples; false = all fixation samples
minFixSamples = 50;        % minimum number of fixation samples needed per channel
nPerm = 1000;              % permutation count for p-values (increase for publication)
maxShiftSec = 5;           % circular shift range for permutation (sec). Must be < recording length.

% pick band edges
switch lower(bandName)
    case 'delta', bandHz = [1 4];
    case 'theta', bandHz = [4 8];
    case 'alpha', bandHz = [8 12];
    case 'beta',  bandHz = [13 30];
    case 'gamma', bandHz = [30 80];
    case 'custom', bandHz = customBand;
    otherwise, error('Unknown bandName.');
end

% -------------------- RESAMPLE FIXATION TO EEG GRID --------------------
fixGaze = fxDat;
Neeg = size(npDat,1);
tFix = (0:numel(fixGaze)-1)'/fsGaze;
tEEG = (0:Neeg-1)'/fsNP;

% step-wise resample (keeps fixation as 0/1 blocks)
fixEEG = interp1(tFix, double(fixGaze), tEEG, 'previous', 'extrap') > 0.5;

% optional: use fixation onsets only
if useOnsets
    fixMask = fixEEG & ~[false; fixEEG(1:end-1)];
else
    fixMask = fixEEG;
end

% -------------------- FILTER DESIGN --------------------
d = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', bandHz(1), ...
    'HalfPowerFrequency2', bandHz(2), ...
    'SampleRate', fsNP);

% -------------------- HELPERS --------------------
circ_mean = @(ang) angle(mean(exp(1j*ang)));
circ_mrl  = @(ang) abs(mean(exp(1j*ang))); % mean resultant length in [0,1]

% Approximate Rayleigh p-value (good for moderate/large N; use permutation for robustness)
rayleigh_p = @(R,N) exp(sqrt(1 + 4*N + 4*(N^2 - (N*R)^2)) - (1 + 2*N));

rng(0); % reproducible permutations

% -------------------- MAIN LOOP --------------------
eegRaw = double(npDat);
[Neeg, nCh] = size(eegRaw);

% results containers
MRL_fix   = nan(nCh,1);
MRL_non   = nan(nCh,1);
PLV_fix   = nan(nCh,1);
meanPhase = nan(nCh,1);   % circular mean phase during fixations (radians)
pRay_fix  = nan(nCh,1);   % Rayleigh p-value (approx)
pPerm_MRL = nan(nCh,1);   % permutation p for MRL (fixation vs shifted mask)

NfixUsed  = zeros(nCh,1);
NnonUsed  = zeros(nCh,1);

% limit for shifts
maxShift = round(maxShiftSec * fsNP);
if maxShift*2 >= Neeg
    maxShift = floor(Neeg/4); % keep safe if very short signal
end

for ch = 1:nCh
    x = eegRaw(:,ch);

    % ----- NaN/Inf handling for filtering -----
    bad = ~isfinite(x);
    if all(bad)
        warning('Channel %d: all NaN/Inf. Skipping.', ch);
        continue;
    end

    % Fill missing just to filter / Hilbert (keep mask to exclude later)
    xFilled = fillmissing(x, 'linear', 'EndValues', 'nearest');
    xFilled(~isfinite(xFilled)) = 0;

    % Filter (zero-phase) and analytic phase
    y = filtfilt(d, xFilled);
    phi = angle(hilbert(y)); % instantaneous phase

    % Exclude originally-bad EEG samples and also exclude where fixMask is NaN (shouldn't be)
    goodEEG = ~bad & isfinite(phi);

    % fixation and non-fixation samples (on EEG grid)
    fixIdx = fixMask & goodEEG;
    nonIdx = ~fixMask & goodEEG;

    angFix = phi(fixIdx);
    angNon = phi(nonIdx);

    NfixUsed(ch) = numel(angFix);
    NnonUsed(ch) = numel(angNon);

    if NfixUsed(ch) < minFixSamples
        % too few fixation samples -> unstable circular stats
        continue;
    end

    % ----- Core circular stats -----
    MRL_fix(ch)   = circ_mrl(angFix);
    MRL_non(ch)   = circ_mrl(angNon);
    meanPhase(ch) = circ_mean(angFix);

    % "PLV" style coupling: fixation mask as weights for unit phasor
    % Equivalent to MRL over fixation phases if weights are 0/1 (so this is redundant but explicit)
    PLV_fix(ch) = abs(sum(exp(1j*phi).*double(fixIdx)) / sum(double(fixIdx)));

    % Rayleigh test (approx); permutation below is the robust one
    R = MRL_fix(ch);
    pRay_fix(ch) = rayleigh_p(R, NfixUsed(ch));

    % ----- Permutation test (circularly shift fixation mask) -----
    % Null: fixations occur at unrelated times -> phase samples at "fixations" look like random time samples.
    % We preserve fixation structure and EEG autocorrelation by circularly shifting fixMask.
    permMRL = nan(nPerm,1);
    for k = 1:nPerm
        shift = randi([-maxShift, maxShift], 1, 1);
        fixPerm = circshift(fixMask, shift) & goodEEG; % also respect valid EEG
        angPerm = phi(fixPerm);

        if numel(angPerm) < minFixSamples
            permMRL(k) = NaN;
        else
            permMRL(k) = circ_mrl(angPerm);
        end
    end

    permMRL = permMRL(isfinite(permMRL));
    if ~isempty(permMRL)
        % one-sided: are fixations MORE phase-concentrated than expected?
        pPerm_MRL(ch) = (sum(permMRL >= MRL_fix(ch)) + 1) / (numel(permMRL) + 1);
    end
end

results = table((1:nCh)', NfixUsed, NnonUsed, MRL_fix, MRL_non, PLV_fix, meanPhase, pRay_fix, pPerm_MRL, ...
    'VariableNames', {'Channel','Nfix','Nnonfix','MRL_fix','MRL_nonfix','PLV_fix','MeanPhase_rad','pRayleighApprox','pPerm_MRL'});
disp(results);

% -------------------- OPTIONAL: QUICK SUMMARY PLOT --------------------
figure('Color','w'); 
subplot(2,1,1);
stem(MRL_fix, 'filled'); ylabel('MRL (fix)'); xlabel('Channel'); title(sprintf('Phase concentration during %s (%g-%g Hz)', bandName, bandHz(1), bandHz(2)));
subplot(2,1,2);
stem(pPerm_MRL, 'filled'); ylabel('Permutation p'); xlabel('Channel'); ylim([0 1]); grid on;




