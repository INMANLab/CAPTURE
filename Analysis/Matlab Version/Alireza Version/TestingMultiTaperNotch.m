% --- Filter Test (Fs = 250 Hz, duration = 2 s) ---
clear; close all; clc;
rng(7);                                % reproducible noise
Fs   = 250;                            % sampling rate (Hz)
T    = 2;                              % duration (s)
t    = (0:1/Fs:T-1/Fs)';               % time vector (N=500)
FCut = 60:63;62.5;                           % target notch frequency (Hz)

% --- Build test signals ---
% 1) Pure tone at the notch frequency
sig1 = 1.0 * sin(2*pi*FCut*t);

% 2) Multi-tone + noise with a component at the notch frequency
sig2 = 0.8*sin(2*pi*10*t) + 0.6*sin(2*pi*30*t) + 0.5*sin(2*pi*FCut*t) ...
     + 0.20*randn(size(t));

% 3) Control (no component at the notch) to test collateral damage
sig3 = 0.8*sin(2*pi*10*t) + 0.6*sin(2*pi*30*t) + 0.20*randn(size(t));

signals = {sig1, sig2, sig3};
names   = {'Pure 62.5 Hz', '10+30+62.5 Hz + noise', '10+30 Hz + noise (no 62.5)'};

% --- Filter each signal ---
filtered = cellfun(@(x) Filter.MultiTaperNotchFilter(x, FCut, Fs, 'n'), signals, 'UniformOutput', false);

% --- Visualization helpers ---
Nfft = 512;  % finer frequency grid for nicer PSD curves
psdFun = @(x) periodogram(x, hamming(length(x)), Nfft, Fs, 'power');
getIdx = @(fvec,f0) find(abs(fvec - f0) == min(abs(fvec - f0)),1,'first');

% --- Figures: time & PSD before/after ---
for k = 1:numel(signals)
    x = signals{k};
    y = filtered{k};
    [Pxx, f] = psdFun(x);
    [Pyy, ~] = psdFun(y);

    figure('Name', ['Signal ' num2str(k) ' - ' names{k}], 'Color', 'w');
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    % Time (first 0.5 s for readability)
    nexttile; 
    plot(t, x, 'LineWidth', 1); xlim([0 0.5]); grid on;
    title(['Time (raw) — ' names{k}]); xlabel('Time (s)'); ylabel('Amp');

    nexttile; 
    plot(t, y, 'LineWidth', 1); xlim([0 0.5]); grid on;
    title('Time (filtered)'); xlabel('Time (s)'); ylabel('Amp');

    % PSD
    nexttile; 
    plot(f, 10*log10(Pxx+eps), 'LineWidth', 1); grid on; xlim([0 120]);
    xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); title('PSD (raw)');
    xline(FCut,'--','Notch f_0','LabelOrientation','horizontal');

    nexttile; 
    plot(f, 10*log10(Pyy+eps), 'LineWidth', 1); grid on; xlim([0 120]);
    xlabel('Frequency (Hz)'); ylabel('PSD (dB)'); title('PSD (filtered)');
    xline(FCut,'--','Notch f_0','LabelOrientation','horizontal');
end

% --- Quantitative sanity check: attenuation at FCut ---
fprintf('--- Attenuation at %.1f Hz (dB) ---\n', FCut);
for k = 1:numel(signals)
    x = signals{k}; y = filtered{k};
    [Pxx,f] = psdFun(x); [Pyy,~] = psdFun(y);
    i0 = getIdx(f, FCut);
    att_db = 10*log10((Pxx(i0)+eps)/(Pyy(i0)+eps));
    fprintf('Signal %d (%s): %.2f dB\n', k, names{k}, att_db);
end

%% frequency-sweep check around the notch ---
FCut = 62.5;
freqsToProbe = FCut + (-4:0.01:4);   % ±4 Hz around notch in 0.5 Hz steps
FCut = 60:.5:63;
amps = zeros(size(freqsToProbe));
for i = 1:numel(freqsToProbe)
    f0 = freqsToProbe(i);
    xt = sin(2*pi*f0*t);
    yt = Filter.MultiTaperNotchFilter(xt, FCut, Fs, 'n');
    [Pyy,f] = psdFun(yt);
    amps(i) = 10*log10(Pyy(getIdx(f, f0)) + eps);
end
figure('Color','w'); plot(freqsToProbe, amps,'-o','LineWidth',1); grid on;
xlabel('Injected tone (Hz)'); ylabel('Output PSD at tone (dB)');
title(sprintf('Notch selectivity around f_0 = %.1f Hz', FCut));
