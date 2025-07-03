% Foof plot test

%% FOOOF Matlab Wrapper Example - Plot a FOOOF model

% Load data
load('npData.mat');


s_rate = 250;
% Calculate a power spectrum with Welch's method
[psd, freqs] = pwelch(np(tidx_box1_np,1), 500, [], [], s_rate);

[psd, freqs] = pwelch(np(tidx_box2_np,1), 500, [], [], s_rate);

% Transpose, to make inputs row vectors
freqs = freqs';
psd = psd';

% FOOOF settings
% settings = struct();  % Use defaults
f_range = [1, 85];

% Run FOOOF, also returning the model
fooof_results = fooof(freqs, psd, f_range, settings, true);

% Plot the resulting model
fooof_plot(fooof_results)