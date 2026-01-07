clc
clear
close all;
%%
f      = 2;        % sine frequency (Hz)
fs     = 250;      % sampling rate (Hz)
tmax   = 2;       % t in [-t, t] (seconds)
sigma = .25;  % Gaussian std (seconds)
waveN = 2*pi*f*sigma;
t = (-tmax : 1/fs : tmax).';   % column vector, length N
%%

WaveletName = "Standard Gaussian";
[morletW,sineW,gaussianW] = GenerateMorlet(f,t,sigma);
figure
subplot 211
Plot_Signals(t,sineW,gaussianW)
title(WaveletName)
subplot 212
Plot_Morlet(t,morletW)

WaveletName = "CalcWavMorlet";
[morletW,sineW,gaussianW] = CalcWavMorlet(f,t,waveN);
figure
subplot 211
Plot_Signals(t,sineW,gaussianW)
title(WaveletName)
subplot 212
Plot_Morlet(t,morletW)

WaveletName = "Standard Morlet";
[morletW,sineW,gaussianW] = WaveMorlet(f,t,waveN);
figure
subplot 211
Plot_Signals(t,sineW,gaussianW)
title(WaveletName)
subplot 212
Plot_Morlet(t,morletW)
%% Functions
function [morletW,sineW,gaussianW] = GenerateMorlet(f,t,sigma)
    % ---- Parameters ----
    % f      % sine frequency (Hz)
    % t   % column vector, length N
    % sigma  % Gaussian variance (seconds)
    sineW = exp(1j*2*pi*f*t);                % complex sine wave
    % g = normpdf(t,0,sigma); % Alternative function from Matlab
    normalizationFactor = 1/sqrt(2*pi*sigma^2);
    gaussianW = normalizationFactor*exp(-(t.^2) ./ (2*sigma^2));  % Gaussian centered at 0 with variance sigma2
    morletW = sineW .* gaussianW;                             % windowed complex sine
    % abs(trapz(t,morletW))

end

function [morletW,sineW,gaussianW] = CalcWavMorlet(f,t,waveN)
    sineW = exp(1i*2*pi*f.*t);
    st = 1./(2*pi*(f/waveN));
    normalizationFactor = 1./sqrt(st*sqrt(pi));
    gaussianW = normalizationFactor * exp(-t.^2/(2*st^2));
    morletW = sineW.*gaussianW; % Morlet wavelet
end

function [morletW,sineW,gaussianW] = WaveMorlet(f,t,waveN)
    sineW = exp(1j*2*pi*f*t);                % complex sine wave
    sigma = 1/(2*pi*(f/waveN));
    normalizationFactor = 1/sqrt(2*pi*sigma^2);
    gaussianW = normalizationFactor * exp(-t.^2/(2*sigma^2));
    morletW = sineW.*gaussianW; % Morlet wavelet
end

% function [morletW,sineW,gaussianW] = WaveMorlet(f,t,waveN)
%     sineW = exp(1j*2*pi*f*t);                % complex sine wave
%     sigma = 1/(2*pi*(F/waveN));
%     normalizationFactor = 1/sqrt(2*pi*sigma^2);
%     gaussianW = normalizationFactor * exp(-t.^2/(2*st^2));
%     morletW = sineW.*gaussianW; % Morlet wavelet
% end

% ---- Plot magnitude ----
function Plot_Signals(t,sineW,gaussianW)
    hold on
    plot(t,real(sineW))
    plot(t,gaussianW)
    xlabel('Time (s)');
    title("Signals")
    legend("Sine  E = "+round(trapz(t,real(sineW)),2),"Gaussian  E = "+round(trapz(t,gaussianW),2))
end

function Plot_Morlet(t,y)
    plot(t, real(y), 'LineWidth', 1.5);
    hold on
    plot(t, imag(y), 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Real part of Morlet');
    title("Energy of the signal = " + round(trapz(t,y),2));
    legend("Real part E = "+round(trapz(t,real(y)),2),"Imaginary part E = "+round(trapz(t,real(y)),2))
end

