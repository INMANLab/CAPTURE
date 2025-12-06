function Calc_Spectrogram(d, fRange, Fs, method, winSize, nOverlap)

[ps, fs] = pwelch(d, winSize, nOverlap, fRange, Fs); % Welch Power

% periodogram(xn, hamming(length(xn)), length(xn), fs, 'power');
[s,f,t] = spectrogram(d,winSize,noverlap,fRange,fs,"psd","yaxis","reassigned"); % STFT


J = calcWavTF(D,F,Fs); %Morlet Wavelet

end

