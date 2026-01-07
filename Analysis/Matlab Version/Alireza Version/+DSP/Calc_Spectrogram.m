function [sp,fp,tp]= Calc_Spectrogram(d, fRange, Fs, winSize, nOverlap, method,parM)
% winSize in seconds
% nOverlap in seconds
winSize = round(Fs*winSize);
nOverlap = round(nOverlap*Fs);
tp = (1:length(d))/Fs;

switch method
    case "STFT"
        win = kaiser(winSize,5);
        [~,fp,tp,sp] = spectrogram(d,win,nOverlap,fRange,Fs,"power","yaxis");
    case "Wavelet_Morse"
        [fp,~,cfs] = morseSpecGram(d',Fs,[fRange(1),fRange(end)],parM);
        sp = abs(cfs).^2;
        sp = sp';
    case "cwt"
        [wtResults,fp] = cwt(d,Fs,FrequencyLimits=[fRange(1),fRange(end)]);
        sp = abs(wtResults).^2;
    case "morlet"
        sp = abs(DSP.calcWavTF(d,fRange,Fs,parM)).^2;
        fp = fRange;
    case "superlet"
        c1 = parM.c1;
        order = parM.order;
        mult = parM.mult;
        sp = faslt(d, Fs, fRange, c1, order, mult);
        % sp = abs(wtResults).^2;
        fp = fRange;
end
% [ps, fs] = pwelch(d, winSize, nOverlap, fRange, Fs); % Welch Power
% periodogram(xn, hamming(length(xn)), length(xn), fs, 'power');
% J = calcWavTF(D,F,Fs); %Morlet Wavelet
end

