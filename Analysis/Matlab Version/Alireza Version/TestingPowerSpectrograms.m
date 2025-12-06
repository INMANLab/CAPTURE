%%
Fs = 250;
t = 0:1/Fs:10;

winSize = 512;round(2*Fs);
nfft = 2^nextpow2(winSize);
g = kaiser(winSize,5);
ol = length(g)-1;
fRange = 0.5:0.5:85;

x1 = 3*sin(2*pi*5*t)+2*sin(2*pi*15*t)+sin(2*pi*50*t);
x2 = vco(sin(2*pi*t).*exp(-t),[10 80],Fs) ...
    + 0.1*sin(2*pi*20*t);

% figure
% plot(t,x2)
% figure
% spectrogram(x2,g,ol,fRange,Fs,"power","yaxis");

x3 = real(exp(1j*pi*sin(4*t)*40));
x4 = chirp(t,10,t(end),80,"linear");
x5 = x2+x3+x4;

%%
figure
subplot 411
plot(t,x1)
xlabel("time")
ylabel("Amplitude")
subplot 412
plot(t,x2)
xlabel("time")
ylabel("Amplitude")
subplot 413
plot(t,x3)
xlabel("time")
ylabel("Amplitude")
subplot 414
plot(t,x4)
xlabel("time")
ylabel("Amplitude")

%%
figure
subplot 411
[sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(tp,fp,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 412
[sp,fp,tp] = spectrogram(x2,g,ol,fRange,Fs,"power","yaxis");
mesh(tp,fp,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 413
[sp,fp,tp] = spectrogram(x3,g,ol,fRange,Fs,"power","yaxis");
mesh(tp,fp,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 414
[sp,fp,tp] = spectrogram(x4,g,ol,fRange,Fs,"power","yaxis");
mesh(tp,fp,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

%%

%%
figure
subplot 211
plot(t,x5)
xlabel("time")
ylabel("Amplitude")

subplot 212
[sp,fp,tp] = spectrogram(x5,g,ol,fRange,Fs,"power","yaxis");
mesh(tp,fp,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar


%% WaveLet
%%
figure
subplot 411
sp = abs(DSP.calcWavTF(x1,fRange,Fs)).^2;
% [sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(t,fRange,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 412
sp = abs(DSP.calcWavTF(x2,fRange,Fs)).^2;
% [sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(t,fRange,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 413
sp = abs(DSP.calcWavTF(x3,fRange,Fs)).^2;
% [sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(t,fRange,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

subplot 414
sp = abs(DSP.calcWavTF(x4,fRange,Fs)).^2;
% [sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(t,fRange,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar


%%
figure
subplot 211
plot(t,x5)
xlabel("time")
ylabel("Amplitude")

subplot 212
sp = abs(DSP.calcWavTF(x5,fRange,Fs)).^2;
% [sp,fp,tp] = spectrogram(x1,g,ol,fRange,Fs,"power","yaxis");
mesh(t,fRange,10*log10(abs(sp)))
title("spectrogram")
view(2), axis tight
xlabel("time")
ylabel("Frequency (Hz)")
colorbar

