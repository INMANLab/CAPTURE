%%
clear;
clc;
% close all;
%% Signal Values
Fs = 250;
t = 1/Fs:1/Fs:10;
effectRange = [2,4];

d.x1 = sin(2*pi*2*t) +... Delta 2 HZ
       sin(2*pi*6*t) +... Theta 6 HZ
       sin(2*pi*10.5*t)+... Alpha 10 Hz
       sin(2*pi*20*t)+... Beta 20 Hz
       sin(2*pi*50*t);    % Gamma 50Hz

effectRange = round(effectRange*Fs);
winEff = zeros(size(d.x1));
winEff(effectRange(1):effectRange(2)) = 1;
d.x1 = d.x1.*winEff;
d.x1 = d.x1+randn(size(d.x1))*.1;

% figure
% plot(t,d.x1)

d.x2 = vco(sin(2*pi*t).*exp(-t),[10 80],Fs) + 0.1*sin(2*pi*20*t);

% figure
% plot(t,x2)
% figure
% spectrogram(x2,g,ol,fRange,Fs,"power","yaxis");

d.x3 = real(exp(1j*pi*sin(4*t)*40));
d.x4 = chirp(t,2,t(end),80,"linear")+randn(size(d.x1))*.1;

d.x5 = d.x4;d.x1 + 2*d.x2 + 2*d.x3 + 2*d.x4;
%% Power Parameters
fRange = 1:0.5:85;
winSize = 2;
nOverlap = winSize-1/Fs;
timeBW =80; 

%% Power STFT
method = "STFT";
figure
for i=[1,5]
    data = d.("x"+i);
    [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method, timeBW);

    nexttile
    plot(t,data)
    xlabel("time")
    ylabel("Amplitude")
    title("Effect Fs="+Fs+" |Winsize="+winSize, Interpreter="none")

    nexttile
    mesh(tp,fp,10*log10(abs(sp)))
    hold on 
    yline(2)
    yline(6)
    yline(10)
    yline(20)
    yline(50)
    xline(effectRange(1)/Fs)
    xline(effectRange(2)/Fs)
    plot3(t,linspace(2,80,length(t)),repmat(max(10*log10(abs(sp)),[],'all'),1,length(t)),'k');
    xlim([tp(1),tp(end)])
    title("spectrogram " + method, Interpreter="none")
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    xlim([tp(1),tp(end)])
    colorbar
end

%% Power Wavelet_Morse
method = "Wavelet_Morse";
figure
for i=[1,5]
    data = d.("x"+i);
    [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method, timeBW);
    % [sp1,fp1,tp]= DSP.Calc_Spectrogram(data, fRange(1:60), Fs, winSize, nOverlap, method, 10);
    % [sp2,fp2,tp]= DSP.Calc_Spectrogram(data, fRange(61:100), Fs, winSize, nOverlap, method, 50);
    % [sp3,fp3,tp]= DSP.Calc_Spectrogram(data, fRange(101:end), Fs, winSize, nOverlap, method, 90);

    nexttile
    plot(t,data)
    xlabel("time")
    ylabel("Amplitude")
    title("Effect Fs="+Fs+" |Winsize="+winSize, Interpreter="none")

    nexttile
    mesh(tp,fp,10*log10(abs(sp)))
    hold on 
    xline(effectRange(1)/Fs)
    xline(effectRange(2)/Fs)
    plot3(t,linspace(2,80,length(t)),repmat(max(10*log10(abs(sp)),[],'all'),1,length(t)),'k');
    xlim([tp(1),tp(end)])
    title("spectrogram " + method, Interpreter="none")
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    xlim([tp(1),tp(end)])
    colorbar
end

%% Power morlet
method = "morlet";
figure
for i=[1,5]
    data = d.("x"+i);
    % [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method, 20);
    [sp1,fp1,~]= DSP.Calc_Spectrogram(data, fRange(fRange<4), Fs, winSize, nOverlap, method, 3);
    [sp2,fp2,~]= DSP.Calc_Spectrogram(data, fRange(fRange>=4 & fRange<8), Fs, winSize, nOverlap, method, 6);
    [sp3,fp3,~]= DSP.Calc_Spectrogram(data, fRange(fRange>=8 & fRange<12), Fs, winSize, nOverlap, method, 12);
    [sp4,fp4,~]= DSP.Calc_Spectrogram(data, fRange(fRange>=12 & fRange<30), Fs, winSize, nOverlap, method, 24);
    [sp5,fp5,tp]= DSP.Calc_Spectrogram(data, fRange(fRange>=30), Fs, winSize, nOverlap, method, 32);
    sp = cat(1,sp1,sp2,sp3,sp4,sp5);
    fp = cat(2,fp1,fp2,fp3,fp4,fp5);

    nexttile
    plot(t,data)
    xlabel("time")
    ylabel("Amplitude")
    title("Effect Fs="+Fs+" |Winsize="+winSize, Interpreter="none")

    nexttile
    mesh(tp,fp,10*log10(abs(sp)))
    hold on 
    xline(effectRange(1)/Fs)
    xline(effectRange(2)/Fs)
    plot3(t,linspace(2,80,length(t)),repmat(max(10*log10(abs(sp)),[],'all'),1,length(t)),'k');    
    title("spectrogram " + method, Interpreter="none")
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    xlim([tp(1),tp(end)])
    colorbar
end
%% Power superlet
method = "superlet";
figure
for i=[1,5]
    data = d.("x"+i);
    parM.c1 = 3;
    parM.order = [1,50];
    parM.mult = 1;
    % [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method, parM);
    [sp,fp,tp] = DSP.Spectrogram_Superlet(data', Fs, [fRange(1),fRange(end)]);
    sp = sp';
    nexttile
    plot(t,data)
    xlabel("time")
    ylabel("Amplitude")
    title("Effect Fs="+Fs+" |Winsize="+winSize, Interpreter="none")

    nexttile
    mesh(tp,fp,10*log10(abs(sp)))
    hold on 
    xline(effectRange(1)/Fs)
    xline(effectRange(2)/Fs)
    plot3(t,linspace(2,80,length(t)),repmat(max(10*log10(abs(sp)),[],'all'),1,length(t)),'k');    
    title("spectrogram " + method, Interpreter="none")
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    xlim([tp(1),tp(end)])
    colorbar
end
%% Time Values
figure
for i=1:5
    subplot(5,1,i)
    plot(t,d.("x"+i))
    xlabel("time")
    ylabel("Amplitude")
end

%% Power STFT
method = "STFT";
figure
for i=1:5
    data = d.("x"+i);
    [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method);
    subplot(5,1,i)
    mesh(tp,fp,10*log10(abs(sp)))
    title("spectrogram " + method)
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    colorbar
end

%% Power Wavelet
method = "Wavelet_Morse";
figure
for i=1:5
    data = d.("x"+i);
    [sp,fp,tp]= DSP.Calc_Spectrogram(data, fRange, Fs, winSize, nOverlap, method);
    subplot(5,1,i)
    mesh(tp,fp,10*log10(abs(sp)))
    title("spectrogram " + method)
    view(2), axis tight
    xlabel("time")
    ylabel("Frequency (Hz)")
    colorbar
end

%% STFT

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

%% Wavelet

figure
subplot 411
[fp,~,cfs] = morseSpecGram(x4',Fs,[fRange(1),fRange(end)]);
sp = abs(cfs).^2;
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

