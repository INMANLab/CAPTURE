files = dir("*.dat");

datAll = [];
for fIdx = 1:length(files)
    fid = fopen(files(fIdx).name,'r');
    dat = fread(fid,'int16');
    fclose(fid);
    dat = reshape(dat,4,[])'-512;
    datAll = cat(1,datAll,dat);
end
%%
A = datAll(215600:291700,1);
A(A>-500)=0;
plot(A)


B = A(6349:6364);
figure
subplot 211
plot(B)
xlabel("sample")
subplot 212
Y = fftshift(fft(B,64));
L = length(Y);
Fs = 250;
plot(Fs/L*(-L/2:L/2-1),abs(Y),"LineWidth",3);
hold on
xline(62.5)
xline(-62.5)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")

B = A(6349:65799);
figure
subplot 211
plot(B)
xlabel("sample")
subplot 212
Y = fftshift(fft(B));
L = length(B);
Fs = 250;
plot(Fs/L*(-L/2:L/2-1),abs(Y),"LineWidth",3);
hold on
xline(62.5)
xline(-62.5)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")


B = A(6349:6364);
plot(B)
% B(10:11)=-512;
figure
subplot 211
stem(B)
xlabel("sample")
subplot 212
Y = fftshift(fft(B,32));
L = length(Y);
Fs = 250;
plot(Fs/L*(-L/2:L/2-1),abs(Y),"LineWidth",3);
hold on
xline(62.5)
xline(-62.5)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")
%%
figure
chIdx = 1;
gDat = npDat(:,chIdx);
plot(t,gDat)
hold on
figure
idx = (0:size(datAll,1)-1)/250;
plot(idx,datAll(:,chIdx))
L = size(datAll,1);
plot(t,datAll(L-length(t)+1:end,chIdx))
%%
% 
figure
for chIdx = 1:4
    idx = 1:size(datAll,1);
    % chunkIdx = [1,idx(datAll(:,chIdx)<=-500)];
    t = [1,idx(datAll(:,chIdx)<=-500)];
    tE = t([1,diff(t)]>1)/(250*60);
    tS = t([diff(t),1]>1)/(250*60);
    % dt = diff(t);
    % dt = dt(dt>1);
    disp("Recording Duration")
    disp(tE-tS);
    disp("Reset Duration")
    disp(tE(1:end-1)-tS(2:end))
    subplot(4,1,chIdx)
    plot(idx/(250*60),datAll(:,chIdx))
    hold on
    plot(t/(250*60),datAll(t),'*')
    xline(t([diff(t),1]>1)/(250*60))
    xline(t([1,diff(t)]>1)/(250*60))
end

%%
figure
for chIdx = 1:4
    idx = 1:size(datAll,1);
    datIdx = datAll(:,chIdx);
    datIdx(datIdx>-500)=0;
    tIdx = reshape([1;find(diff(datIdx));length(datIdx)],2,[]);
    % tE = t([1,diff(t)]>1)/(250*60);
    % tS = t([diff(t),1]>1)/(250*60);
    % dt = diff(t);
    % dt = dt(dt>1);
    disp("==========Recording Durations (Minutes):")
    disp(diff(tIdx)*-1/(250*60));
    disp("----------Reset Duration (Seconds)")
    disp(tIdx(1,2:end)-tIdx(2,1:end-1))

    subplot(4,1,chIdx)
    plot(idx/(250*60),datAll(:,chIdx))
    hold on
    plot(idx(datIdx~=0)/(250*60),datAll(datIdx~=0,chIdx),'*')
    xline(tIdx(1,:)/(250*60),'r')
    xline(tIdx(2,:)/(250*60),'b')
end



%%
tEnd = 10;
t = 1/250:1/250:tEnd;
t = t*1000;
tN = 1/62.5:1/62.5:tEnd;
tN = tN *1000;

figure
subplot 211
plot(t,dat(1:length(t),1))
xlabel("time (ms)")
hold on 
xline(tN)
subplot 212
pwelch(dat(1:length(t)),[],[],[],250)

wo = 62.5/(250/2);  
bw = wo/35;
[b62,a62] = iirnotch(wo, bw); % 60Hz IIR filter
datF = filtfilt(b62,a62,dat);

figure
subplot 211
plot(t,datF(1:length(t),1))
xlabel("time (ms)")
hold on 
xline(tN)
subplot 212
pwelch(datF(1:length(t)),[],[],[],250)

figure
plot(t,dat(1:length(t),1))
hold on
plot(t,datF(1:length(t),1))

%%
[n,fo,mo,w] = firpmord([50-1 50 120 120+1], [0 1 0], [1 1 1]*0.01, 250);
b = firpm(n,fo,mo,w);
a=1;
%------------- Apply the filter
datF2 = filtfilt(b,a, dat);
figure
plot(t,datF2(1:length(t),1))
hold on
plot(t,sin(2*pi*62.5*t/1000)*10)


%% 
y = sin(2*pi*10*t/1000)+sin(2*pi*15*t/1000);
y(64:64:length(t))= y(64:64:length(t))+2;
% plot(t,y)
pwelch(y,[],[],[],250)