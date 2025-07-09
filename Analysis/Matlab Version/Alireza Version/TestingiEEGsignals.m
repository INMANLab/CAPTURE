%run from the Rasterplot GUI
dat = app.RWdat.allVars;
datEEG = dat.d_np;
% t = dat.ntp_np;
fs = 250;

t = (0:size(datEEG,1)-1)/fs;

winSize = 1; % in seconds
overlap = .75;

%%
ch = 1;
d = datEEG(:,ch);
d(isnan(d)) = 0;

[n,fo,mo,w] = firpmord([.5 1 50 50+.5], [0 1 0], [.01 .01 .01], fs);
b = firpm(n,fo,mo,w);
a=1;
%------------- Apply the filter
df = filtfilt(b,a, d);
d = df;
%%
figure
subplot 311
hold on

plot(t,d);

xVals = tArray/1000*60;
xline(xVals(evTimeSeries(:,1)==1),'r')
% xline(xVals(evTimeSeries(:,2)==1),'k')
legend("iEEG_CH"+ch,evUnique{1},evUnique{2})
xlabel("time(sec)")
axis tight
subplot 312
spectrogram(d,winSize*fs,round(winSize*fs*overlap),1:50,fs,'yaxis');
% ylim([0,10])
colorbar off
axis tight
[Ss,Ff,Tt] = spectrogram(d,winSize*fs,round(winSize*fs*overlap),[],fs,'yaxis');
ftoUse = Ff>4 & Ff<8;
thetaS = Ss(ftoUse,:);
thetaS = sum(thetaS);
subplot 313
plot(Tt,db(thetaS))
axis tight

%%
for tEventIdx = [4,7,11,13,15]
    tEvent = xVals(evTimeSeries(:,1)==1);
    
    figure
    subplot 311
    hold on
    idx = t>=(tEvent(tEventIdx)-10) & t<=(tEvent(tEventIdx)+10);
    plot(t(idx),d(idx));
    xline(tEvent(tEventIdx),'r')
    legend("iEEG_CH"+ch,evUnique{1},evUnique{2})
    xlabel("time(sec)")
    axis tight
    subplot 312
    spectrogram(d(idx),winSize*fs,round(winSize*fs*overlap),1:50,fs,'yaxis');
    [Ss,Ff,Tt] = spectrogram(d(idx),winSize*fs,round(winSize*fs*overlap),[],fs,'yaxis');
    ftoUse = Ff>4 & Ff<8;
    thetaS = Ss(ftoUse,:);
    thetaS = sum(thetaS);
    colorbar off
    axis tight
    subplot 313
    plot(Tt,db(thetaS))
    axis tight
end

