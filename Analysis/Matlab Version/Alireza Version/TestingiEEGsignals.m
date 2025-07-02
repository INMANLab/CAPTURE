%run from the Rasterplot GUI
dat = app.RWdat.allVars;
datEEG = dat.d_np;
% t = dat.ntp_np;
fs = 250;

t = (0:size(datEEG,1)-1)/fs;

winSize = 2; % in seconds
overlap = .5;


ch = 1;
d = datEEG(:,ch);
d(isnan(d)) = 0;
figure
subplot 211
plot(t,d);
subplot 212
spectrogram(d,winSize*fs,round(winSize*fs*overlap),[],fs);

