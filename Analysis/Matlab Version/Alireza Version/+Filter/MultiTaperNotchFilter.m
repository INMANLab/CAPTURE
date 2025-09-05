function sigF = MultiTaperNotchFilter(sig,fCut,Fs,plotFlag)

params.Fs = Fs;

winSize = length(sig);
T = winSize/Fs; % in Seconds
W = 3; % in Hz, this is the W in multitaper language [-W W]

% movingwin = [5, 2.5];
TW = T*W;
[e,v] = dpss(winSize,TW,50);
k = find(v>0.99,1,'last');
if(isempty(k))
    k=1;
end
% figure
% lv = length(v);
% stem(1:lv,v,"filled")
% title("Proportion of Energy in [-w,w] of k-th Slepian Sequence")
% 


params.Fs = Fs;
params.tapers = [TW k];
params.pad = 1;
% [S1,f]=mtspectrumc(eegDatF,params);
% datafit = fitlinesc(eegDat,params,0,'y',62.5);
% eegDatF = eegDat;
% eegDatF = eegDatF-2*datafit;
sigF= rmlinesc(sig, params, 0, plotFlag, fCut);

% params.tapers = [TW k];
% params.pad = 1;
% params.trialave = 0;
% tau = 10;
% 
% eegDatF= rmlinesmovingwinc(eegDat, movingwin, tau, params, 0,'',62.5,true);


% figure
% eegDatF= rmlinesmovingwinc(eegDat, movingwin, tau, params, 0,'y',62.5,true);
% title("RW"+pIdx+"_Walk"+wIdx+" Channel: "+chIdx+" Chunk:"+chunk)


% params.fpass = [62,63];
% [S1,f]=mtspectrumc(eegDat,params);

% params.Fs = Fs;
% movingwin = [5, 2.5];
% params.fpass = [62,63];
% params.tapers = [3 5];
% params.pad = 1;
% [S1,f]=mtspectrumc(eegDatF,params);
% [datafit,Amps]=fitlinesc(eegDat,params,0,'y',62.5);
% eegDatF = eegDat;
% eegDatF = eegDatF-2*datafit;
% eegDatF= rmlinesc(eegDat, params, 0,'y',62.5);
% figure;plot(datafit)