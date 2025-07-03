

pID = 1;
wID = 1;

pData = GUIhandle.MultTable{pID};
ntpDat = pData.ntp_np{1};
% RasterPlot = RasterPlotGUI;

% events = GUIhandle.MultTrans.Evnt;
% walks = GUIhandle.MultTrans.Walk;
% patient = GUIhandle.MultTrans.patient;
% kdIdx = GUIhandle.MultTrans.DT_kd_idx;  
%%  Extract Events

patient = GUIhandle.MultTrans.Patient;
walks = GUIhandle.MultTrans.Walk(patient == pID);
evTable = GUIhandle.MultTrans.Evnt(patient == pID);
kdIdx = GUIhandle.MultTrans.DT_kd_idx(patient == pID);
npIdx = GUIhandle.MultTrans.DT_np_idx(patient == pID);
walkNums = unique(walks);
for i=1:max(walkNums)
    wI = walks(walks==i);
    eI = evTable(walks==i);
    nI = npIdx(walks==i);
    kI = kdIdx(walks==i);
    EvList = unique(eI);
    for evIdx = 1:10
    end
end

%% Use the individual Data

files = dir("*.mat");
fileNames = string({files.name})';
cell2mat(regexp(fileNames,'.*RW.*_Walk.*','start'));


files = dir('*.mat');
fileNames = {files.name};
tokens = regexp(fileNames, 'RW(\d+)_Walk(\d+)\.mat', 'tokens');
tokens = vertcat(tokens{:});
patientNumbers = str2double(cellfun(@(x) x(1), tokens));
walkNumbers = str2double(cellfun(@(x) x(2), tokens));  



dir("*RW.*_Walk.*.mat")
clear;
clc;
pID = 1;
wID = 1;
minCounts = 8;
dat = load("RWNApp_RW"+pID+"_Walk"+wID+".mat");
evTable = dat.evnts_tbl;

tArray = (evTable.NTP-evTable.NTP(1))/60;
evCounts = groupsummary(evTable,"Event");
evUnique = evCounts.Event(evCounts.GroupCount>minCounts);

evTable.Description2 = regexprep(evTable.Description, '.*?\[(.*?)\].*?', '[$1]');
allLabels = regexp(evTable.Description, '\[([^\]]+)\]', 'tokens');
allLabels = [allLabels{:}];
uniqueLabels = unique(allLabels);

figure
hold on
for evIdx = 1:length(evUnique)
    evName = evUnique(evIdx);
    evT = tArray(strcmp(evTable.Event,evName));
    plot(evT,ones(1,length(evT))*evIdx,"|",LineWidth=1);
end
yticks(1:length(evUnique))
yticklabels(evUnique)
xlabel("time(min)")
axis tight
%%





plotwidth=1;     % spike thickness
plotcolor='k';   % spike color
trialgap=1.5;    % distance between trials
defaultfs=1000;  % default sampling rate
showtimescale=1; % display timescale
showlabels=1;    % display x and y labels

figure;
hresp=gca;
fs=defaultfs;

t=[10 250 9000 1300,1600,2405,2900];
times = t;
numtrials = 3;
triallen = 1000;

trials=ceil(times/triallen);
reltimes=mod(times,triallen);
reltimes(~reltimes)=triallen;
numspikes=length(times);
xx=ones(3*numspikes,1)*nan;
yy=ones(3*numspikes,1)*nan;
yy(1:3:3*numspikes)=(trials-1)*trialgap;
yy(2:3:3*numspikes)=yy(1:3:3*numspikes)+1;

%scale the time axis to ms
xx(1:3:3*numspikes)=reltimes*1000/fs;
xx(2:3:3*numspikes)=reltimes*1000/fs;
xlim=[1,triallen*1000/fs];


axes(hresp);
h=plot(xx, yy, plotcolor, 'linewidth',plotwidth);
axis ([xlim,0,(numtrials)*1.5]);  


%%
A = squareform(pdist(vals','correlation'));
corrMat = corr(vals);
sumCorr = sum(corrMat,2);
[~, sortIdx] = sort(sumCorr, 'descend');
tree = linkage(vals',"average",'correlation');
D = pdist(vals','correlation');
leafOrder = optimalleaforder(tree,D);
dendrogram(tree,Reorder=leafOrder)