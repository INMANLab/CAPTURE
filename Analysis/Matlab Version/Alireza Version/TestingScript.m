clear;
clc;

GUIhandle = RWAnalysis2;
% GUIhandle.AnalysisFile = 'D:\CAPTURE Project\GUI\Data\RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;


% h = GUIhandle.openMultFRankingGUI;
% 
% GUIhandle.openMultSegGLMEGUI;

MultTransSpecGramGUI_AK(GUIhandle);
GUIhandle.openMultTransSpecGramGUI;

%%
pHTemp = pH;
darkGrayShades = ["#333333", "#2F4F4F", "#2E2E2E", "#282828", "#232323", "#1C1C1C", "#171717", "#0F0F0F"];
lineInfo = struct();
for mIdx=1:size(pHTemp,2)
    yDat = pHTemp(3,mIdx).YData;
    mVal = mean(yDat);
    sVal = std(yDat);
    pHTemp(3,mIdx).YData = (yDat-mVal)/sVal;
    pHTemp(3,mIdx).Color = darkGrayShades(mIdx);
    lineInfo.mVal(mIdx) = mVal;
    lineInfo.sdVal(mIdx) = sVal;
    lineInfo.maxVal(mIdx) = max(pHTemp(3,mIdx).YData);
    lineInfo.minVal(mIdx) = min(pHTemp(3,mIdx).YData);
end
orderPlot = str2num(app.MeasureIndexesEditField.Value);
for idx = (length(orderPlot)-1):-1:1
    mIdx0 = orderPlot(idx+1);
    mIdx = orderPlot(idx);
    mAmplitude = lineInfo.maxVal(mIdx)-lineInfo.minVal(mIdx);
    pHTemp(3,mIdx).YData = pHTemp(3,mIdx).YData+lineInfo.maxVal(mIdx0)+mAmplitude;%abs(lineInfo.minVal(mIdx));
end
yLimMax = max(pHTemp(3,mIdx).YData);
yLimMin = min(pHTemp(3,orderPlot(end)).YData);
set(pH(3,:),'visible','off');
set(pH(3,[7,8,4]),'visible','on');
ylim(app.aH,[yLimMin,yLimMax])


% yyaxis(app.aH,'left');
% xlim(app.aH,str2num(app.XlimEditField.Value));