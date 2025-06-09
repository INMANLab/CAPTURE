clear;
clc;

GUIhandle = RWAnalysis2;
% GUIhandle.AnalysisFile = 'D:\CAPTURE Project\GUI\Data\RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;


MultTransSpecGramGUI(GUIhandle);