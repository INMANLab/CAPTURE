%% Initialize
clear;
clc;
close all;

%% Load the data
GUIhandle = RWAnalysis2;
GUIhandle.AnalysisFile = 'D:\CAPTURE Project\GUI\Data\RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;

%% Call the desired GUI
MultTransSpecGramGUI(GUIhandle);