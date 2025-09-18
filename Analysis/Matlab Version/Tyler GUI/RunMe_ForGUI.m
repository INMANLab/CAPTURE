%% Initialize
clear;
clc;
close all;

%% Load the data
GUIhandle = RWAnalysis2;
GUIhandle.AnalysisFile = 'Z:\Tyler\RealWorld\RWAnalysis2_MedianNorm_2sec.mat';
GUIhandle.loadMultData;

%% Call the desired GUI
MultTransSpecGramGUI(GUIhandle);