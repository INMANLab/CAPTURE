clear;
clc;

% Example script for the permutation package API.
% Requires: RWA_PowerBox_Example.mat with variables: patient, pwr_box1, pwr_box2
load RWA_PowerBox_Example.mat;

patientList = patient(:);
group1 = pwr_box1(:);
group2 = pwr_box2(:);

% Reproducible randomization for demo runs.
rng(7);
nPerm = 1000;

permutation.FeasibilityCheck('oneSample',patientList, nPerm);
% One-sample vs 0 (trial-level sign flip within participant)
[nullDist1, p1, stat1, info1] = permutation.OneSampleVsZero( ...
    group1, patientList, nPerm=nPerm, Seed=42);

permutation.PlotNullDistribution(nullDist1, stat1, PValue=p1);
permutation.PlotParticipantAveragesBoxplot(p1, stat1, info1, ...
    Title="One-sample vs zero");


% One-sample vs 0 (trial-level sign flip within participant)
[nullDist1, p1, stat1, info1] = permutation.OneSampleVsZero( ...
    group2, patientList, nPerm=nPerm, Seed=42);

permutation.PlotNullDistribution(nullDist1, stat1, PValue=p1);
permutation.PlotParticipantAveragesBoxplot(p1, stat1, info1, ...
    Title="One-sample vs zero");

% Two-sample (pool group1 & group2 trials per participant, random equal split)
permutation.FeasibilityCheck('twoSample',patientList, nPerm);
[nullDist2, p2, stat2, info2] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=nPerm, Seed=42);

permutation.PlotNullDistribution(nullDist2, stat2, PValue=p2);
permutation.PlotParticipantAveragesBoxplot(p2, stat2, info2, ...
    GroupLabels=["Group 1", "Group 2"], ...
    Title="Two-sample");

% General Feasibility check for two-sample label-shuffle null
out = permutation.Feasibility('twoSample', patientList, 1000, group1, group2);
disp(out.message);