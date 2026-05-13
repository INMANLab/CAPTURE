clear;
clc;

% Example script for the permutation package API.
% Requires: RWA_PowerBox_Example.mat with variables: patient, pwr_box1, pwr_box2
load RWA_PowerBox_Example.mat;

patientList = patient(:);
group1 = 10*log10(pwr_box1(:));
group2 = 10*log10(pwr_box2(:));
patientChannel = channels(:);

% Reproducible randomization for demo runs.
rng(7);
nPerm = 1000;

%%
permutation.FeasibilityCheck('oneSample', patientList, nPerm);
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
permutation.FeasibilityCheck('twoSample', patientList, nPerm);
[nullDist2, p2, stat2, info2] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=nPerm, Seed=42);

permutation.PlotNullDistribution(nullDist2, stat2, PValue=p2);
permutation.PlotParticipantAveragesBoxplot(p2, stat2, info2, ...
    GroupLabels=["Group 1", "Group 2"], ...
    Title="Two-sample");

%%
% Nested (participant x channel): permute within each stratum; aggregate statistic
permutation.FeasibilityCheck('oneSample', patientList, nPerm, patientChannel);
[nullNest1, pNest1, statNest1, infoNest1] = permutation.OneSampleVsZero( ...
    group1, patientList, nPerm=nPerm, Seed=42, ...
    patientChannel=patientChannel, Aggregation="participant");
permutation.PlotNullDistribution(nullNest1, statNest1, PValue=pNest1);
permutation.PlotParticipantAveragesBoxplot(pNest1, statNest1, infoNest1, ...
    Title="One-sample nested (participant x channel)",PlotPerChannelStratumPoints=true);



permutation.FeasibilityCheck('oneSample', patientList, nPerm, patientChannel);
[nullNest2, pNest2, statNest2, infoNest2] = permutation.OneSampleVsZero( ...
    group2, patientList, nPerm=nPerm, Seed=42, ...
    patientChannel=patientChannel, Aggregation="participant");
permutation.PlotNullDistribution(nullNest2, statNest2, PValue=pNest2);
permutation.PlotParticipantAveragesBoxplot(pNest2, statNest2, infoNest2, ...
    Title="One-sample nested (participant x channel)",PlotPerChannelStratumPoints=true);

permutation.FeasibilityCheck('twoSample', patientList, nPerm, patientChannel);
[nullNest2, pNest2, statNest2, infoNest2] = permutation.TwoSample( ...
    group1, group2, patientList, nPerm=nPerm, Seed=42, ...
    patientChannel=patientChannel, Aggregation="participant");
permutation.PlotNullDistribution(nullNest2, statNest2, PValue=pNest2);
permutation.PlotParticipantAveragesBoxplot(pNest2, statNest2, infoNest2, ...
    GroupLabels=["Group 1", "Group 2"], ...
    Title="Two-sample nested (participant x channel)", ...
    PlotPerChannelStratumPoints=true, ParticipantMarkerCycle=["^", "s", "d", "v"]);

% General Feasibility check for two-sample label-shuffle null
out = permutation.Feasibility('twoSample', patientList, 1000, group1, group2);
disp(out.message);
outNested = permutation.Feasibility('twoSample', patientList, 1000, group1, group2, patientChannel);
disp(outNested.message);
