clear;
clc;

% This script demonstrates the modular permutation testing API.
% It includes deterministic synthetic smoke tests and an optional example
% using RWA_PowerBox_Example.mat when that dataset is present.
load RWA_PowerBox_Example.mat
patientList= patient;
group1 = pwr_box1;
group2 = pwr_box2;


rng(7); % Reproducible randomization for smoke tests.
nPerm = 1000;
statType = 'mean';
tail = 'two-sided';


resWithin = permutation.OneSample_ParticipantLevel( ...
    patientList, group1, nPerm, statType, tail, true, 0.05);
disp(resWithin);
permutation.plotPermutationResult(resWithin, group1, [], ...
    'One-Sample vs Zero (Synthetic Effect)');
resAcross = permutation.OneSample_AcrossParticipants( ...
    patientList, group1, nPerm, statType, tail, true, 0.05);




oneSampleRes = permutation.oneSampleVsZero_test(group1, nPerm, statType, tail);
disp('One-sample-vs-zero result:');
disp(oneSampleRes);
permutation.plotPermutationResult(oneSampleRes, group1, [], ...
    'One-Sample vs Zero (Synthetic Effect)');

