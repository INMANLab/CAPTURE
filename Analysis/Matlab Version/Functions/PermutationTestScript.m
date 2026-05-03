clear;
clc;

% This script demonstrates the modular permutation testing API.
% It includes deterministic synthetic smoke tests and an optional example
% using RWA_PowerBox_Example.mat when that dataset is present.

rng(7); % Reproducible randomization for smoke tests.
nPerm = 1000;
statType = 'mean';
tail = 'two-sided';

%% Synthetic data setup (multiple observations per participant)
nPatients = 10;
nObsPerPatient = 12;
N = nPatients * nObsPerPatient;
patientList = repelem((1:nPatients)', nObsPerPatient);

% Known-effect paired data: group1 has positive shift relative to group2.
group1_effect = 0.5 + randn(N, 1);
group2_effect = randn(N, 1);

% Null data: no systematic shift.
group1_null = randn(N, 1);
group2_null = randn(N, 1);

%% 1) One-sample versus zero
oneSampleRes = permutation.oneSampleVsZero_test(group1_effect, nPerm, statType, tail);
disp('One-sample-vs-zero result:');
disp(oneSampleRes);
permutation.plotPermutationResult(oneSampleRes, group1_effect, [], ...
    'One-Sample vs Zero (Synthetic Effect)');

%% 2) Paired two-sample
twoSampleRes = permutation.twoSamplePaired_test(group1_effect, group2_effect, nPerm, statType, tail);
disp('Two-sample paired result (known effect):');
disp(twoSampleRes);
permutation.plotPermutationResult(twoSampleRes, group1_effect, group2_effect, ...
    'Paired Two-Sample (Synthetic Effect)');

%% 3) Participant-level multiple tests + FDR
participantRes = permutation.runParticipantPermutation( ...
    group1_effect, group2_effect, patientList, ...
    'two-sample-paired', nPerm, statType, tail, true, 0.05);
disp('Participant-level paired test results (FDR-adjusted):');
disp(participantRes);
permutation.plotParticipantPermutationSummary( ...
    participantRes, group1_effect, group2_effect, patientList, ...
    'Participant-Level Summary (Synthetic Effect)');

% Collect participant-level p-values for compact inspection.
pRaw = arrayfun(@(s) s.pValue, participantRes.perPatient);
pAdj = arrayfun(@(s) s.pValueAdj, participantRes.perPatient);
disp(table((1:numel(pRaw))', pRaw(:), pAdj(:), 'VariableNames', ...
    {'PatientIdx','PValueRaw','PValueFDR'}));

%% 4) Across-participant aggregated test
acrossRes = permutation.runAcrossParticipantsPermutation( ...
    group1_effect, group2_effect, patientList, ...
    'two-sample-paired', nPerm, statType, tail, true, 0.05);
disp('Across-participant aggregated paired test:');
disp(acrossRes);

%% 5) Null sanity checks (expect larger p-values on average)
nullOneSampleRes = permutation.oneSampleVsZero_test(group1_null, nPerm, statType, tail);
nullTwoSampleRes = permutation.twoSamplePaired_test(group1_null, group2_null, nPerm, statType, tail);
disp('Null sanity checks:');
disp(struct('oneSampleP', nullOneSampleRes.pValue, 'twoSampleP', nullTwoSampleRes.pValue));

%% Optional legacy dataset demo
if exist('RWA_PowerBox_Example.mat', 'file') == 2
    load('RWA_PowerBox_Example.mat', 'pwr_box1', 'pwr_box2', 'patient');
    if exist('pwr_box1', 'var') && exist('pwr_box2', 'var') && exist('patient', 'var')
        legacyRes = permutation.runParticipantPermutation( ...
            pwr_box1(:), pwr_box2(:), patient(:), ...
            'two-sample-paired', nPerm, statType, tail, true, 0.05);
        disp('Legacy dataset participant-level paired result:');
        disp(legacyRes);
    end
end