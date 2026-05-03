clear;
clc;

% Trial-within-participant permutation; cohort statistic = mean across
% participants of per-participant summaries. See README.md.

rng(42);
nPerm = 500;
tail = 'two-sided';

%% Smoke: identity permutation (column 1) reproduces observed cohort stat
nPatients = 4;
nTrials = 6;
N = nPatients * nTrials;
patientList = repelem((1:nPatients)', nTrials);
group1_smoke = randn(N, 1);
resSmoke = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
    group1_smoke, patientList, nPerm, tail);
assert(abs(resSmoke.nullDist(1) - resSmoke.obsStat) < 1e-10, ...
    'Identity permutation: first null should equal observed statistic.');

%% Tiny hand-traceable layout: 2 participants, 3 trials each
patientTiny = [1;1;1;2;2;2];
g1tiny = [1; 2; 3; 10; 11; 12];
resTiny = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
    g1tiny, patientTiny, 50, tail);
disp('Tiny one-sample (2 participants x 3 trials):');
disp(resTiny);

%% Synthetic cohort: known effect
nPatients = 10;
nObsPerPatient = 12;
N = nPatients * nObsPerPatient;
patientList = repelem((1:nPatients)', nObsPerPatient);
group1_effect = 0.35 + randn(N, 1);
group2_effect = randn(N, 1);
group1_null = randn(N, 1);
group2_null = randn(N, 1);

oneRes = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
    group1_effect, patientList, nPerm, tail);
disp('One-sample vs zero (trial perm, cohort mean of participant means):');
disp(oneRes);

twoRes = permutation.twoGroupPaired_trialPermutedCohortMean( ...
    group1_effect, group2_effect, patientList, nPerm, tail);
disp('Two-group paired (trial swap, cohort mean of participant mean diffs):');
disp(twoRes);

permutation.plotPermutationResult(oneRes, group1_effect, [], ...
    'One-sample cohort test', patientList);
permutation.plotPermutationResult(twoRes, group1_effect, group2_effect, ...
    'Two-group cohort test', patientList);

%% Across-participant entry (same as cohort trial-perm functions on raw trials)
acrossOne = permutation.runAcrossParticipantsPermutation( ...
    group1_effect, [], patientList, 'one-sample-vs-zero', nPerm, [], tail, false, 0.05);
acrossTwo = permutation.runAcrossParticipantsPermutation( ...
    group1_effect, group2_effect, patientList, 'two-sample-paired', nPerm, [], tail, false, 0.05);
disp('runAcrossParticipantsPermutation (delegates to cohort trial-perm):');
disp(acrossOne.overall.pValue);
disp(acrossTwo.overall.pValue);

%% Null sanity (expect less extreme stats / higher p-values on average)
nullOne = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
    group1_null, patientList, nPerm, tail);
nullTwo = permutation.twoGroupPaired_trialPermutedCohortMean( ...
    group1_null, group2_null, patientList, nPerm, tail);
disp('Null cohort checks:');
disp(struct('oneSampleP', nullOne.pValue, 'twoSampleP', nullTwo.pValue));

%% Optional example dataset
if exist('RWA_PowerBox_Example.mat', 'file') == 2
    load('RWA_PowerBox_Example.mat', 'pwr_box1', 'pwr_box2', 'patient');
    if exist('pwr_box1', 'var') && exist('pwr_box2', 'var') && exist('patient', 'var')
        cohortLegacy = permutation.twoGroupPaired_trialPermutedCohortMean( ...
            pwr_box1(:), pwr_box2(:), patient(:), nPerm, tail);
        disp('Example dataset cohort paired test:');
        disp(cohortLegacy);
    end
end
