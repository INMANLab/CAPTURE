function result = oneSampleVsZero_test(group1, patientList, nPerm, tail)
%ONESAMPLEVSZERO_TEST One-sample vs zero (trial-level perm, cohort mean).
%
%   This is a stable alias for oneSampleVsZero_trialPermutedCohortMean.
%   Signature: (group1, patientList, nPerm, tail)

    if nargin < 3 || isempty(nPerm), nPerm = 1000; end
    if nargin < 4 || isempty(tail), tail = 'two-sided'; end

    result = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
        group1, patientList, nPerm, tail);
end
