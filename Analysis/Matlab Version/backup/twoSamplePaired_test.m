function result = twoSamplePaired_test(group1, group2, patientList, nPerm, tail)
%TWOSAMPLEPAIRED_TEST Two-group paired (trial-level perm, cohort mean).
%
%   Stable alias for twoGroupPaired_trialPermutedCohortMean.
%   Signature: (group1, group2, patientList, nPerm, tail)

    if nargin < 4 || isempty(nPerm), nPerm = 1000; end
    if nargin < 5 || isempty(tail), tail = 'two-sided'; end

    result = permutation.twoGroupPaired_trialPermutedCohortMean( ...
        group1, group2, patientList, nPerm, tail);
end
