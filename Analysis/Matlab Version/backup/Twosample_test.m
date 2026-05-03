function result = Twosample_test(group1, group2, patientList, nPerm, tail)
%TWOSAMPLE_TEST Backward-compatible name for cohort paired test.
%   result = permutation.Twosample_test(group1, group2, patientList, nPerm, tail)
%   forwards to twoGroupPaired_trialPermutedCohortMean.

    if nargin < 4 || isempty(nPerm), nPerm = 1000; end
    if nargin < 5 || isempty(tail), tail = 'two-sided'; end

    result = permutation.twoGroupPaired_trialPermutedCohortMean( ...
        group1, group2, patientList, nPerm, tail);
end
