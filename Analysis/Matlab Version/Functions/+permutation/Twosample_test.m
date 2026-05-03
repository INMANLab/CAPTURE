function result = Twosample_test(group1, group2, nPerm, statType, tail)
%TWOSAMPLE_TEST Backward-compatible wrapper for paired two-sample testing.
%   This legacy name forwards to permutation.twoSamplePaired_test.

    if nargin < 3 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 4 || isempty(statType)
        statType = 'mean';
    end
    if nargin < 5 || isempty(tail)
        tail = 'two-sided';
    end

    result = permutation.twoSamplePaired_test(group1, group2, nPerm, statType, tail);
end