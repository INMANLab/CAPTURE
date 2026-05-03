function result = twoSamplePaired_test(group1, group2, nPerm, statType, tail)
%TWOSAMPLEPAIRED_TEST Paired two-sample permutation test.
%   result = permutation.twoSamplePaired_test(group1, group2, nPerm, statType, tail)
%   tests group1 vs group2 with within-pair label swaps.

    arguments
        group1 (:,1) double
        group2 (:,1) double
        nPerm (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPerm, 1)} = 1000
        statType (1,:) char = 'mean'
        tail (1,:) char = 'two-sided'
    end

    statType = validatestring(lower(statType), {'mean','median','t'});
    tail = validatestring(lower(tail), {'two-sided','right','left'});

    if numel(group1) ~= numel(group2)
        error('permutation:twoSamplePaired_test:DimensionMismatch', ...
            'group1 and group2 must have the same number of observations.');
    end

    validMask = ~isnan(group1) & ~isnan(group2);
    g1 = group1(validMask);
    g2 = group2(validMask);
    nObs = numel(g1);

    if nObs == 0
        error('permutation:twoSamplePaired_test:NoData', ...
            'No valid paired observations after removing NaNs.');
    end

    swapMask = permutation.generatePermutationLabels(nObs, nPerm, 'paired_swap');
    nullDist = zeros(nPerm, 1);

    for iPerm = 1:nPerm
        currentSwap = logical(swapMask(:, iPerm));
        permG1 = g1;
        permG2 = g2;

        tmp = permG1(currentSwap);
        permG1(currentSwap) = permG2(currentSwap);
        permG2(currentSwap) = tmp;

        nullDist(iPerm) = permutation.computeStatistic(permG1, permG2, statType);
    end

    obsStat = permutation.computeStatistic(g1, g2, statType);
    pValue = localPermutationPValue(nullDist, obsStat, tail);

    result = struct( ...
        'nullDist', nullDist, ...
        'obsStat', obsStat, ...
        'pValue', pValue, ...
        'pValueAdj', NaN, ...
        'isSignificant', false, ...
        'criticalP', NaN, ...
        'statType', statType, ...
        'tail', tail, ...
        'nPerm', nPerm, ...
        'testType', 'two-sample-paired', ...
        'metadata', struct('nObs', nObs));
end

function pVal = localPermutationPValue(nullDist, obsStat, tail)
    switch tail
        case 'two-sided'
            pVal = (sum(abs(nullDist) >= abs(obsStat)) + 1) / (numel(nullDist) + 1);
        case 'right'
            pVal = (sum(nullDist >= obsStat) + 1) / (numel(nullDist) + 1);
        case 'left'
            pVal = (sum(nullDist <= obsStat) + 1) / (numel(nullDist) + 1);
    end
end
