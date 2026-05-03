function result = oneSampleVsZero_test(group1, nPerm, statType, tail)
%ONESAMPLEVSZERO_TEST One-sample permutation test versus zero.
%   result = permutation.oneSampleVsZero_test(group1, nPerm, statType, tail)
%   computes an observed statistic and permutation-based null distribution
%   using sign flips.

    arguments
        group1 (:,1) double
        nPerm (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(nPerm, 1)} = 1000
        statType (1,:) char = 'mean'
        tail (1,:) char = 'two-sided'
    end

    statType = validatestring(lower(statType), {'mean','median','t'});
    tail = validatestring(lower(tail), {'two-sided','right','left'});

    cleanData = group1(~isnan(group1));
    nObs = numel(cleanData);
    if nObs == 0
        error('permutation:oneSampleVsZero_test:NoData', ...
            'group1 has no valid (non-NaN) observations.');
    end

    labels = permutation.generatePermutationLabels(nObs, nPerm, 'onesample_signflip');
    nullDist = zeros(nPerm, 1);
    for iPerm = 1:nPerm
        permuted = cleanData .* labels(:, iPerm);
        nullDist(iPerm) = permutation.computeStatistic(permuted, [], statType);
    end

    obsStat = permutation.computeStatistic(cleanData, [], statType);
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
        'testType', 'one-sample-vs-zero', ...
        'metadata', struct('nObs', nObs));
end

function pVal = localPermutationPValue(nullDist, obsStat, tail)
% Conservative +1 correction avoids zero p-values.
    switch tail
        case 'two-sided'
            pVal = (sum(abs(nullDist) >= abs(obsStat)) + 1) / (numel(nullDist) + 1);
        case 'right'
            pVal = (sum(nullDist >= obsStat) + 1) / (numel(nullDist) + 1);
        case 'left'
            pVal = (sum(nullDist <= obsStat) + 1) / (numel(nullDist) + 1);
    end
end
