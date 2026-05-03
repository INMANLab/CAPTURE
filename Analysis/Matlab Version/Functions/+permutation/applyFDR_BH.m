function [pAdj, isSignificant, criticalP] = applyFDR_BH(pValues, q)
%APPLYFDR_BH Benjamini-Hochberg false discovery rate correction.
%   [pAdj, isSignificant, criticalP] = permutation.applyFDR_BH(pValues, q)
%   performs BH correction on a vector of p-values.

    arguments
        pValues (:,1) double
        q (1,1) double {mustBeGreaterThan(q, 0), mustBeLessThan(q, 1)} = 0.05
    end

    validMask = ~isnan(pValues);
    pValid = pValues(validMask);
    m = numel(pValid);

    pAdj = nan(size(pValues));
    isSignificant = false(size(pValues));
    criticalP = NaN;

    if m == 0
        return;
    end

    [pSorted, sortIdx] = sort(pValid);
    ranks = (1:m)';

    thresholds = (ranks / m) * q;
    below = pSorted <= thresholds;
    if any(below)
        maxIdx = find(below, 1, 'last');
        criticalP = pSorted(maxIdx);
    end

    % Monotone adjusted p-values from BH procedure.
    pAdjSorted = (m ./ ranks) .* pSorted;
    pAdjSorted = min(1, pAdjSorted);
    for i = m-1:-1:1
        pAdjSorted(i) = min(pAdjSorted(i), pAdjSorted(i+1));
    end

    invSort(sortIdx) = 1:m; %#ok<AGROW>
    pAdjValid = pAdjSorted(invSort(:));
    sigValid = pAdjValid <= q;

    pAdj(validMask) = pAdjValid;
    isSignificant(validMask) = sigValid;
end
