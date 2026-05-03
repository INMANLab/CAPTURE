function statValue = computeStatistic(x, y, statType)
%COMPUTESTATISTIC Compute scalar statistics for permutation tests.
%   statValue = permutation.computeStatistic(x, y, statType) computes a
%   scalar statistic according to statType.
%
%   Inputs
%   ------
%   x        : Numeric vector (required)
%   y        : Numeric vector or [].
%              - If empty, one-sample statistic is computed on x.
%              - If provided, statistic is computed on paired differences x-y.
%   statType : 'mean' | 'median' | 't'
%
%   Output
%   ------
%   statValue: Scalar observed statistic.

    arguments
        x (:,1) double
        y (:,1) double = []
        statType (1,:) char = 'mean'
    end

    statType = validatestring(lower(statType), {'mean','median','t'});

    if isempty(y)
        data = x;
    else
        if numel(x) ~= numel(y)
            error('permutation:computeStatistic:DimensionMismatch', ...
                'x and y must have the same length for paired statistics.');
        end
        data = x - y;
    end

    data = data(~isnan(data));
    if isempty(data)
        statValue = NaN;
        return;
    end

    switch statType
        case 'mean'
            statValue = mean(data);
        case 'median'
            statValue = median(data);
        case 't'
            n = numel(data);
            sampleStd = std(data, 0);
            if n < 2 || sampleStd == 0
                statValue = 0;
            else
                statValue = mean(data) / (sampleStd / sqrt(n));
            end
    end
end
