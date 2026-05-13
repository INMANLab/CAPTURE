function statisticValue = aggregateStratumMeans(muStratum, stratumParticipantIds, ...
    stratumChannelIds, aggregation)
%AGGREGATESTRATUMMEANS Collapse stratum-level means to one scalar test statistic.
%
%   aggregation "participant" : mean_p( mean_c( mu_{p,c} | same p ) )
%   aggregation "channel"     : mean_c( mean_p( mu_{p,c} | same c ) )

arguments
    muStratum (:,1) {mustBeNumeric}
    stratumParticipantIds (:,1)
    stratumChannelIds (:,1)
    aggregation (1,1) string
end

agg = validatestring(lower(aggregation), {'participant', 'channel'});

if numel(muStratum) ~= numel(stratumParticipantIds) || ...
        numel(muStratum) ~= numel(stratumChannelIds)
    error("aggregateStratumMeans:stratum vectors must match muStratum length.");
end

if strcmp(agg, 'participant')
    [~, ~, gp] = unique(stratumParticipantIds, 'stable');
    Up = max(gp);
    barByPrimary = zeros(Up, 1);
    for k = 1:Up
        barByPrimary(k) = mean(muStratum(gp == k), 'omitnan');
    end
    statisticValue = mean(barByPrimary, 'omitnan');
else
    [~, ~, gc] = unique(stratumChannelIds, 'stable');
    Uc = max(gc);
    barByPrimary = zeros(Uc, 1);
    for k = 1:Uc
        barByPrimary(k) = mean(muStratum(gc == k), 'omitnan');
    end
    statisticValue = mean(barByPrimary, 'omitnan');
end

end
