function [nullDistribution, pValue, statisticValue, info] = TwoSample( ...
    group1, group2, patientList, args)
%TWOSAMPLE Participant-level two-condition permutation test.
%
%   [nullDistribution, pValue, statisticValue] = permutation.TwoSample( ...
%       group1, group2, patientList)
%   [...] = permutation.TwoSample(group1, group2, patientList, nPerm=2000, Seed=42)
%   [..., info] = permutation.TwoSample(...)
%
%   group1, group2, patientList must be Nx1 with the same N (paired trial rows).
%   Without patientChannel: per participant, pool both conditions and randomly split
%   into two equal-sized sets; statistic is mean across participants of (mean_g1 -
%   mean_g2).
%
%   With patientChannel: same pooling/split within each (participant, channel)
%   stratum. Aggregation collapses stratum-wise (mean_g1 - mean_g2) to one scalar
%   ("participant" or "channel"); see OneSampleVsZero help.

arguments
    group1 (:,1) {mustBeNumeric}
    group2 (:,1) {mustBeNumeric}
    patientList (:,1)
    args.nPerm (1,1) double {mustBePositive, mustBeInteger} = 1000
    args.Seed (1,1) double = NaN
    args.patientChannel (:,1) = []
    args.Aggregation (1,1) string = "participant"
end

nPerm = args.nPerm;
if ~isnan(args.Seed)
    validateattributes(args.Seed, {'double', 'single'}, {'scalar', 'integer', 'nonnegative'});
    rng(args.Seed, 'twister');
end

if numel(group1) ~= numel(group2) || numel(group1) ~= numel(patientList)
    error("TwoSample:group1, group2, and patientList must have the same length.");
end

useNested = ~isempty(args.patientChannel);
if useNested && numel(args.patientChannel) ~= numel(group1)
    error("TwoSample:patientChannel must have the same length as group1.");
end

rngStream = RandStream.getGlobalStream();

if useNested
    [stratumParticipantIds, stratumChannelIds, ~] = ...
        permutation.participantChannelStrata(patientList, args.patientChannel);
    mu1s = permutation.stratumMeansByParticipantChannel( ...
        group1, patientList, args.patientChannel);
    mu2s = permutation.stratumMeansByParticipantChannel( ...
        group2, patientList, args.patientChannel);
    muDiff = mu1s - mu2s;
    statisticValue = permutation.aggregateStratumMeans( ...
        muDiff, stratumParticipantIds, stratumChannelIds, args.Aggregation);
    mu1 = [];
    mu2 = [];
else
    stratumParticipantIds = [];
    stratumChannelIds = [];
    muDiff = [];
    mu1 = permutation.participantMeansById(group1, patientList);
    mu2 = permutation.participantMeansById(group2, patientList);
    statisticValue = mean(mu1 - mu2, 'omitnan');
end

nullDistribution = zeros(nPerm, 1);
for p = 1:nPerm
    if useNested
        [g1p, g2p] = permutation.permuteLabelsWithinParticipants( ...
            'twoSample', group1, patientList, group2, rngStream, args.patientChannel);
        m1s = permutation.stratumMeansByParticipantChannel( ...
            g1p, patientList, args.patientChannel);
        m2s = permutation.stratumMeansByParticipantChannel( ...
            g2p, patientList, args.patientChannel);
        nullDistribution(p) = permutation.aggregateStratumMeans( ...
            m1s - m2s, stratumParticipantIds, stratumChannelIds, args.Aggregation);
    else
        [g1p, g2p] = permutation.permuteLabelsWithinParticipants( ...
            'twoSample', group1, patientList, group2, rngStream);
        m1 = permutation.participantMeansById(g1p, patientList);
        m2 = permutation.participantMeansById(g2p, patientList);
        nullDistribution(p) = mean(m1 - m2, 'omitnan');
    end
end

pValue = permutation.twoSidedPValue(statisticValue, nullDistribution);

if nargout > 3
    info = struct();
    info.nPerm = nPerm;
    info.feasibility = permutation.Feasibility( ...
        'twoSample', patientList, nPerm, group1, group2, args.patientChannel);
    [info.participantIds, ~] = permutation.participantTrialIndices(patientList);
    if useNested
        info.stratumMeansObservedG1 = mu1s;
        info.stratumMeansObservedG2 = mu2s;
        info.stratumParticipantIds = stratumParticipantIds;
        info.stratumChannelIds = stratumChannelIds;
        info.Aggregation = args.Aggregation;
        info.participantMeansObservedG1 = permutation.stratumToParticipantMeans( ...
            mu1s, stratumParticipantIds, info.participantIds);
        info.participantMeansObservedG2 = permutation.stratumToParticipantMeans( ...
            mu2s, stratumParticipantIds, info.participantIds);
    else
        info.participantMeansObservedG1 = mu1;
        info.participantMeansObservedG2 = mu2;
    end
end

end
