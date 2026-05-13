function [nullDistribution, pValue, statisticValue, info] = OneSampleVsZero( ...
    group1, patientList, args)
%PERMUTATIONONESAMPLEVSZERO Participant-mean test vs 0 with trial-level sign-flip null.
%
%   [nullDistribution, pValue, statisticValue] = permutation.OneSampleVsZero( ...
%       group1, patientList)
%   [...] = permutation.OneSampleVsZero(group1, patientList, nPerm=2000, Seed=42)
%   [..., info] = permutation.OneSampleVsZero(...)
%
%   Without patientChannel: sign-flip i.i.d. within each participant; statistic is
%   mean across participants of mean trial value (per participant).
%
%   With patientChannel (Nx1, aligned with trials): sign-flip within each
%   (participant, channel) stratum. Use Aggregation to collapse stratum means to
%   one scalar: "participant" (mean over channels per participant, then over
%   participants) or "channel" (mean over participants per channel, then over
%   channels).
%
%   Optional name-values (after first two arguments):
%     nPerm           (default 1000)
%     Seed            nonnegative integer; sets rng before permutations (default: NaN)
%     patientChannel  Nx1 optional; nested permutation strata
%     Aggregation     "participant" | "channel" (default "participant"; ignored if no channel)

arguments
    group1 (:,1) {mustBeNumeric}
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

if numel(group1) ~= numel(patientList)
    error('OneSampleVsZero:group1 and patientList must have the same length.');
end

useNested = ~isempty(args.patientChannel);
if useNested && numel(args.patientChannel) ~= numel(group1)
    error('OneSampleVsZero:patientChannel must have the same length as group1.');
end

rngStream = RandStream.getGlobalStream();

if useNested
    [stratumParticipantIds, stratumChannelIds, ~] = ...
        permutation.participantChannelStrata(patientList, args.patientChannel);
    muStratumObs = permutation.stratumMeansByParticipantChannel( ...
        group1, patientList, args.patientChannel);
    statisticValue = permutation.aggregateStratumMeans( ...
        muStratumObs, stratumParticipantIds, stratumChannelIds, args.Aggregation);
else
    muStratumObs = [];
    stratumParticipantIds = [];
    stratumChannelIds = [];
    muObs = permutation.participantMeansById(group1, patientList);
    statisticValue = mean(muObs, 'omitnan');
end

nullDistribution = zeros(nPerm, 1);
for p = 1:nPerm
    if useNested
        g1p = permutation.permuteLabelsWithinParticipants( ...
            'oneSample', group1, patientList, [], rngStream, args.patientChannel);
        muP = permutation.stratumMeansByParticipantChannel( ...
            g1p, patientList, args.patientChannel);
        nullDistribution(p) = permutation.aggregateStratumMeans( ...
            muP, stratumParticipantIds, stratumChannelIds, args.Aggregation);
    else
        g1p = permutation.permuteLabelsWithinParticipants( ...
            'oneSample', group1, patientList, [], rngStream);
        muP = permutation.participantMeansById(g1p, patientList);
        nullDistribution(p) = mean(muP, 'omitnan');
    end
end

pValue = permutation.twoSidedPValue(statisticValue, nullDistribution);

if nargout > 3
    info = struct();
    info.nPerm = nPerm;
    info.feasibility = permutation.Feasibility( ...
        'oneSample', patientList, nPerm, group1, [], args.patientChannel);
    [info.participantIds, ~] = permutation.participantTrialIndices(patientList);
    if useNested
        info.stratumMeansObserved = muStratumObs;
        info.stratumParticipantIds = stratumParticipantIds;
        info.stratumChannelIds = stratumChannelIds;
        info.Aggregation = args.Aggregation;
        info.participantMeansObserved = permutation.stratumToParticipantMeans( ...
            muStratumObs, stratumParticipantIds, info.participantIds);
    else
        info.participantMeansObserved = muObs;
    end
end

end
