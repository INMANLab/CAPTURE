function [nullDistribution, pValue, statisticValue, info] = permutationOneSampleVsZero( ...
    group1, patientList, args)
%PERMUTATIONONESAMPLEVSZERO Participant-mean test vs 0 with trial-level sign-flip null.
%
%   [nullDistribution, pValue, statisticValue] = permutationOneSampleVsZero( ...
%       group1, patientList)
%   [...] = permutationOneSampleVsZero(group1, patientList, nPerm)
%   [...] = permutationOneSampleVsZero(group1, patientList, ..., nPerm=2000, Seed=42)
%   [..., info] = permutationOneSampleVsZero(...)
%
%   Trial-level observations are sign-flipped i.i.d. within each participant;
%   the statistic is mean across participants of mean trial value (per participant).
%
%   Optional name-values (after first two arguments):
%     nPerm  (default 1000)
%     Seed   nonnegative integer; sets rng before permutations (default: do not reset)

arguments
    group1 (:,1) {mustBeNumeric}
    patientList (:,1)
    args.nPerm (1,1) double {mustBePositive, mustBeInteger} = 1000
    args.Seed (1,1) double = NaN
end

nPerm = args.nPerm;
if ~isnan(args.Seed)
    validateattributes(args.Seed, {'double', 'single'}, {'scalar', 'integer', 'nonnegative'});
    rng(args.Seed, 'twister');
end

if numel(group1) ~= numel(patientList)
    error('permutationOneSampleVsZero:group1 and patientList must have the same length.');
end

rngStream = RandStream.getGlobalStream();

muObs = participantMeansById(group1, patientList);
statisticValue = mean(muObs, 'omitnan');

nullDistribution = zeros(nPerm, 1);
for p = 1:nPerm
    g1p = permuteLabelsWithinParticipants('oneSample', group1, [], patientList, rngStream);
    muP = participantMeansById(g1p, patientList);
    nullDistribution(p) = mean(muP, 'omitnan');
end

pValue = twoSidedPValue(statisticValue, nullDistribution);

if nargout > 3
    info = struct();
    info.feasibility = permutationFeasibility('oneSample', patientList, nPerm, group1);
    info.participantMeansObserved = muObs;
    [info.participantIds, ~] = participantTrialIndices(patientList);
end

end
