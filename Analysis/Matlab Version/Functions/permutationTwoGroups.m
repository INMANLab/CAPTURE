function [nullDistribution, pValue, statisticValue, info] = permutationTwoGroups( ...
    group1, group2, patientList, args)
%PERMUTATIONTWOGROUPS Participant-level two-sample test; trial-level label shuffle null.
%
%   [nullDistribution, pValue, statisticValue] = permutationTwoGroups( ...
%       group1, group2, patientList)
%   [...] = permutationTwoGroups(group1, group2, patientList, nPerm)
%   [...] = permutationTwoGroups(group1, group2, patientList, nPerm=2000, Seed=42)
%   [..., info] = permutationTwoGroups(...)
%
%   group1, group2, patientList must be Nx1 with the same N (paired trial rows).
%   Per participant, trial values from both conditions are pooled and randomly
%   split into two sets of equal size (same row count), then participant means
%   and the mean across participants of (mean_g1 - mean_g2) form the statistic.

arguments
    group1 (:,1) {mustBeNumeric}
    group2 (:,1) {mustBeNumeric}
    patientList (:,1)
    args.nPerm (1,1) double {mustBePositive, mustBeInteger} = 1000
    args.Seed (1,1) double = NaN
end

nPerm = args.nPerm;
if ~isnan(args.Seed)
    validateattributes(args.Seed, {'double', 'single'}, {'scalar', 'integer', 'nonnegative'});
    rng(args.Seed, 'twister');
end

if numel(group1) ~= numel(group2) || numel(group1) ~= numel(patientList)
    error("permutationTwoGroups:group1, group2, and patientList must have the same length.");
end

rngStream = RandStream.getGlobalStream();

mu1 = participantMeansById(group1, patientList);
mu2 = participantMeansById(group2, patientList);
statisticValue = mean(mu1 - mu2, 'omitnan');

nullDistribution = zeros(nPerm, 1);
for p = 1:nPerm
    [g1p, g2p] = permuteLabelsWithinParticipants( ...
        'twoSample', group1, group2, patientList, rngStream);
    m1 = participantMeansById(g1p, patientList);
    m2 = participantMeansById(g2p, patientList);
    nullDistribution(p) = mean(m1 - m2, 'omitnan');
end

pValue = twoSidedPValue(statisticValue, nullDistribution);

if nargout > 3
    info = struct();
    info.feasibility = permutationFeasibility('twoSample', patientList, nPerm, group1, group2);
    info.participantMeansObservedG1 = mu1;
    info.participantMeansObservedG2 = mu2;
    [info.participantIds, ~] = participantTrialIndices(patientList);
end

end
