function result = twoGroupPaired_trialPermutedCohortMean(group1, group2, patientList, nPerm, tail)
%TWOGROUPPAIRED_TRIALPERMUTEDCOHORTMEAN Two-group paired cohort test with trial permutations.
%
%   Participant summaries (paired trials within participant)
%   ----------------------------------------------------------
%       d_p = mean( group1_p - group2_p , 'omitnan' )   trial-wise within p.
%
%   Observed cohort statistic
%   -------------------------
%       T_obs = mean_p( d_p , 'omitnan' ).
%
%   Null (each b = 1..nPerm)
%   ------------------------
%   1) Within each participant, independently at each trial, optionally swap
%      group1 and group2 labels (trial-level swap mask).
%   2) Recompute d_p^(b) from permuted trial-level vectors.
%   3) T_b = mean_p( d_p^(b) , 'omitnan' ).
%
%   P-value: permutationPValue.
%
%   Requires group1, group2, patientList same length; trials row-matched.

    arguments
        group1 (:,1) double
        group2 (:,1) double
        patientList (:,1)
        nPerm (1,1) double {mustBeInteger, mustBePositive} = 1000
        tail (1,:) char = 'two-sided'
    end

    tail = validatestring(lower(char(tail)), {'two-sided','right','left'});

    if numel(group1) ~= numel(group2) || numel(group1) ~= numel(patientList)
        error('permutation:twoGroupPaired_trialPermutedCohortMean:LengthMismatch', ...
            'group1, group2, and patientList must have the same length.');
    end

    labels = permutation.generateTrialLevelPermutationLabels(patientList, nPerm, 'two_group_swap');

    aggObs = permutation.aggregateParticipantPairedMeanDiff(group1, group2, patientList);
    obsStat = permutation.cohortMeanStatistic(aggObs.participantDiff);

    nullDist = zeros(nPerm, 1);
    for b = 1:nPerm
        [g1p, g2p] = permutation.applyTrialLevelPermutation( ...
            group1, group2, patientList, labels, b);
        aggB = permutation.aggregateParticipantPairedMeanDiff(g1p, g2p, patientList);
        nullDist(b) = permutation.cohortMeanStatistic(aggB.participantDiff);
    end

    pValue = permutation.permutationPValue(nullDist, obsStat, tail);

    result = struct( ...
        'nullDist', nullDist, ...
        'obsStat', obsStat, ...
        'pValue', pValue, ...
        'tail', tail, ...
        'nPerm', nPerm, ...
        'testType', 'two-group-paired-trialPerm-cohortMean', ...
        'statType', 'cohortMeanOfParticipantMeanPairedDiff', ...
        'labelStruct', labels, ...
        'metadata', struct( ...
            'participantDiffObserved', aggObs.participantDiff, ...
            'patientIDs', aggObs.patientIDs, ...
            'formula', 'T = mean_p mean_trial(g1-g2|participant p)'));
end
