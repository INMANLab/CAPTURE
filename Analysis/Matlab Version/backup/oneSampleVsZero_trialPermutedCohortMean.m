function result = oneSampleVsZero_trialPermutedCohortMean(group1, patientList, nPerm, tail)
%ONESAMPLEVSZERO_TRIALPERMUTEDCOHORTMEAN One-sample vs zero with trial-level permutation.
%
%   Participant summaries
%   -----------------------
%   For each participant p (unique patientList, stable order),
%       m_p = mean( trial values of group1 within p , 'omitnan' ).
%
%   Observed cohort statistic
%   -------------------------
%       T_obs = mean_p( m_p , 'omitnan' ).
%
%   Null (each b = 1..nPerm)
%   ------------------------
%   1) Draw independent trial-wise sign flips (+/-1) within each participant
%      (column b of generateTrialLevelPermutationLabels).
%   2) Apply to trial-level group1 -> group1^(b).
%   3) Recompute m_p^(b) from group1^(b).
%   4) T_b = mean_p( m_p^(b) , 'omitnan' ).
%
%   P-value
%   -------
%   See permutationPValue (two-sided or one-sided tail).
%
%   Outputs: nullDist (nPerm x 1), obsStat, pValue, metadata.

    arguments
        group1 (:,1) double
        patientList (:,1)
        nPerm (1,1) double {mustBeInteger, mustBePositive} = 1000
        tail (1,:) char = 'two-sided'
    end

    tail = validatestring(lower(char(tail)), {'two-sided','right','left'});

    if numel(group1) ~= numel(patientList)
        error('permutation:oneSampleVsZero_trialPermutedCohortMean:LengthMismatch', ...
            'group1 and patientList must have the same length.');
    end

    labels = permutation.generateTrialLevelPermutationLabels(patientList, nPerm, 'one_sample_signflip');

    aggObs = permutation.aggregateParticipantMeans(group1, patientList);
    obsStat = permutation.cohortMeanStatistic(aggObs.participantMeans);

    nullDist = zeros(nPerm, 1);
    for b = 1:nPerm
        [g1p, ~] = permutation.applyTrialLevelPermutation( ...
            group1, [], patientList, labels, b);
        aggB = permutation.aggregateParticipantMeans(g1p, patientList);
        nullDist(b) = permutation.cohortMeanStatistic(aggB.participantMeans);
    end

    pValue = permutation.permutationPValue(nullDist, obsStat, tail);

    result = struct( ...
        'nullDist', nullDist, ...
        'obsStat', obsStat, ...
        'pValue', pValue, ...
        'tail', tail, ...
        'nPerm', nPerm, ...
        'testType', 'one-sample-vs-zero-trialPerm-cohortMean', ...
        'statType', 'cohortMeanOfParticipantMeans', ...
        'labelStruct', labels, ...
        'metadata', struct( ...
            'participantMeansObserved', aggObs.participantMeans, ...
            'patientIDs', aggObs.patientIDs, ...
            'formula', 'T = mean_p mean_trial(group1|participant p)'));
end
