function result = TwoSample_ParticipantLevel(patient, group1, group2, nPerm, statType, tail, doFDR, q)
%TWOSAMPLE_PARTICIPANTLEVEL Wrapper for cohort two-group paired test.
%   Forwards to twoGroupPaired_trialPermutedCohortMean. statType, doFDR, q
%   are accepted for API compatibility but ignored.

    if nargin < 4 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 6 || isempty(tail)
        tail = 'two-sided';
    end

    result = permutation.twoGroupPaired_trialPermutedCohortMean( ...
        group1, group2, patient, nPerm, tail);
end
