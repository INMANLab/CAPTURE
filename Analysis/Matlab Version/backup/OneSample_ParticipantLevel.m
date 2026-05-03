function result = OneSample_ParticipantLevel(patient, group1, nPerm, statType, tail, doFDR, q)
%ONESAMPLE_PARTICIPANTLEVEL Wrapper for cohort one-sample vs zero.
%   Forwards to oneSampleVsZero_trialPermutedCohortMean. statType, doFDR, and q
%   are accepted for API compatibility but ignored.

    if nargin < 3 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 5 || isempty(tail)
        tail = 'two-sided';
    end

    result = permutation.oneSampleVsZero_trialPermutedCohortMean( ...
        group1, patient, nPerm, tail);
end
