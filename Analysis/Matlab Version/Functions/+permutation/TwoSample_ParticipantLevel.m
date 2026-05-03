function result = TwoSample_ParticipantLevel(patient, group1, group2, nPerm, statType, tail, doFDR, q)
%TWOSAMPLE_PARTICIPANTLEVEL Legacy wrapper for participant-level testing.
%   This forwards to permutation.runParticipantPermutation using paired mode.

    if nargin < 4 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 5 || isempty(statType)
        statType = 'mean';
    end
    if nargin < 6 || isempty(tail)
        tail = 'two-sided';
    end
    if nargin < 7 || isempty(doFDR)
        doFDR = true;
    end
    if nargin < 8 || isempty(q)
        q = 0.05;
    end

    result = permutation.runParticipantPermutation( ...
        group1, group2, patient, 'two-sample-paired', nPerm, statType, tail, doFDR, q);
end