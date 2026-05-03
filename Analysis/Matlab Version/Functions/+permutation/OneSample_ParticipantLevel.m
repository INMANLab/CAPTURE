function result = OneSample_ParticipantLevel(patient, group1, nPerm, statType, tail, doFDR, q)
%ONESAMPLE_PARTICIPANTLEVEL Run one-sample-vs-zero test within participants.
%   result = permutation.OneSample_ParticipantLevel(patient, group1, ...)
%   is a convenience wrapper around runParticipantPermutation.

    if nargin < 3 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 4 || isempty(statType)
        statType = 'mean';
    end
    if nargin < 5 || isempty(tail)
        tail = 'two-sided';
    end
    if nargin < 6 || isempty(doFDR)
        doFDR = true;
    end
    if nargin < 7 || isempty(q)
        q = 0.05;
    end

    result = permutation.runParticipantPermutation( ...
        group1, [], patient, 'one-sample-vs-zero', nPerm, statType, tail, doFDR, q);
end
