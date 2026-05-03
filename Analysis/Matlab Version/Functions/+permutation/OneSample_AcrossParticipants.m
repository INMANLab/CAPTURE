function result = OneSample_AcrossParticipants(patient, group1, nPerm, statType, tail, doFDR, q)
%ONESAMPLE_ACROSSPARTICIPANTS Run one-sample-vs-zero test across participants.
%   result = permutation.OneSample_AcrossParticipants(patient, group1, ...)
%   first aggregates each participant (mean within participant) and then
%   applies one-sample permutation testing across participant summaries.

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

    result = permutation.runAcrossParticipantsPermutation( ...
        group1, [], patient, 'one-sample-vs-zero', nPerm, statType, tail, doFDR, q);
end
