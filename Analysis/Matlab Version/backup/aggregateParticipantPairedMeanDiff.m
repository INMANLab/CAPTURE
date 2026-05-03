function out = aggregateParticipantPairedMeanDiff(group1, group2, patientList)
%AGGREGATEPARTICIPANTPAIREDMEANDIFF Per-participant mean paired difference.
%
%   For each participant p, on trials where both group1 and group2 are finite,
%
%       d_p = mean( group1_p - group2_p , 'omitnan' )
%
%   evaluated trial-wise within that participant. If a participant has no
%   valid paired trials, an error is raised.

    arguments
        group1 (:,1) double
        group2 (:,1) double
        patientList (:,1)
    end

    if numel(group1) ~= numel(group2) || numel(group1) ~= numel(patientList)
        error('permutation:aggregateParticipantPairedMeanDiff:LengthMismatch', ...
            'group1, group2, and patientList must have the same length.');
    end

    [patientIDs, ~, idx] = unique(patientList, 'stable');
    P = numel(patientIDs);
    participantDiff = nan(P, 1);

    for k = 1:P
        mask = idx == k;
        d = group1(mask) - group2(mask);
        d = d(isfinite(d));
        if isempty(d)
            error('permutation:aggregateParticipantPairedMeanDiff:NoValidTrials', ...
                'Participant %s has no valid paired finite trials.', ...
                string(patientIDs(k)));
        end
        participantDiff(k) = mean(d);
    end

    out = struct( ...
        'patientIDs', patientIDs, ...
        'participantDiff', participantDiff, ...
        'nPatients', P);
end
