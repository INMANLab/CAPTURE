function out = aggregateParticipantMeans(group1, patientList)
%AGGREGATEPARTICIPANTMEANS Per-participant mean of trial-level observations.
%
%   For each participant id in unique(patientList,'stable'), in trial order
%   preserved by row index of group1,
%
%       m_p = mean(group1(trials of p), 'omitnan')
%
%   Participants with zero valid (finite) trials raise an error.

    arguments
        group1 (:,1) double
        patientList (:,1)
    end

    if numel(group1) ~= numel(patientList)
        error('permutation:aggregateParticipantMeans:LengthMismatch', ...
            'group1 and patientList must have the same length.');
    end

    [patientIDs, ~, idx] = unique(patientList, 'stable');
    P = numel(patientIDs);
    participantMeans = nan(P, 1);

    for k = 1:P
        mask = idx == k;
        vals = group1(mask);
        vals = vals(isfinite(vals));
        if isempty(vals)
            error('permutation:aggregateParticipantMeans:NoValidTrials', ...
                'Participant %s has no valid (finite) trials in group1.', ...
                string(patientIDs(k)));
        end
        participantMeans(k) = mean(vals);
    end

    out = struct( ...
        'patientIDs', patientIDs, ...
        'participantMeans', participantMeans, ...
        'nPatients', P);
end
