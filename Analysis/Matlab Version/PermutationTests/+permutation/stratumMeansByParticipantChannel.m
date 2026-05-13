function muStratum = stratumMeansByParticipantChannel(group, patientList, patientChannel)
%STRATUMMEANSBYPARTICIPANTCHANNEL Mean trial value per (participant, channel) stratum.

arguments
    group (:,1) {mustBeNumeric}
    patientList (:,1)
    patientChannel (:,1)
end

[~, ~, indexCell] = permutation.participantChannelStrata(patientList, patientChannel);
S = numel(indexCell);
muStratum = zeros(S, 1);
for k = 1:S
    idx = indexCell{k};
    muStratum(k) = mean(group(idx), 'omitnan');
end

end
