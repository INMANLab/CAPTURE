function muByParticipant = participantMeansById(group, patientList)
%PARTICIPANTMEANSBYID Mean trial value per participant (column vector Ux1).

arguments
    group (:,1) {mustBeNumeric}
    patientList (:,1)
end

[~, indexCell] = permutation.participantTrialIndices(patientList);
U = numel(indexCell);
muByParticipant = zeros(U, 1);
for k = 1:U
    idx = indexCell{k};
    muByParticipant(k) = mean(group(idx), 'omitnan');
end

end
