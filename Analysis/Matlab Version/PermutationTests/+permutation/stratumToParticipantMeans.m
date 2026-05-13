function v = stratumToParticipantMeans(muStratum, stratumParticipantIds, participantIds)
%STRATUMTOPARTICIPANTMEANS Mean of stratum values over channels, per participant (stable order).
%
%   participantIds : Ux1 from permutation.participantTrialIndices (defines order).

arguments
    muStratum (:,1) {mustBeNumeric}
    stratumParticipantIds (:,1)
    participantIds (:,1)
end

v = zeros(numel(participantIds), 1);
sp = string(stratumParticipantIds(:));
for k = 1:numel(participantIds)
    v(k) = mean(muStratum(sp == string(participantIds(k))), 'omitnan');
end

end
