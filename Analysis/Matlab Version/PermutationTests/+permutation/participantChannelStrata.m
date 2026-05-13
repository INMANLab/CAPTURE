function [stratumParticipantIds, stratumChannelIds, indexCell] = participantChannelStrata( ...
    patientList, patientChannel)
%PARTICIPANTCHANNELSTRATA Unique (participant, channel) strata and trial indices.
%
%   [stratumParticipantIds, stratumChannelIds, indexCell] = ...
%       permutation.participantChannelStrata(patientList, patientChannel)
%
%   patientList, patientChannel : Nx1, same length (trial-level IDs).
%   Strata are unique rows in first-occurrence order (stable).
%   indexCell{k} : linear indices of trials for stratum k.

arguments
    patientList (:,1)
    patientChannel (:,1)
end

if numel(patientList) ~= numel(patientChannel)
    error("participantChannelStrata:patientList and patientChannel must have the same length.");
end

pl = patientList(:);
pc = patientChannel(:);

if iscategorical(pl)
    pl = string(pl);
end
if iscategorical(pc)
    pc = string(pc);
end

T = table(pl, pc, 'VariableNames', {'participant', 'channel'});
[U, ~, ic] = unique(T, 'rows', 'stable');
stratumParticipantIds = U.participant;
stratumChannelIds = U.channel;
S = height(U);
indexCell = cell(S, 1);
for k = 1:S
    indexCell{k} = find(ic == k);
end

end
