function [participantIds, indexCell] = participantTrialIndices(patientList)
%PARTICIPANTTRIALINDICES Map each trial row to a participant (stable order).
%
%   [participantIds, indexCell] = participantTrialIndices(patientList)
%
%   patientList : Nx1 trial-level participant IDs (numeric or categorical)
%   participantIds : Ux1 unique IDs in first-occurrence order
%   indexCell{k}   : linear indices of trials for participantIds(k)

arguments
    patientList (:,1)
end

if iscategorical(patientList)
    patientList = string(patientList);
end

[participantIds, ~, ic] = unique(patientList, "stable");
U = numel(participantIds);
indexCell = cell(U, 1);
for k = 1:U
    indexCell{k} = find(ic == k);
end

end
