function [riskFlag,validationRisk] = FeasibilityCheck(mode, patientList, nPerm)
%FEASIBILITY Unique null draws vs requested permutations.
% Returns logical Zero if there is no risk
arguments
    mode
    patientList (:,1)
    nPerm (1,1) {mustBePositive, mustBeInteger}
end

mode = validatestring(mode, {'twoSample', 'oneSample'});

[~, indexCell] = permutation.participantTrialIndices(patientList);
U = numel(indexCell);

validationRisk = zeros(U,1);
if strcmp(mode, 'twoSample')
    for k = 1:U
        nObs = length(indexCell{k});
        nTotal = factorial(2*nObs)/factorial(nObs)^2;
        if(nTotal<=(10*nPerm))
            validationRisk(k)=1;
            warning("Permutation May not be valid, Targeted number of permutations = "+nPerm+" Number of possible permutations = "+nTotal)
        end
    end
else
    for k = 1:U
        nObs = length(indexCell{k});
        nTotal = 2^nObs;
        if(nTotal<=(10*nPerm))
            validationRisk(k)=1;
            warning("Permutation May not be valid, Targeted number of permutations = "+nPerm+" Number of possible permutations = "+nTotal)
        end
    end
end
riskFlag = logical(sum(validationRisk));
end




