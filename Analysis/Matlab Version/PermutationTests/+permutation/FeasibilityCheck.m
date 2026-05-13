function [riskFlag, validationRisk] = FeasibilityCheck(mode, patientList, nPerm, patientChannel)
%FEASIBILITYCHECK Per-stratum check: warn if possible null draws are scarce vs nPerm.
%
%   [riskFlag, validationRisk] = permutation.FeasibilityCheck(mode, patientList, nPerm)
%   [riskFlag, validationRisk] = permutation.FeasibilityCheck(mode, patientList, nPerm, patientChannel)
%
%   validationRisk is Ux1 (one per participant, or one per (participant, channel) stratum).

arguments
    mode
    patientList (:,1)
    nPerm (1,1) {mustBePositive, mustBeInteger}
    patientChannel (:,1) = []
end

mode = validatestring(mode, {'twoSample', 'oneSample'});

if isempty(patientChannel)
    [~, indexCell] = permutation.participantTrialIndices(patientList);
else
    if numel(patientChannel) ~= numel(patientList)
        error("FeasibilityCheck:patientChannel must have the same length as patientList.");
    end
    [~, ~, indexCell] = permutation.participantChannelStrata(patientList, patientChannel);
end
U = numel(indexCell);

validationRisk = zeros(U, 1);
logThreshold = log(10) + log(double(nPerm));

if strcmp(mode, 'twoSample')
    for k = 1:U
        nObs = numel(indexCell{k});
        if nObs == 0
            continue
        end
        logNumAssignments = gammaln(2 * nObs + 1) - 2 * gammaln(nObs + 1);
        if logNumAssignments <= logThreshold
            validationRisk(k) = 1;
            if logNumAssignments < log(realmax("double"))
                nTotal = round(exp(logNumAssignments));
            else
                nTotal = Inf;
            end
            warning("Permutation May not be valid, Targeted number of permutations = %d; " + ...
                "log count of possible label splits in this stratum = %.3g (approx count = %s).", ...
                nPerm, logNumAssignments, char(string(nTotal)));
        end
    end
else
    for k = 1:U
        nObs = numel(indexCell{k});
        if nObs == 0
            continue
        end
        logNumAssignments = nObs * log(2);
        if logNumAssignments <= logThreshold
            validationRisk(k) = 1;
            if logNumAssignments < log(realmax("double"))
                nTotal = round(exp(logNumAssignments));
            else
                nTotal = Inf;
            end
            warning("Permutation May not be valid, Targeted number of permutations = %d; " + ...
                "log count of possible sign patterns in this stratum = %.3g (approx count = %s).", ...
                nPerm, logNumAssignments, char(string(nTotal)));
        end
    end
end
riskFlag = logical(sum(validationRisk));
end
