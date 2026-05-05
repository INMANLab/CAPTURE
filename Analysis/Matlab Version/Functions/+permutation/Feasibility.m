function out = Feasibility(mode, patientList, nPerm, group1, group2)
%FEASIBILITY Unique null draws vs requested permutations (log scale).
%
%   out = permutation.Feasibility("twoSample", patientList, nPerm, group1, group2)
%   out = permutation.Feasibility("oneSample", patientList, nPerm, group1)
%
%   Estimates whether there are enough distinct label assignments to use
%   nPerm unique permutations without replacement. For typical data the space
%   is enormous; Monte Carlo sampling (with replacement) is indicated.

arguments
    mode
    patientList (:,1)
    nPerm (1,1) {mustBePositive, mustBeInteger}
    group1 (:,1) {mustBeNumeric}
    group2 (:,1) {mustBeNumeric} = zeros(0, 1)
end

mode = validatestring(mode, {'twoSample', 'oneSample'});

[~, indexCell] = permutation.participantTrialIndices(patientList);
U = numel(indexCell);

if strcmp(mode, 'twoSample')
    if isempty(group2)
        error("Feasibility:twoSample requires group2.");
    end
    if ~isequal(size(group1(:)), size(group2(:)))
        error("Feasibility:group1 and group2 must be the same size.");
    end
    logNum = 0;
    for k = 1:U
        idx = indexCell{k};
        nPair = numel(idx);
        if nPair == 0
            continue
        end
        v1 = group1(idx);
        v2 = group2(idx);
        if any(isnan(v1)) || any(isnan(v2))
            error("Feasibility:NaN trial values are not supported.");
        end
        n = 2 * nPair;
        m = nPair;
        logNum = logNum + logMultinomialCoefficient(n, m, n - m);
    end
else
    if ~isempty(group2)
        error("Feasibility:oneSample does not use group2.");
    end
    logNum = 0;
    for k = 1:U
        idx = indexCell{k};
        nk = numel(idx);
        if nk == 0
            continue
        end
        v = group1(idx);
        if any(isnan(v))
            error("Feasibility:NaN trial values are not supported.");
        end
        logNum = logNum + nk * log(2);
    end
end

out = struct();
out.logNumUnique = logNum;
out.useMonteCarlo = true;

logMax = log(realmax("double"));
if logNum > logMax
    out.numUnique = Inf;
    out.canDrawUniqueWithoutReplacement = true;
    out.message = sprintf( ...
        "Distinct assignments exceed floating range (log count ~= %.3g). Using Monte Carlo with nPerm=%d.", ...
        logNum, nPerm);
else
    out.numUnique = exp(logNum);
    out.canDrawUniqueWithoutReplacement = out.numUnique >= nPerm;
    if out.canDrawUniqueWithoutReplacement
        out.message = sprintf( ...
            "Approximately %.3g distinct assignments; nPerm=%d is achievable without repetition.", ...
            out.numUnique, nPerm);
    else
        out.message = sprintf( ...
            "Approximately %.3g distinct assignments; fewer than nPerm=%d. Monte Carlo still valid (sampling with replacement).", ...
            out.numUnique, nPerm);
    end
end

end

function y = logMultinomialCoefficient(n, a, b)
if a < 0 || b < 0 || a + b ~= n
    error("Feasibility:invalid multinomial indices.");
end
m = min(a, b);
y = gammaln(n + 1) - gammaln(m + 1) - gammaln(n - m + 1);
end
