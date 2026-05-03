function [g1Out, g2Out] = applyTrialLevelPermutation(group1, group2, patientList, labelStruct, permIdx)
%APPLYTRIALLEVELPERMUTATION Apply one permutation replicate at trial level.
%
%   group2 may be [] for one-sample mode. See generateTrialLevelPermutationLabels.

    if nargin < 5
        error('permutation:applyTrialLevelPermutation:MissingInputs', ...
            'Requires group1, group2, patientList, labelStruct, permIdx.');
    end
    if nargin < 2 || isempty(group2)
        group2 = [];
    end

    validateattributes(group1, {'numeric'}, {'column'}, mfilename, 'group1');
    validateattributes(patientList, {'numeric','logical','char','string','categorical'}, ...
        {'column'}, mfilename, 'patientList');
    validateattributes(permIdx, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'permIdx');

    if permIdx > labelStruct.nPerm
        error('permutation:applyTrialLevelPermutation:InvalidPermIdx', ...
            'permIdx must be <= nPerm (%d).', labelStruct.nPerm);
    end

    g1Out = group1;
    if isempty(group2)
        g2Out = [];
    else
        validateattributes(group2, {'numeric'}, {'column'}, mfilename, 'group2');
        g2Out = group2;
    end

    mode = labelStruct.mode;
    P = numel(labelStruct.patientIDs);

    for k = 1:P
        rows = labelStruct.trialRowIndices{k};
        lab = labelStruct.labelsPerPatient{k}(:, permIdx);

        switch mode
            case 'one_sample_signflip'
                g1Out(rows) = group1(rows) .* lab;

            case 'two_group_swap'
                if isempty(group2)
                    error('permutation:applyTrialLevelPermutation:NeedGroup2', ...
                        'two_group_swap requires non-empty group2.');
                end
                sw = logical(lab);
                r = rows(sw);
                tmp = g1Out(r);
                g1Out(r) = g2Out(r);
                g2Out(r) = tmp;
        end
    end
end
