function labelStruct = generateTrialLevelPermutationLabels(patientList, nPerm, mode)
%GENERATETRIALLEVELPERMUTATIONLABELS Unique trial-level labels per permutation.
%
%   For each participant separately, builds an (nTrials_p x nPerm) matrix.
%   Column 1 is always the identity permutation (no randomization).
%
%   Modes
%   -----
%   'one_sample_signflip'
%       Labels are +1 or -1 applied trial-wise to group1 within participant.
%       Column 1 is all +1.
%
%   'two_group_swap'
%       Labels are 0 (keep) or 1 (swap group1 and group2 at that trial).
%       Column 1 is all 0.
%
%   Traceability
%   ------------
%   patientIDs = unique(patientList,'stable'). trialRowIndices{k} lists row
%   indices into the original Nx1 vectors for participant k, in dataset row
%   order (so labels align with group1(trialRowIndices{k}), ...).

    if nargin < 2 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 3
        error('permutation:generateTrialLevelPermutationLabels:MissingMode', ...
            'mode is required: one_sample_signflip or two_group_swap.');
    end

    validateattributes(nPerm, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nPerm');
    mode = validatestring(lower(char(mode)), {'one_sample_signflip','two_group_swap'});

    [patientIDs, ~, idx] = unique(patientList, 'stable');
    P = numel(patientIDs);

    trialRowIndices = cell(P, 1);
    labelsPerPatient = cell(P, 1);

    for k = 1:P
        trialRowIndices{k} = find(idx == k);
        nT = numel(trialRowIndices{k});
        L = zeros(nT, nPerm);
        if strcmp(mode, 'one_sample_signflip')
            L(:, 1) = 1;
            if nPerm > 1
                L(:, 2:end) = (rand(nT, nPerm - 1) > 0.5) * 2 - 1;
            end
        else
            L(:, 1) = 0;
            if nPerm > 1
                L(:, 2:end) = rand(nT, nPerm - 1) > 0.5;
            end
        end
        labelsPerPatient{k} = L;
    end

    labelStruct = struct( ...
        'mode', mode, ...
        'nPerm', nPerm, ...
        'patientIDs', patientIDs, ...
        'trialRowIndices', {trialRowIndices}, ...
        'labelsPerPatient', {labelsPerPatient}, ...
        'metadata', struct( ...
            'description', 'Column 1 is identity; labels apply within participant only.', ...
            'ordering', 'unique(patientList,''stable'') matches aggregate* functions'));
end
