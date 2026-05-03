function labels = generatePermutationLabels(nObs, nPerm, mode)
%GENERATEPERMUTATIONLABELS Generate reproducible permutation label matrices.
%   labels = permutation.generatePermutationLabels(nObs, nPerm, mode)
%   returns an nObs x nPerm matrix used by downstream permutation routines.
%
%   mode:
%     - 'onesample_signflip': first column all +1 (identity), remaining
%       columns random +/-1 sign flips.
%     - 'paired_swap': first column all 0 (identity/no-swap), remaining
%       columns random 0/1 swap masks for paired observations.

    if nargin < 2 || isempty(nPerm)
        nPerm = 1000;
    end
    if nargin < 3
        error('permutation:generatePermutationLabels:MissingMode', ...
            'mode is required: onesample_signflip or paired_swap.');
    end

    validateattributes(nObs, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'nObs');
    validateattributes(nPerm, {'numeric'}, {'scalar','integer','>=',1}, mfilename, 'nPerm');
    if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
        error('permutation:generatePermutationLabels:InvalidModeType', ...
            'mode must be a char vector or scalar string.');
    end

    mode = validatestring(lower(char(mode)), {'onesample_signflip','paired_swap'});

    labels = zeros(nObs, nPerm);

    switch mode
        case 'onesample_signflip'
            labels(:, 1) = 1;
            if nPerm > 1
                labels(:, 2:end) = (rand(nObs, nPerm - 1) > 0.5) * 2 - 1;
            end
        case 'paired_swap'
            labels(:, 1) = 0;
            if nPerm > 1
                labels(:, 2:end) = rand(nObs, nPerm - 1) > 0.5;
            end
    end
end
