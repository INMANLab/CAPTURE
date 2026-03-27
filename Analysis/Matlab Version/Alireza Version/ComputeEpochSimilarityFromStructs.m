function simOut = ComputeEpochSimilarityFromStructs(MTd1, MTd2, idx1, idx2, method)
% Compute one epoch-by-epoch similarity matrix for each common measure
% stored in two MATLAB structures.
%
% INPUTS:
%   MTd1, MTd2 : structures with fields named like d_measureName
%                Each field must be a matrix of size [samples x epochs]
%   method     : similarity method
%                'corr'   -> Pearson correlation similarity using pdist
%                'cosine' -> cosine similarity using pdist
%
% OUTPUT:
%   simOut : structure with one field per measure, each containing:
%       .S           -> similarity matrix across epochs
%       .fieldName   -> original field name
%       .nEpochs1    -> number of valid epochs from MTd1
%       .nEpochs2    -> number of valid epochs from MTd2
%       .nEpochs     -> total number of epochs in combined matrix
%
% NOTES:
%   For each common field:
%     1. Remove epochs whose column mean is NaN
%     2. Concatenate epochs from MTd1 and MTd2:
%            X = [X1 X2]
%        where X is [samples x total_epochs]
%     3. Compute epoch-to-epoch similarity
%
%   Since pdist works on rows, epochs are transposed to rows before use.

    if nargin < 3 || isempty(method)
        method = 'correlation';
    end

    f1 = fieldnames(MTd1);
    f2 = fieldnames(MTd2);

    f1 = f1(startsWith(f1, 'd_'));
    f2 = f2(startsWith(f2, 'd_'));

    commonFields = intersect(f1, f2);

    if isempty(commonFields)
        error('No common fields starting with "d_" were found in both structures.');
    end

    simOut = struct;

    for iF = 1:numel(commonFields)
        fname = commonFields{iF};

        X1 = MTd1.(fname);
        X2 = MTd2.(fname);

        X1 = X1(:,idx1);
        X2 = X2(:,idx2);

        % Remove invalid epochs
        X1 = X1(:, ~isnan(mean(X1, 1)));
        X2 = X2(:, ~isnan(mean(X2, 1)));

        if ~ismatrix(X1) || ~ismatrix(X2)
            warning('Skipping %s because one input is not a 2D matrix.', fname);
            continue;
        end

        if size(X1,1) ~= size(X2,1)
            warning('Skipping %s because sample counts differ (%d vs %d).', ...
                fname, size(X1,1), size(X2,1));
            continue;
        end



        nEp1 = size(X1, 2);
        nEp2 = size(X2, 2);

        if nEp1 == 0 || nEp2 == 0
            warning('Skipping %s because one structure has no valid epochs after NaN removal.', fname);
            continue;
        end

        % Combine epochs from both structures
        X = cat(2, X1, X2);   % [samples x total_epochs]
        Y = X.';              % [total_epochs x samples], epochs as rows


        D = pdist(Y, method);
        S = 1 - squareform(D);


        % Ensure diagonal = 1
        S(1:size(S,1)+1:end) = 1;

        % --- Index blocks ---
        idx1 = 1:nEp1;
        idx2 = nEp1+1:nEp1+nEp2;

        % --- Within X1 ---
        if nEp1 > 1
            W1 = S(idx1, idx1);
            mask1 = triu(true(nEp1),1); % exclude diagonal
            within_X1 = mean(W1(mask1), 'omitnan');
        else
            within_X1 = NaN;
        end

        % --- Within X2 ---
        if nEp2 > 1
            W2 = S(idx2, idx2);
            mask2 = triu(true(nEp2),1);
            within_X2 = mean(W2(mask2), 'omitnan');
        else
            within_X2 = NaN;
        end

        % --- Between X1 and X2 ---
        B = S(idx1, idx2);
        between_X1X2 = mean(B(:), 'omitnan');

        % --- Store ---
        simOut.(fname).S = S;
        simOut.(fname).within_X1 = within_X1;
        simOut.(fname).within_X2 = within_X2;
        simOut.(fname).between_X1X2 = between_X1X2;

        simOut.(fname).nEpochs1 = nEp1;
        simOut.(fname).nEpochs2 = nEp2;
    end
end