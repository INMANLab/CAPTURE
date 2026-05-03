function figHandle = plotPermutationResult(result, group1, group2, plotTitle)
%PLOTPERMUTATIONRESULT Visualize permutation test and observed data.
%   figHandle = permutation.plotPermutationResult(result, group1, group2, plotTitle)
%   creates a two-panel figure:
%     (1) histogram of nullDist with observed statistic marker
%     (2) boxplot of observed data with test results annotation
%
%   Inputs
%   ------
%   result    : struct returned by permutation test functions.
%   group1    : numeric column vector.
%   group2    : numeric column vector (optional; use [] for one-sample).
%   plotTitle : custom figure title (optional).

    if nargin < 2
        error('permutation:plotPermutationResult:MissingInputs', ...
            'At least result and group1 are required.');
    end
    if nargin < 3
        group2 = [];
    end
    if nargin < 4
        plotTitle = '';
    end

    validateattributes(group1, {'numeric'}, {'column'}, mfilename, 'group1');
    if ~isempty(group2)
        validateattributes(group2, {'numeric'}, {'column'}, mfilename, 'group2');
        if numel(group2) ~= numel(group1)
            error('permutation:plotPermutationResult:DimensionMismatch', ...
                'group1 and group2 must have equal length for paired plotting.');
        end
    end

    if ~isfield(result, 'nullDist') || ~isfield(result, 'obsStat') || ~isfield(result, 'pValue')
        error('permutation:plotPermutationResult:InvalidResult', ...
            'result must contain nullDist, obsStat, and pValue fields.');
    end

    figHandle = figure('Color', 'w', 'Name', 'Permutation Test Visualization');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel 1: Null distribution and observed statistic.
    nexttile;
    nullDist = result.nullDist(:);
    histogram(nullDist, 'FaceAlpha', 0.75);
    hold on;
    xline(result.obsStat, 'r-', 'LineWidth', 2);
    hold off;
    xlabel('Statistic value');
    ylabel('Count');
    title('Null Distribution');
    legend({'Null dist.', 'Observed stat'}, 'Location', 'best');
    grid on;

    % Panel 2: Boxplot-style summary of observations.
    nexttile;
    if isempty(group2)
        clean1 = group1(~isnan(group1));
        boxchart(ones(size(clean1)), clean1);
        set(gca, 'XTick', 1, 'XTickLabel', {'Group1'});
        yline(0, 'k--', 'LineWidth', 1);
        scatter(1, mean(clean1, 'omitnan'), 60, 'filled', 'MarkerFaceColor', [0.2 0.6 0.2]);
    else
        validMask = ~isnan(group1) & ~isnan(group2);
        clean1 = group1(validMask);
        clean2 = group2(validMask);
        xVals = [ones(size(clean1)); 2 * ones(size(clean2))];
        yVals = [clean1; clean2];
        boxchart(xVals, yVals);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Group1', 'Group2'});
        hold on;
        % Overlay group means for quick interpretation.
        scatter([1 2], [mean(clean1) mean(clean2)], 60, 'filled', ...
            'MarkerFaceColor', [0.2 0.6 0.2]);
        hold off;
    end
    ylabel('Observed values');
    title('Observed Data');
    grid on;

    % Superimposed test summary annotation.
    statsText = sprintf(['testType: %s\nstatType: %s\nobsStat: %.4f\n', ...
        'pValue: %.4g\nnPerm: %d'], ...
        localField(result, 'testType', 'n/a'), ...
        localField(result, 'statType', 'n/a'), ...
        result.obsStat, result.pValue, localField(result, 'nPerm', NaN));

    if isfield(result, 'pValueAdj') && ~isnan(result.pValueAdj)
        statsText = sprintf('%s\npValueAdj: %.4g', statsText, result.pValueAdj);
    end

    annotation(figHandle, 'textbox', [0.69 0.72 0.28 0.23], ...
        'String', statsText, 'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', 'EdgeColor', [0.6 0.6 0.6]);

    if isempty(plotTitle)
        sgtitle('Permutation Test Result');
    else
        sgtitle(plotTitle);
    end
end

function value = localField(s, fieldName, fallback)
%LOCALFIELD Return struct field value or fallback.
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = fallback;
    end
end
