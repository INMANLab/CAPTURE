function figHandle = plotPermutationResult(result, group1, group2, plotTitle, patientList)
%PLOTPERMUTATIONRESULT Visualize permutation null and observed data summaries.
%
%   If patientList is provided (non-empty), the right panel shows participant-
%   level summaries (mean per participant, or paired mean difference) that
%   match the cohort test construction. Otherwise trial-level values are shown.

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
    if nargin < 5
        patientList = [];
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

    figHandle = figure('Color', 'w', 'Name', 'Permutation Test Visualization', ...
        'Visible', 'off');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    nullDist = result.nullDist(:);
    histogram(nullDist, 'FaceAlpha', 0.75);
    hold on;
    xline(result.obsStat, 'r-', 'LineWidth', 2);
    hold off;
    xlabel('Cohort statistic (null replicates)');
    ylabel('Count');
    title('Null Distribution');
    legend({'Null dist.', 'Observed'}, 'Location', 'best');
    grid on;

    nexttile;
    useParticipant = ~isempty(patientList) && numel(patientList) == numel(group1);

    if useParticipant && isfield(result, 'testType') && ...
            contains(result.testType, 'two-group', 'IgnoreCase', true) && ~isempty(group2)
        agg = permutation.aggregateParticipantPairedMeanDiff(group1, group2, patientList);
        P = numel(agg.patientIDs);
        boxchart(ones(P, 1), agg.participantDiff);
        set(gca, 'XTick', 1, 'XTickLabel', {'mean(g1-g2) per participant'});
        yline(0, 'k--', 'LineWidth', 1);
        ylabel('Participant paired mean difference');
        title('Participant summaries (observed)');
    elseif useParticipant && isempty(group2)
        agg = permutation.aggregateParticipantMeans(group1, patientList);
        P = numel(agg.patientIDs);
        boxchart(ones(P, 1), agg.participantMeans);
        set(gca, 'XTick', 1, 'XTickLabel', {'mean(group1) per participant'});
        yline(0, 'k--', 'LineWidth', 1);
        ylabel('Participant mean');
        title('Participant summaries (observed)');
    elseif isempty(group2)
        clean1 = group1(~isnan(group1));
        boxchart(ones(size(clean1)), clean1);
        set(gca, 'XTick', 1, 'XTickLabel', {'Group1 trials'});
        yline(0, 'k--', 'LineWidth', 1);
        scatter(1, mean(clean1, 'omitnan'), 60, 'filled', 'MarkerFaceColor', [0.2 0.6 0.2]);
        ylabel('Trial values');
        title('Observed trials');
    else
        validMask = ~isnan(group1) & ~isnan(group2);
        clean1 = group1(validMask);
        clean2 = group2(validMask);
        xVals = [ones(size(clean1)); 2 * ones(size(clean2))];
        yVals = [clean1; clean2];
        boxchart(xVals, yVals);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Group1', 'Group2'});
        hold on;
        scatter([1 2], [mean(clean1) mean(clean2)], 60, 'filled', ...
            'MarkerFaceColor', [0.2 0.6 0.2]);
        hold off;
        ylabel('Trial values');
        title('Observed trials');
    end
    grid on;

    statLabel = 'cohortMean';
    if isfield(result, 'statType')
        statLabel = char(result.statType);
    end
    statsText = sprintf(['testType: %s\nstatistic: %s\nobsStat: %.4f\n', ...
        'pValue: %.4g\nnPerm: %d'], ...
        localField(result, 'testType', 'n/a'), statLabel, ...
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
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = fallback;
    end
end
