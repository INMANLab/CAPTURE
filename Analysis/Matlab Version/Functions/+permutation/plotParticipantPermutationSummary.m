function figHandle = plotParticipantPermutationSummary(result, group1, group2, patientList, plotTitle)
%PLOTPARTICIPANTPERMUTATIONSUMMARY Visualize per-participant outcomes.
%   figHandle = permutation.plotParticipantPermutationSummary( ...
%       result, group1, group2, patientList, plotTitle)
%
%   Creates:
%     (1) Boxplot-style chart of observations grouped by participant
%     (2) Bar/scatter panel with per-participant p-values and FDR markers

    if nargin < 4
        error('permutation:plotParticipantPermutationSummary:MissingInputs', ...
            'Inputs required: result, group1, group2, patientList.');
    end
    if nargin < 5
        plotTitle = 'Participant-Level Permutation Summary';
    end

    validateattributes(group1, {'numeric'}, {'column'}, mfilename, 'group1');
    validateattributes(patientList, {'numeric','logical','char','string','categorical'}, ...
        {'column'}, mfilename, 'patientList');
    if ~isempty(group2)
        validateattributes(group2, {'numeric'}, {'column'}, mfilename, 'group2');
        if numel(group2) ~= numel(group1)
            error('permutation:plotParticipantPermutationSummary:DimensionMismatch', ...
                'group1 and group2 must match in length.');
        end
    end
    if numel(patientList) ~= numel(group1)
        error('permutation:plotParticipantPermutationSummary:PatientLengthMismatch', ...
            'patientList length must match group1 length.');
    end
    if ~isfield(result, 'perPatient') || ~isfield(result, 'patientIDs')
        error('permutation:plotParticipantPermutationSummary:InvalidResult', ...
            'result must come from runParticipantPermutation.');
    end

    perPatient = result.perPatient(:);
    patientIDs = result.patientIDs(:);
    nPatients = numel(patientIDs);

    pRaw = nan(nPatients, 1);
    pAdj = nan(nPatients, 1);
    sigFDR = false(nPatients, 1);
    for i = 1:nPatients
        pRaw(i) = perPatient(i).pValue;
        if isfield(perPatient(i), 'pValueAdj')
            pAdj(i) = perPatient(i).pValueAdj;
        end
        if isfield(perPatient(i), 'isSignificant')
            sigFDR(i) = logical(perPatient(i).isSignificant);
        end
    end

    [~, ~, idxPatient] = unique(patientList, 'stable');

    figHandle = figure('Color', 'w', 'Name', 'Participant Permutation Summary');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Panel 1: Observed values by participant.
    nexttile;
    if isempty(group2)
        valid = ~isnan(group1);
        boxchart(idxPatient(valid), group1(valid));
        ylabel('Group1 values');
        title('Group1 by Participant');
    else
        valid1 = ~isnan(group1);
        valid2 = ~isnan(group2);
        hold on;
        boxchart(idxPatient(valid1) - 0.18, group1(valid1), ...
            'BoxFaceColor', [0.2 0.4 0.8]);
        boxchart(idxPatient(valid2) + 0.18, group2(valid2), ...
            'BoxFaceColor', [0.85 0.4 0.2]);
        hold off;
        ylabel('Observed values');
        title('Group1/Group2 by Participant');
        legend({'Group1','Group2'}, 'Location', 'best');
    end
    xlabel('Participant');
    set(gca, 'XTick', 1:nPatients, 'XTickLabel', string(patientIDs));
    xtickangle(45);
    grid on;

    % Panel 2: p-values and significance markers.
    nexttile;
    b = bar(1:nPatients, pRaw, 'FaceAlpha', 0.65, 'FaceColor', [0.4 0.4 0.8]);
    b.EdgeColor = 'none';
    hold on;
    scatter(1:nPatients, pRaw, 25, 'k', 'filled');
    if any(~isnan(pAdj))
        scatter(1:nPatients, pAdj, 55, [0.1 0.6 0.1], 'LineWidth', 1.5);
    end

    yline(0.05, '--r', 'alpha=0.05', 'LineWidth', 1);
    sigIdx = find(sigFDR);
    if ~isempty(sigIdx)
        scatter(sigIdx, pRaw(sigIdx), 90, 'pentagram', ...
            'MarkerFaceColor', [0.9 0.2 0.2], 'MarkerEdgeColor', 'k');
    end
    hold off;

    xlabel('Participant');
    ylabel('p-value');
    title('Per-Participant p-values');
    set(gca, 'XTick', 1:nPatients, 'XTickLabel', string(patientIDs), 'YLim', [0 1]);
    xtickangle(45);
    grid on;
    legendItems = {'Raw p', 'Raw p (points)'};
    if any(~isnan(pAdj))
        legendItems{end+1} = 'FDR-adjusted p';
    end
    legendItems{end+1} = 'alpha=0.05';
    if ~isempty(sigIdx)
        legendItems{end+1} = 'FDR significant';
    end
    legend(legendItems, 'Location', 'best');

    sgtitle(plotTitle);
end
