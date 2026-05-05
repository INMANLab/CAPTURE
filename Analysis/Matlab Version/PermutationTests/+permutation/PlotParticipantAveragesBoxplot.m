function h = PlotParticipantAveragesBoxplot(pValue, statisticValue, info, options)
%PLOTPARTICIPANTAVERAGESBOXPLOT Boxplot of participant averages with test summary overlay.
%
%   permutation.PlotParticipantAveragesBoxplot(pValue, statisticValue, info)
%   h = permutation.PlotParticipantAveragesBoxplot(..., Name=Value)
%
%   info must come from permutation.OneSampleVsZero or permutation.TwoSample.
%   Recognized fields:
%     participantMeansObserved      (one-sample)
%     participantMeansObservedG1    (two-sample group 1)
%     participantMeansObservedG2    (two-sample group 2)
%     nPerm                         number of permutations (recommended)

arguments
    pValue (1,1) {mustBeNumeric}
    statisticValue (1,1) {mustBeNumeric}
    info (1,1) struct
    options.GroupLabels (1,:) string = ["Group 1", "Group 2"]
    options.YLabel (1,1) string = "Participant average"
    options.Title (1,1) string = "Participant-level averages"
    options.FigureHandle = []
end

if isempty(options.FigureHandle)
    h.fig = figure('Color', 'w');
else
    h.fig = options.FigureHandle;
    figure(h.fig);
end

hasOneSample = isfield(info, 'participantMeansObserved');
hasTwoSample = isfield(info, 'participantMeansObservedG1') && isfield(info, 'participantMeansObservedG2');

if hasOneSample
    participantAvg1 = info.participantMeansObserved(:);
    participantAvg2 = zeros(0, 1);
elseif hasTwoSample
    participantAvg1 = info.participantMeansObservedG1(:);
    participantAvg2 = info.participantMeansObservedG2(:);
else
    error(['PlotParticipantAveragesBoxplot: info must include either ', ...
        'participantMeansObserved or participantMeansObservedG1/participantMeansObservedG2.']);
end

hasSecondGroup = ~isempty(participantAvg2);
if hasSecondGroup
    allData = [participantAvg1; participantAvg2];
    grp = [ones(numel(participantAvg1), 1); 2 * ones(numel(participantAvg2), 1)];
    boxplot(allData, grp, 'Labels', cellstr(options.GroupLabels(1:2)));
else
    boxplot(participantAvg1, ones(numel(participantAvg1), 1), 'Labels', {'Group 1'});
end

hold on;

% Overlay participant points with light horizontal jitter.
if hasSecondGroup
    x1 = 1 + 0.08 * (rand(numel(participantAvg1), 1) - 0.5);
    x2 = 2 + 0.08 * (rand(numel(participantAvg2), 1) - 0.5);
    scatter(x1, participantAvg1, 24, 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(x2, participantAvg2, 24, 'filled', 'MarkerFaceAlpha', 0.6);
else
    x1 = 1 + 0.08 * (rand(numel(participantAvg1), 1) - 0.5);
    scatter(x1, participantAvg1, 24, 'filled', 'MarkerFaceAlpha', 0.6);
end

grid on;
ylabel(options.YLabel);
title(options.Title);

nPerm = NaN;
if isfield(info, 'nPerm')
    nPerm = info.nPerm;
end
summaryText = buildSummaryText(pValue, nPerm, statisticValue);
if strlength(summaryText) > 0
    ax = gca;
    xl = xlim(ax);
    yl = ylim(ax);
    xText = xl(1) + 0.02 * (xl(2) - xl(1));
    yText = yl(2) - 0.04 * (yl(2) - yl(1));
    text(xText, yText, summaryText, ...
        'VerticalAlignment', 'top', ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', [0.8 0.8 0.8], ...
        'Margin', 6);
end

h.ax = gca;
hold off;

end

function txt = buildSummaryText(pValue, nPerm, statisticValue)
parts = strings(0, 1);
if ~isnan(pValue)
    parts(end + 1) = sprintf('p-value: %.4g', pValue); %#ok<AGROW>
end
if ~isnan(nPerm)
    parts(end + 1) = sprintf('nPerm: %d', round(nPerm)); %#ok<AGROW>
end
if ~isnan(statisticValue)
    parts(end + 1) = sprintf('statistic: %.4g', statisticValue); %#ok<AGROW>
end

if isempty(parts)
    txt = "";
else
    txt = strjoin(parts, newline);
end
end
