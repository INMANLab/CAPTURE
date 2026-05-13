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
%
%   When tests used patientChannel (nested strata), participantMeans* are
%   means of stratum means across channels per participant (for visualization).
%   The overlay statistic is the aggregated test value (see Aggregation in
%   OneSampleVsZero / TwoSample); boxplots remain participant-centric.
%
%   PlotPerChannelStratumPoints (default false):
%     If true and info includes stratum means and IDs, overlays one marker per
%     (participant, channel) stratum (y = stratum mean). Each participant uses
%     one marker shape from ParticipantMarkerCycle (reused for all channels of
%     that participant), e.g. two triangles for two channels, four squares for
%     four channels. X within each box group is spread by channel; requires nested info.
%   StratumMarkerSize (default 36): marker size for stratum scatter.
%   ShowParticipantLegend (default true): legend for participant / marker mapping.
%   ParticipantMarkerCycle (default ^, s, d, ...): marker string per participant
%     index mod cycle length (participant 1 -> first, participant 2 -> second, ...).

arguments
    pValue (1,1) {mustBeNumeric}
    statisticValue (1,1) {mustBeNumeric}
    info (1,1) struct
    options.GroupLabels (1,:) string = ["Group 1", "Group 2"]
    options.YLabel (1,1) string = "Participant average"
    options.Title (1,1) string = "Participant-level averages"
    options.FigureHandle = []
    options.PlotPerChannelStratumPoints (1,1) logical = false
    options.StratumMarkerSize (1,1) double {mustBePositive} = 36
    options.ShowParticipantLegend (1,1) logical = true
    options.ParticipantMarkerCycle (1, :) string = ["^", "s", "d", "v", ">", "<", "p", "h", "o", "x"]
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

useStratumScatter = options.PlotPerChannelStratumPoints;
if useStratumScatter
    if hasTwoSample
        hasStratum = isfield(info, 'stratumMeansObservedG1') && ...
            isfield(info, 'stratumMeansObservedG2') && ...
            isfield(info, 'stratumParticipantIds') && isfield(info, 'stratumChannelIds') && ...
            isfield(info, 'participantIds');
    else
        hasStratum = isfield(info, 'stratumMeansObserved') && ...
            isfield(info, 'stratumParticipantIds') && isfield(info, 'stratumChannelIds') && ...
            isfield(info, 'participantIds');
    end
    if ~hasStratum
        warning("PlotParticipantAveragesBoxplot:PlotPerChannelStratumPoints requires " + ...
            "nested test output (stratum* fields in info). Using default participant scatter.");
        useStratumScatter = false;
    end
end

if useStratumScatter
    participantIds = info.participantIds(:);
    h.legend = plotStratumPointsByParticipant( ...
        info, hasTwoSample, participantIds, options.ParticipantMarkerCycle, ...
        options.StratumMarkerSize, options.ShowParticipantLegend);
else
    h.legend = [];
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
end

grid on;
ylabel(options.YLabel);
title(options.Title);

nPerm = NaN;
if isfield(info, 'nPerm')
    nPerm = info.nPerm;
end
summaryText = buildSummaryText(pValue, nPerm, statisticValue, info);
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

function leg = plotStratumPointsByParticipant(info, hasTwoSample, participantIds, ...
    defaultMarkers, markerSize, showLegend)

leg = [];
markers = defaultMarkers;
nMark = numel(markers);
U = numel(participantIds);

if hasTwoSample
    y1 = info.stratumMeansObservedG1(:);
    y2 = info.stratumMeansObservedG2(:);
    sp = info.stratumParticipantIds(:);
    sc = info.stratumChannelIds(:);
    if numel(y1) ~= numel(sp) || numel(y2) ~= numel(sp)
        error("PlotParticipantAveragesBoxplot:stratum vectors length mismatch.");
    end
    x1 = stratumXInBox(sp, sc, participantIds, 1);
    x2 = stratumXInBox(sp, sc, participantIds, 2);
    [leg, ~] = scatterStrataByParticipant(x1, y1, sp, participantIds, markers, nMark, markerSize, showLegend, true);
    scatterStrataByParticipant(x2, y2, sp, participantIds, markers, nMark, markerSize, false, false);
else
    y = info.stratumMeansObserved(:);
    sp = info.stratumParticipantIds(:);
    sc = info.stratumChannelIds(:);
    if numel(y) ~= numel(sp)
        error("PlotParticipantAveragesBoxplot:stratum vectors length mismatch.");
    end
    x = stratumXInBox(sp, sc, participantIds, 1);
    [leg, ~] = scatterStrataByParticipant(x, y, sp, participantIds, markers, nMark, markerSize, showLegend, true);
end

end

function x = stratumXInBox(stratumParticipantIds, stratumChannelIds, participantIds, baseX)
% Horizontal spread by channel rank within participant; deterministic.
S = numel(stratumParticipantIds);
x = zeros(S, 1);
width = 0.14;

for k = 1:numel(participantIds)
    pid = participantIds(k);
    mask = participantIdMask(stratumParticipantIds, pid);
    if ~any(mask)
        continue
    end
    chanThis = stratumChannelIds(mask);
    uc = sortChannelsUnique(chanThis);
    nC = numel(uc);
    idxList = find(mask);
    for ii = 1:numel(idxList)
        s = idxList(ii);
        cid = stratumChannelIds(s);
        rank = channelRank(cid, uc);
        x(s) = baseX + (rank - (nC + 1) / 2) * width;
    end
end
end

function uc = sortChannelsUnique(chanThis)
uc = unique(chanThis, 'stable');
if isnumeric(uc)
    uc = sort(uc);
elseif isstring(uc) || ischar(uc)
    uc = sort(string(uc(:)));
else
    try
        [~, ord] = sort(string(uc(:)));
        uc = uc(ord);
    catch %#ok<CTCH>
        % keep stable unique order
    end
end
end

function r = channelRank(cid, sortedUniqueChannels)
sc = string(sortedUniqueChannels(:));
r = find(sc == string(cid), 1);
if isempty(r)
    r = 1;
end
end

function mask = participantIdMask(stratumParticipantIds, pid)
mask = string(stratumParticipantIds(:)) == string(pid);
end

function [leg, hSc] = scatterStrataByParticipant(x, y, stratumParticipantIds, participantIds, ...
    markers, nMark, markerSize, showLegend, includeInLegend)

U = numel(participantIds);
hSc = gobjects(U, 1);
for k = 1:U
    pid = participantIds(k);
    mask = participantIdMask(stratumParticipantIds, pid);
    if ~any(mask)
        continue
    end
    mk = char(markers(mod(k - 1, nMark) + 1));
    col = linesColor(k);
    if includeInLegend && showLegend
        hSc(k) = scatter(x(mask), y(mask), markerSize, ...
            'Marker', mk, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col * 0.55, ...
            'MarkerFaceAlpha', 0.85, ...
            'DisplayName', char(participantLegendLabel(pid)));
    else
        hSc(k) = scatter(x(mask), y(mask), markerSize, ...
            'Marker', mk, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col * 0.55, ...
            'MarkerFaceAlpha', 0.85, ...
            'HandleVisibility', 'off');
    end
end
nonEmpty = false(U, 1);
for k = 1:U
    nonEmpty(k) = isgraphics(hSc(k));
end
if showLegend && includeInLegend && any(nonEmpty)
    leg = legend(hSc(nonEmpty), 'Location', 'bestoutside', 'AutoUpdate', 'off');
else
    leg = [];
end
end

function c = linesColor(k)
ax = gca;
pal = ax.ColorOrder;
if isempty(pal)
    pal = lines(7);
end
c = pal(mod(k - 1, size(pal, 1)) + 1, :);
end

function lbl = participantLegendLabel(pid)
lbl = "P: " + string(pid);
end

function txt = buildSummaryText(pValue, nPerm, statisticValue, info)
arguments
    pValue (1,1) {mustBeNumeric}
    nPerm (1,1) {mustBeNumeric}
    statisticValue (1,1) {mustBeNumeric}
    info = []
end
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
if ~isempty(info) && isstruct(info) && isfield(info, 'Aggregation')
    parts(end + 1) = sprintf('Aggregation: %s', char(string(info.Aggregation))); %#ok<AGROW>
end
if ~isempty(info) && isstruct(info) && ...
        (isfield(info, 'stratumMeansObserved') || isfield(info, 'stratumMeansObservedG1'))
    parts(end + 1) = "nested: (participant x channel)"; %#ok<AGROW>
end

if isempty(parts)
    txt = "";
else
    txt = strjoin(parts, newline);
end
end
