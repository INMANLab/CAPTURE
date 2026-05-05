function h = PlotNullDistribution(nullDistribution, statisticValue, options)
%PLOTNULLDISTRIBUTION Plot null distribution and observed statistic.
%
%   permutation.PlotNullDistribution(nullDistribution, statisticValue)
%   h = permutation.PlotNullDistribution(..., Name=Value)
%
%   Name-value pairs:
%     Title        char/string (default auto from p-value if provided)
%     PValue       scalar (shown in title when provided)
%     XLabel       default "Null statistic"
%     YLabel       default "Count"
%     LineColor    default [0.85 0.33 0.10]
%     FigureHandle []

arguments
    nullDistribution (:,1) {mustBeNumeric}
    statisticValue (1,1) {mustBeNumeric}
    options.Title = ""
    options.PValue = []
    options.XLabel (1,1) string = "Null statistic"
    options.YLabel (1,1) string = "Count"
    options.LineColor (1,3) double = [0.8500 0.3250 0.0980]
    options.FigureHandle = []
end

if isempty(options.FigureHandle)
    h.fig = figure('Color', 'w');
else
    h.fig = options.FigureHandle;
    figure(h.fig);
end

nd = double(nullDistribution(:));
histogram(nd, 'FaceColor', [0.6 0.6 0.75], 'EdgeColor', 'none');
hold on;
xline(double(statisticValue), 'Color', options.LineColor, 'LineWidth', 2);
hold off;
grid on;
xlabel(options.XLabel);
ylabel(options.YLabel);

ttl = string(options.Title);
if strlength(ttl) == 0 && ~isempty(options.PValue)
    ttl = sprintf('Permutation null vs observed (p = %.4g)', options.PValue);
elseif strlength(ttl) == 0
    ttl = "Permutation null vs observed";
end
title(char(ttl));

h.ax = gca;

end
