function AnnotateDataRange(tRanges,yOffset,headSize)

for k = 1:size(tRanges,2)
    x1 = tRanges(1,k);
    x2 = tRanges(2,k);

    % Draw main horizontal line
    plot([x1 x2], [yOffset yOffset], 'k-', 'LineWidth', 2);

    % Left arrowhead
    plot([x1 x1+headSize], [yOffset yOffset+headSize/5], 'k-', 'LineWidth', 2);
    plot([x1 x1+headSize], [yOffset yOffset-headSize/5], 'k-', 'LineWidth', 2);

    % Right arrowhead
    plot([x2 x2-headSize], [yOffset yOffset+headSize/5], 'k-', 'LineWidth', 2);
    plot([x2 x2-headSize], [yOffset yOffset-headSize/5], 'k-', 'LineWidth', 2);
end