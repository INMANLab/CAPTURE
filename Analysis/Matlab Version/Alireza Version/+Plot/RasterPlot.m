function RasterPlot(ax, evenNames, dataTable, yVals, xVals, xLabel,evTimeSeries)
    
offset = 0;
    for evIdx = 1:length(evenNames)
        % evName = evenNames(evIdx);
        % evT = xVals(evTimeSeries(:,evIdx)==1);
        % plot(ax, evT,ones(1,length(evT))*evIdx,"|",LineWidth=1.5);
        evT = evTimeSeries(:,evIdx)+offset;
        offset = offset+1.1;
        plot(ax, xVals,evT);
        hold on
    end
    % for evIdx = 1:length(evenNames)
    %     evName = evenNames(evIdx);
    %     evT = xVals(strcmp(dataTable.Event,evName));
    %     plot(ax, evT,ones(1,length(evT))*evIdx,"|",LineWidth=1.5);
    %     hold on
    % end
    ylim([1,length(evenNames)+1])
    yticks(ax,yVals)
    yticklabels(ax,evenNames)
    xlabel(ax,xLabel)
end