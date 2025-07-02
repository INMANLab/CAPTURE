function PlotAverageCorrelationMatrix(evUnique,Walks,Patients,ax)
if (~exist("ax","var"))
    fig = figure;
    ax = axes('Parent', fig);
end
%---------- Load and compute Correlation Matrix
corrMat = zeros(length(evUnique),length(evUnique),length(Patients),length(Walks));
for wID = Walks
    for pID = Patients
        evTable = LoadNewData(pID,wID);
        corrMat(:,:,pID,wID) = ComputeCorrelationMatrix(evTable,evUnique);
        % FisherZ transform the rho values
        % corrMat(:,:,pID,wID) = atanh(corrMat);
    end
end

%--------- Plot Correlation Matrix
corrMat = mean(corrMat,4,"omitmissing"); % Collapse across walks
corrMat = mean(corrMat,3,"omitmissing"); % Collapse patients
corrMat(eye(size(corrMat))==1)=NaN; %remove the diagonal values
imagesc(ax,corrMat); % Display the matrix as an image
colorbar; % Add a colorbar for reference
% caxis([-1, 1]); % Set color limits for correlation values (-1 to 1)
colormap(ax,winter); % Set the colormap

% Customize the axes
xticks(ax,1:length(evUnique)); % Set x-axis ticks
yticks(ax,1:length(evUnique)); % Set y-axis ticks
xticklabels(ax,evUnique); % Label x-axis ticks with variable names
yticklabels(ax,evUnique); % Label y-axis ticks with variable names
title(ax,'Average Correlation Matrix');


end

function evTable = LoadNewData(pID,wID)
    dat = load("RWNApp_RW"+pID+"_Walk"+wID+".mat");
    evTable = dat.evnts_tbl;
end

function corrMat = ComputeCorrelationMatrix(evTable,evUnique)
    tStart = evTable.NTP(1);
    timeVals = (evTable.NTP-tStart)/60*1000;
    tArray = (0:timeVals(end))';
    tIdx = knnsearch(tArray,timeVals);

    evTimeSeries = zeros(length(tArray),length(evUnique));
    for evIdx = 1:length(evUnique)
        evName = string(evUnique{evIdx});
        evTimeSeries(tIdx,evIdx) = strcmp(evTable.Event,evName);
        evTimeSeries(:,evIdx) = conv(evTimeSeries(:,evIdx),hamming(30),'same');
    end
    corrMat = atanh(corr(evTimeSeries));
end
