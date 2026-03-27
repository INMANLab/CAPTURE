function PlotSimilarityBarsFromTable(allTbl)
% Plot grouped bar graph with:
% - mean across participants
% - SEM error bars
% - participant dots (same color as bars, no legend entry)

reqVars = {'pIdx','measure','within_X1','within_X2','between'};
for i = 1:numel(reqVars)
    if ~ismember(reqVars{i}, allTbl.Properties.VariableNames)
        error('allTbl must contain "%s".', reqVars{i});
    end
end

if iscellstr(allTbl.measure) || ischar(allTbl.measure)
    allTbl.measure = string(allTbl.measure);
end

measures = unique(allTbl.measure,'stable');
pList    = unique(allTbl.pIdx,'stable');

nM = numel(measures);
nP = numel(pList);

vals1 = nan(nP,nM);
vals2 = nan(nP,nM);
vals3 = nan(nP,nM);

for m = 1:nM
    for p = 1:nP
        idx = allTbl.measure==measures(m) & allTbl.pIdx==pList(p);
        if any(idx)
            vals1(p,m) = mean(allTbl.within_X1(idx),'omitnan');
            vals2(p,m) = mean(allTbl.within_X2(idx),'omitnan');
            vals3(p,m) = mean(allTbl.between(idx),'omitnan');
        end
    end
end

Y = [mean(vals1,1,'omitnan')' ...
     mean(vals2,1,'omitnan')' ...
     mean(vals3,1,'omitnan')'];

E = [localSEM(vals1)' ...
     localSEM(vals2)' ...
     localSEM(vals3)'];

figure;
b = bar(Y,'grouped'); hold on; box off;

% keep bar colors
c1 = b(1).FaceColor;
c2 = b(2).FaceColor;
c3 = b(3).FaceColor;

xticks(1:nM);
xticklabels(measures);
xtickangle(45);

ylabel('Similarity');
xlabel('Measure');

legend({'Within X1','Within X2','Between X1-X2'},'Location','best');

% bar centers
xPos = nan(nM,3);
for k = 1:3
    xPos(:,k) = b(k).XEndPoints(:);
end

% error bars
for k = 1:3
    errorbar(xPos(:,k),Y(:,k),E(:,k), ...
        'k','linestyle','none','linewidth',1,'capsize',8);
end

% participant dots (same color, no legend)
jit = 0.06;

for m = 1:nM

    % within1
    y = vals1(:,m);
    x = xPos(m,1) + (rand(size(y))-0.5)*2*jit;
    plot(x(~isnan(y)),y(~isnan(y)),'.', ...
        'Color',c1,'MarkerSize',14,'HandleVisibility','off');

    % within2
    y = vals2(:,m);
    x = xPos(m,2) + (rand(size(y))-0.5)*2*jit;
    plot(x(~isnan(y)),y(~isnan(y)),'.', ...
        'Color',c2,'MarkerSize',14,'HandleVisibility','off');

    % between
    y = vals3(:,m);
    x = xPos(m,3) + (rand(size(y))-0.5)*2*jit;
    plot(x(~isnan(y)),y(~isnan(y)),'.', ...
        'Color',c3,'MarkerSize',14,'HandleVisibility','off');

end

hold off
end

function s = localSEM(X)
n = sum(~isnan(X),1);
s = std(X,0,1,'omitnan')./sqrt(n);
s(n<=1) = NaN;
end