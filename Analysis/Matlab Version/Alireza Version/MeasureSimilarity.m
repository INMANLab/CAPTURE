% 
% MTd1 = rmfield(MTd1,'d_wv');
% MTd2 = rmfield(MTd2,'d_wv');
%%
% MTd1 = rmfield(MTd1,'d_wv');
% MTd2 = rmfield(MTd2,'d_wv');
Measures = fields(MTd1);
Measures = Measures(startsWith(Measures, 'd_'));
MTd1.d_All = [];
MTd2.d_All = [];
for f = 1:length(Measures)
    MTd1.d_All = cat(1,MTd1.All, MTd1.(Measures{f}));
    MTd2.d_All = cat(1,MTd2.All, MTd2.(Measures{f}));
end

pId1 = MT1.patient;
pId2 = MT2.patient;

allTbl = table;

for pIdx = 1:5
    idx1 = pId1 == pIdx;
    idx2 = pId2 == pIdx;

    simOut = ComputeEpochSimilarityFromStructs(MTd1, MTd2, idx1, idx2, 'cosine');

    flds = fieldnames(simOut);

    for f = 1:numel(flds)
        fname = flds{f};

        if ~isfield(simOut.(fname),'within_X1')
            continue
        end

        measure = erase(fname,'d_');

        T = table;
        T.pIdx        = pIdx;
        T.measure     = string(measure);
        T.within_X1   = simOut.(fname).within_X1;
        T.within_X2   = simOut.(fname).within_X2;
        T.between     = simOut.(fname).between_X1X2;
        T.nEpochs1    = simOut.(fname).nEpochs1;
        T.nEpochs2    = simOut.(fname).nEpochs2;

        allTbl = [allTbl; T];
    end
end


PlotSimilarityBarsFromTable(allTbl);
%%
% allTbl must have:
% pIdx, measure, within_X1, within_X2, between

figure;

measures = unique(allTbl.measure, 'stable');
nM = numel(measures);


% average across pIdx for each measure
Y = nan(nM, 3);

for i = 1:nM
    idx = allTbl.measure == measures(i);
    Y(i,1) = mean(allTbl.within_X1(idx), 'omitnan');
    Y(i,2) = mean(allTbl.within_X2(idx), 'omitnan');
    Y(i,3) = mean(allTbl.between(idx),   'omitnan');
end

bar(Y, 'grouped');
xticks(1:nM);
xticklabels(measures);
xtickangle(45);

ylabel('Similarity');
xlabel('Measure');
legend({'Within X1','Within X2','Between X1-X2'}, 'Location', 'best');
title('Average Similarity Across Measures');
box off;
%%
pId1 = MT1.patient;
pId2 = MT2.patient;
wId1 = MT1.walk;
wId2 = MT2.walk;

pIdx = 5;
wIdx =1;
idx1 = pId1 == pIdx; %& wId1 == wIdx;
idx2 = pId2 == pIdx; %& wId2 == wIdx;

% idx1 = true(size(pId1));
% idx2 = true(size(pId2));

Measures = fields(MTd1);
Measures = Measures(startsWith(Measures, 'd_'));
MTd1.d_All = [];
MTd2.d_All = [];
for f = 1:length(Measures)
    MTd1.d_All = cat(1,MTd1.All, MTd1.(Measures{f}));
    MTd2.d_All = cat(1,MTd2.All, MTd2.(Measures{f}));
end

simOut = ComputeEpochSimilarityFromStructs(MTd1, MTd2,idx1,idx2, 'cosine');
Measures = fields(simOut);
figure
for f = 1:length(Measures)
    subplot(3,3,f)
    imagesc(simOut.(Measures{f}).S);
    axis image;
    title(Measures{f},'Interpreter', 'none')
    hold on
    xline(simOut.(Measures{f}).nEpochs1,'r')
    yline(simOut.(Measures{f}).nEpochs1,'r')
end

PlotSimilarityBars(simOut)

% similarity matrix for measure d_alpha
S_alpha = simOut.d_alpha.S;

imagesc(S_alpha);
axis square;
colorbar;
title('Epoch similarity: d\_alpha');