function PlotSimilarityBars(simOut)
% Plot bar graph of similarity scores across measures
% Bars:
%   within_X1 , within_X2 , between_X1X2
%
% INPUT:
%   simOut : output from ComputeEpochSimilarityFromStructs

fields = fieldnames(simOut);
n = numel(fields);

within1  = nan(n,1);
within2  = nan(n,1);
between  = nan(n,1);

labels = strings(n,1);

for i = 1:n
    f = fields{i};
    
    if ~isfield(simOut.(f),'within_X1')
        continue
    end
    
    within1(i) = simOut.(f).within_X1;
    within2(i) = simOut.(f).within_X2;
    between(i) = simOut.(f).between_X1X2;
    
    % remove d_ prefix for display
    labels(i) = erase(f,"d_");
end

Y = [within1 within2 between];

figure;
bar(Y,'grouped');
box off;

xticks(1:n);
xticklabels(labels);
xtickangle(45);

ylabel('Similarity');
xlabel('Measure');

legend({'Within X1','Within X2','Between X1-X2'}, ...
       'Location','best');

title('Epoch Similarity Across Measures');

end