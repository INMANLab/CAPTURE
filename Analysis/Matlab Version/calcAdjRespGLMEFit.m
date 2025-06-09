function Adjusted = calcAdjRespGLMEFit(F,GLME,Pred)
%remove fixed and random effects from response in glme model
%ModelStr = 'ccep ~ 1 + pupil2d + ntrialslist + (1|pair);

idx = ismember(GLME.PredictorNames,Pred); %predictor of interest (to keep)
x = F.(GLME.PredictorNames{find(idx,1)}); 
y = F.(GLME.ResponseName); 

rm_names = GLME.PredictorNames(~idx);
rm_vals = nan(length(x),length(rm_names));
for k=1:length(rm_names)
    rm_vals(:,k) = F.(rm_names{k});
end

nan_idx = isnan(x)|isnan(y)|any(isnan(rm_vals),2);
x(nan_idx) = [];
y(nan_idx) = [];
rm_vals(nan_idx,:) = [];

[fE,fENames,fEStats] = fixedEffects(GLME);
[rE,rENames,rEStats] = randomEffects(GLME);

%removing random effects
yadj = y;
for m=1:length(rE)
    val = rE(m);
    grp = rENames.Group{m};
    lvl = str2double(rENames.Level{m});
    idx = strcmp(rm_names,grp);
    idx2 = (lvl==rm_vals(:,idx));
    yadj(idx2) = yadj(idx2) - val;
end

% removing fixed effects and intercept
for m=1:length(fE)
    val = fE(m);
    nam = regexprep(fENames.Name{m},'_1',''); %regexprep is a temporary fix for categorical variables with 2 categories -> need to find a better fix
    if contains(nam,'Intercept')
        yadj = yadj - val;
    else
        if ~strcmp(nam,Pred)
            idx = strcmp(nam,rm_names); 
            yadj = yadj - (val.*rm_vals(:,idx));
        end
    end
end

% rELevels = cellfun(@str2double,rENames.Level);
% for m=1:length(rm_names)
%     if ~isempty(regexp(GLME.Formula.char,rm_names{m},'once'))
%         reidxs = ismember(rENames.Group,rm_names{m});
%         relvls = rELevels(reidxs);
%         re = rE(reidxs);
%         for k=1:length(relvls)
%             idx = (rm_vals(:,m)==relvls(k));
%             if any(idx)
%                 yadj(idx) = y(idx) - re(k);
%             end
%         end
%     end
% end
% 
% %removing fixed effects
% for m=1:length(rm_names)
%     if ~isempty(regexp(GLME.Formula.char,rm_names{m},'once'))
%         idx = strcmp(regexprep(fENames.Name,'_1',''),rm_names{m}); %regexprep is a temporary fix for categorical variables with 2 categories -> need to find a better fix
%         if any(idx)
%             yadj = yadj - (fE(idx).*rm_vals(:,m));
%         end
%     end
% end

warning('off','curvefit:fit:nonDoubleXData');
warning('off','curvefit:fit:noStartPoint');
% [fobj,gof] = fit(x,yadj,'poly1','Normalize','off');
[fobj,gof] = fit(x,yadj,fittype('p1*x'));
[r,p] = corrcoef(x,yadj);
[rorig,porig] = corrcoef(x,y);

Adjusted.Response = yadj;
Adjusted.ResponseName = GLME.ResponseName;
Adjusted.Predictor = x;
Adjusted.PredictorName = Pred;
Adjusted.PolyFit = fobj;
Adjusted.R = r(1,2);
Adjusted.P = p(1,2);
Adjusted.Rorig = rorig(1,2); %uncorrected correlation 
Adjusted.Porig = porig(1,2);
Adjusted.R2 = gof.rsquare;
Adjusted.AdjR2 = gof.adjrsquare;
Adjusted.NaNIdx = nan_idx;
end

% figure;
% plot(fobj,x,yadj);
% title(regexprep(evalc('fobj'),'[A-Za-z=\s]+Poly1:','','once'))
% xlabel(Pred)
% ylabel(GLME.ResponseName)