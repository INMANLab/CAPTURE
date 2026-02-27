function fixations = CreateFixationTable(fixRuns,xS,yS,t)


nF = size(fixRuns,1);
fixStart = zeros(nF,1); fixEnd = zeros(nF,1); fixDur = zeros(nF,1);
fixX = zeros(nF,1); fixY = zeros(nF,1);
for i = 1:nF
    i0 = fixRuns(i,1); i1 = fixRuns(i,2);
    fixStart(i) = t(i0);
    fixEnd(i)   = t(i1);
    fixDur(i)   = fixEnd(i) - fixStart(i);
    fixX(i)     = xS(i);
    fixY(i)     = yS(i);
end
fixations = table(fixStart, fixEnd, fixDur, fixX, fixY, ...
    fixRuns(:,1), fixRuns(:,2), ...
    'VariableNames', {'start_t','end_t','dur_s','x_mean','y_mean','idx0','idx1'});

end