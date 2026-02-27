function report = compareFixTables(fixTableA, fixTableB, time, Fs)

N = length(time);

% --- convert to masks ---
maskA = false(N,1);
maskB = false(N,1);

for k = 1:height(fixTableA)
    maskA(fixTableA.idx0(k):fixTableA.idx1(k)) = true;
end

for k = 1:height(fixTableB)
    maskB(fixTableB.idx0(k):fixTableB.idx1(k)) = true;
end

% --- sample-level confusion ---
TP = sum(maskA & maskB);
TN = sum(~maskA & ~maskB);
FP = sum(maskA & ~maskB);
FN = sum(~maskA & maskB);

precision = TP / (TP + FP + eps);
recall    = TP / (TP + FN + eps);
f1        = 2 * precision * recall / (precision + recall + eps);
acc       = (TP + TN) / N;
iou       = TP / (TP + FP + FN + eps);

% --- fixation duration totals ---
durA = sum(maskA) / Fs;
durB = sum(maskB) / Fs;
durOverlap = sum(maskA & maskB) / Fs;

% --- event-level overlap matching ---
overlapMatrix = eventOverlapMatrix(fixTableA, fixTableB);

report.sample.TP = TP;
report.sample.FP = FP;
report.sample.FN = FN;
report.sample.TN = TN;

report.metrics.precision = precision;
report.metrics.recall    = recall;
report.metrics.f1        = f1;
report.metrics.accuracy  = acc;
report.metrics.iou       = iou;

report.duration.A_sec = durA;
report.duration.B_sec = durB;
report.duration.overlap_sec = durOverlap;

report.events.overlapMatrix = overlapMatrix;

end


function M = eventOverlapMatrix(A, B)

M = zeros(height(A), height(B));

for i = 1:height(A)
    a0 = A.idx0(i); a1 = A.idx1(i);

    for j = 1:height(B)
        b0 = B.idx0(j); b1 = B.idx1(j);

        inter = max(0, min(a1,b1) - max(a0,b0) + 1);
        union = (a1-a0+1) + (b1-b0+1) - inter;

        if union > 0
            M(i,j) = inter / union;   % IoU
        end
    end
end

end

