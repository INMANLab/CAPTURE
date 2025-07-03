function PV = calcRWAPVal(npwr,pm)

PV = nan(size(npwr));
n = size(pm,1);
parfor (k=1:size(npwr,1),8)
    pvals = sum(npwr(k,:)>=pm,1)./n;
    pvals = min(pvals,1-pvals);
    pvals(pvals<=0) = 1/n;
    pvals = 2*pvals; %correct for 2-tail test
    PV(k,:) = pvals;
end