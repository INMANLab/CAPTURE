function p = twoSidedPValue(statisticValue, nullDistribution)
%TWOSIDEDPVALUE Liberal two-sided Monte Carlo p-value: (1 + #{|T*| >= |T|}) / (B + 1).

nullDistribution = nullDistribution(:);
B = numel(nullDistribution);
absT = abs(statisticValue);
absNull = abs(nullDistribution);
count = sum(absNull >= absT);
p = (1 + count) / (B + 1);

end
