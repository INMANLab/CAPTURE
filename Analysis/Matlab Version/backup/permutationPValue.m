function pVal = permutationPValue(nullDist, obsStat, tail)
%PERMUTATIONPVALUE Two-sided or one-sided permutation p-value.
%
%   Uses the usual conservative counting rule with +1 continuity:
%       two-sided: (1 + sum(|T_b| >= |T_obs|)) / (B + 1)
%       right:     (1 + sum(T_b >= T_obs)) / (B + 1)
%       left:      (1 + sum(T_b <= T_obs)) / (B + 1)
%
%   where B = numel(nullDist).

    arguments
        nullDist (:,1) double
        obsStat (1,1) double
        tail (1,:) char = 'two-sided'
    end

    tail = validatestring(lower(char(tail)), {'two-sided','right','left'});
    B = numel(nullDist);

    switch tail
        case 'two-sided'
            pVal = (sum(abs(nullDist) >= abs(obsStat)) + 1) / (B + 1);
        case 'right'
            pVal = (sum(nullDist >= obsStat) + 1) / (B + 1);
        case 'left'
            pVal = (sum(nullDist <= obsStat) + 1) / (B + 1);
    end
end
