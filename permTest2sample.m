
function [pVal, obsStat, permStats] = permTest2sample(array1, array2, numPerms)
% permTest2sample - Two-sample permutation test (unequal n allowed)
%
% Inputs:
%   array1   - vector (n1 x 1 or 1 x n1)
%   array2   - vector (n2 x 1 or 1 x n2)
%   numPerms - number of permutations
%
% Outputs:
%   pVal      - two-sided permutation p-value
%   obsStat  - observed test statistic (mean difference)
%   permStats - permutation distribution of the statistic

    % Ensure column vectors
    array1 = array1(:);
    array2 = array2(:);

    n1 = numel(array1);
    n2 = numel(array2);

    % Observed statistic
    obsStat = mean(array1) - mean(array2);

    % Pool data
    pooled = [array1; array2];
    N = n1 + n2;

    permStats = zeros(numPerms,1);

    for i = 1:numPerms
        idx = randperm(N);
        grp1 = pooled(idx(1:n1));
        grp2 = pooled(idx(n1+1:end));
        permStats(i) = mean(grp1) - mean(grp2);
    end

    % Two-sided p-value
    pVal = mean(abs(permStats) >= abs(obsStat));
end