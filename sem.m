function [SE] = sem(X, dim)
    % sem() Calculates the Standard Error from pop X

    % Evaluate inputs to function
    if nargin < 2 || isempty(dim)
        dim = find(size(X)~=1,1,'first');
        if isempty(dim) % check for scalar input
            dim = 1; 
        end
    end

    % Count non-NaN values along dimension
    n = sum(~isnan(X),dim);
    % standard deviation (ignoring NaNs)
    s = std(X,0,dim,'omitnan');
    % standard error
    SE = s ./ sqrt(n);

    % % Old Version 
    % stdv = std(X,'omitnan');
    % sqrtNum = sqrt(sum(~isnan(X)));
    % SE = (stdv/sqrtNum);

end