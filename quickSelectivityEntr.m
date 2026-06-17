% quickSelectivityEntr.m
% new selectivity metric
% from Manny 2026/04/17

% (i.e. frequency), returns the entropy value that characterized the
% response sharpness i.e. tuning strength of the neuron
function [sel] = quickSelectivityEntr(Ys)
    y_baseline = Ys - min(Ys);
    % Making into prob dist
    if sum(y_baseline) == 0
        yprime = 0;
    else
        yprime = y_baseline ./ sum(y_baseline);
    end
    % Calculating entropy
    entropyVal = (-1 * sum(yprime .* log2(yprime), 'omitnan'));
    % normalizing entropy from 0 to 1 (helps with generalizing across studies, and experiments)
    % note the higher the value the less selective (0 = perfectly tuned)
    entropyVal = entropyVal ./ log2(length(Ys));
    sel = 1 - entropyVal;
end