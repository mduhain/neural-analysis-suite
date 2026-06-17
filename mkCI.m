% mkCI.m
%
% make Confidence Intervals for a vector (N x 1)
%

function [ci] = mkCI(Data,alpha)

% Confidence level
if nargin == 1
    alpha = 0.05; % 95% confidence interval
end

% Calculate mean and standard error
mean_data = mean(Data,'omitnan');
std_error = std(Data,'omitnan') / sqrt(length(Data));

% Calculate the critical value from the t-distribution
t_critical = tinv(1 - alpha/2, length(Data) - 1);

% Calculate the margin of error
margin_of_error = t_critical * std_error;

% Calculate confidence intervals
ci = zeros(2,1);
ci(1) = mean_data - margin_of_error;
ci(2) = mean_data + margin_of_error;

% Display results
% fprintf('Mean: %.4f\n', mean_data);
% fprintf('95%% Confidence Interval: [%.4f, %.4f]\n', lower_bound, upper_bound);

end