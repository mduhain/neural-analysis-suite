%figurem.m
% minimal color corrected figure;
%
r = version('-release');
if contains(r,'2025') 
    figure('theme','light','color',[1 1 1]);
else
    figure('Color',[1,1,1]);
end