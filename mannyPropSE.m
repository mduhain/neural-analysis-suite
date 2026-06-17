function [SE] = mannyPropSE(x,n)
%mannyPropSTD : Calculate Standard error from proportional values, from MANNY
%   Inputs: x is the sub population, n is the total population.
%   Outputs: SE is the standard error following Bernouli's Princ.

SE = sqrt(((x/n)*(1-(x/n)))/n);

end