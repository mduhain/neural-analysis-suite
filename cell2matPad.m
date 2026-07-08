% cell2matPad.m
%   cell2mat with NaN padding!
%   [locMatFixed] = cell2matPad(locMat)
%
%   Accepts 1D (1xN or Nx1) cell arrays of vectors (M x 1) of different sizes, M
%   Returns an [NxM] matrix padded with NaNs
%
%


function [locMatFixed] = cell2matPad(locMat)
    maxLength = max(cellfun(@length,locMat));
    paddedCellArray = cell(size(locMat)); 
    for ni = 1:length(locMat)
        locVector = locMat{ni};
        numNaNs = maxLength - length(locVector);
        paddedCellArray{ni} = [locVector; NaN(numNaNs, 1)];
    end
    locMatFixed = cell2mat(paddedCellArray);
end