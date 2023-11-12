function [percentError] = percenterror(groundTruth, array)
% PERCENTERROR Average absolute pixelwise error
    diff = abs(groundTruth - array);
    groundTruth(groundTruth == 0) = 1;
    diff = abs(diff ./ groundTruth);
    percentError = mean(diff, "all");
end