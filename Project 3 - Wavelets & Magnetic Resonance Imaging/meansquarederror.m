function [mse] = meansquarederror(groundTruth, array)
% MEANSQUAREDERROR Mean Squared Error of two arrays
    diff = groundTruth - array;
    diff = diff .^ 2;
    mse = mean(diff, "all");
end