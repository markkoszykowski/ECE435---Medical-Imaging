function [array] = simplecompression(array, s)
% SIMPLECOMPRESSION Returns the input array zeroed below a provided
% percentile
    threshold = prctile(abs(array), s, "all");
    array(abs(array) < threshold) = 0;
end