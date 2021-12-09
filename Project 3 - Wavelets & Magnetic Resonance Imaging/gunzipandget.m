function [unzipList] = gunzipandget(folder)
% GUNZIPANDGET Returns list of unzipped '.gz' files from a specified
% directory
    unzipList = {dir(fullfile(folder, "*.gz")).name};
    for i = 1:length(unzipList)
        unzipList(i) = gunzip(fullfile(folder, unzipList(i)));
    end
end