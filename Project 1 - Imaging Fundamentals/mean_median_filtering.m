function [mean_stack, median_stack] = mean_median_filtering(image_s, filter_sizes)
%MEAN_MEDIAN_FILTERING Perform mean and median filtering on image(s) for
%given filter sizes
    image_s = squeeze(image_s);
    
    mean_stack = zeros([size(image_s) length(filter_sizes)], class(image_s));
    median_stack = zeros([size(image_s) length(filter_sizes)], class(image_s));

    for i = 1:length(filter_sizes)
        filter_size = filter_sizes(i);
        
        otherdims = repmat({':'}, 1, ndims(image_s));
        
        mean_filter = (1 / filter_size^2) * ones(filter_size, filter_size);
        mean_stack(otherdims{:}, i) = imfilter(image_s, mean_filter);
        
        median_filter = [filter_size filter_size];
        try
            median_stack(otherdims{:}, i) = medfilt2(image_s, median_filter);
        catch
            median_stack(otherdims{:}, i) = medfilt3(image_s, [median_filter 1]);
        end
    end
    mean_stack = squeeze(mean_stack);
    median_stack = squeeze(median_stack);
end