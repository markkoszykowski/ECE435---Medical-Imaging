function [stack] = import_thorax(zip_path, dir_name)
%IMPORT_THORAX Import images to 3D stack
    unzip(zip_path, dir_name);
    
    files = dir(dir_name);
    files = files(contains({files.name}, ".jpg"));
    
    first_img = imread(fullfile(dir_name, files(1).name));
    
    stack = zeros([size(first_img, 1:2) length(files)], class(first_img));
    for i = 1:length(files)
        file_path = fullfile(dir_name, files(i).name);
        
        rgb_img = imread(file_path);
        bw_img = rgb2gray(rgb_img);
        
        stack(:, :, i) = bw_img;
    end
end