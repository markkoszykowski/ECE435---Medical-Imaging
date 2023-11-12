%% Project 1 ECE435 Mark Koszykowski

clc;
clear;
close all;
%% I


zip_file = "thoraxCT.zip";
unzip_directory = "thorax";

stack = import_thorax(zip_file, unzip_directory);

x_slice = round(size(stack, 1) / 2);
y_slice = round(size(stack, 2) / 2);
z_slice = round(size(stack, 3) / 2);

figure;

pos = [0.1, 0.5, 0.3, 0.5];
subplot("Position", pos)
imshow(squeeze(stack(:, x_slice, :))');
title("X-Slice ({\itx}=" + x_slice + ")");
xlabel("{\ity}");
ylabel("{\itz}");

pos = [0.6, 0.5, 0.3, 0.5];
subplot("Position", pos);
imshow(squeeze(stack(y_slice, :, :))');
title("Y-Slice ({\ity}=" + (size(stack, 2) - y_slice) + ")");
xlabel("{\itx}");
ylabel("{\itz}");

pos = [0.3, 0.1, 0.4, 0.4];
subplot("Position", pos);
imshow(squeeze(stack(:, :, z_slice)));
title("Z-Slice ({\itz}=" + (size(stack, 3) - z_slice) + ")");
xlabel("{\itx}");
ylabel("{\ity}");

%% II


filter_sizes = [3 5];

[means, medians] = mean_median_filtering(stack(:, :, z_slice), filter_sizes);

figure;

for i = 1:length(filter_sizes)
    subplot(length(filter_sizes), 2, i);
    imshow(squeeze(means(:, :, i)));
    title("Mean filtering ({\itn}=" + filter_sizes(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
    
    subplot(length(filter_sizes), 2, i+2);
    imshow(squeeze(medians(:, :, i)));
    title("Median filtering ({\itn}=" + filter_sizes(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
end

p_values = [.01 .03 .05];

snp = salt_and_pepper(stack(:, :, z_slice), p_values);

snp_filter_size = filter_sizes(1);

[snp_means, snp_medians] = mean_median_filtering(snp, snp_filter_size);

figure;

for i = 1:length(p_values)
    subplot(length(p_values), 3, 3*(i-1)+1);
    imshow(squeeze(snp(:, :, i)));
    title("Salt-N-Pepa Noise ({\itp}=" + p_values(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
    
    subplot(length(p_values), 3, 3*(i-1)+2);
    imshow(squeeze(snp_means(:, :, i)));
    title("Mean Filtering ({\itp}=" + p_values(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
   
    subplot(length(p_values), 3, 3*(i-1)+3);
    imshow(squeeze(snp_medians(:, :, i)));
    title("Median Filtering ({\itp}=" + p_values(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
end

%% III


gaussian_noise_stds = [.01, .1];
gaussian_noise_stack = zeros([size(squeeze(stack(:, :, z_slice))) length(gaussian_noise_stds)], class(stack));

figure;

for i = 1:length(gaussian_noise_stds)
    subplot(1, length(gaussian_noise_stds), i);
    gaussian_noise_stack(:, :, i) = imnoise(squeeze(stack(:, :, z_slice)), "gaussian", gaussian_noise_stds(i));
    imshow(squeeze(gaussian_noise_stack(:, :, i)))
    title("Gaussian Noise ({\it\sigma}=" + gaussian_noise_stds(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
end


gaussian_filter_stds = [.75 1];
gaussian_filter_stack = zeros(size(gaussian_noise_stack), class(stack));

figure;

for i = 1:length(gaussian_filter_stds)
    subplot(1, length(gaussian_filter_stds), i);
    gaussian_filter_stack(:, :, i) = imgaussfilt(squeeze(gaussian_noise_stack(:, :, i)), gaussian_filter_stds(i));
    imshow(squeeze(gaussian_filter_stack(:, :, i)))
    title("Gaussian Noise ({\it\sigma}=" + gaussian_noise_stds(i) + ") with Gaussian Filtering (\it{\sigma}=" + gaussian_filter_stds(i) + ")");
    xlabel("{\itx}");
    ylabel("{\ity}");
end


pixel_spacing = 0.703;
slice_spacing = 0.625;

delta_x = 2;
delta_z = 1;

[xyz_psf_filter, xy_psf_filter, xz_psf_filter, yz_psf_fitler] = psf_filters(delta_x, delta_z, pixel_spacing, slice_spacing);

figure;

pos = [0.1, 0.5, 0.3, 0.5];
subplot("Position", pos)
imshow(cast(conv2(squeeze(stack(:, x_slice, :)), yz_psf_fitler, "same"), class(stack))');
title("X-Slice ({\itx}=" + x_slice + ") with PSF");
xlabel("{\ity}");
ylabel("{\itz}");

pos = [0.6, 0.5, 0.3, 0.5];
subplot("Position", pos);
imshow(cast(conv2(squeeze(stack(y_slice, :, :)), xz_psf_filter, "same"), class(stack))');
title("Y-Slice ({\ity}=" + (size(stack, 2) - y_slice) + ") with PSF");
xlabel("{\itx}");
ylabel("{\itz}");

pos = [0.3, 0.1, 0.4, 0.4];
subplot("Position", pos);
imshow(cast(conv2(squeeze(stack(:, :, z_slice)), xy_psf_filter, "same"), class(stack)));
title("Z-Slice ({\itz}=" + (size(stack, 3) - z_slice) + ") with PSF");
xlabel("{\itx}");
ylabel("{\ity}");

%% IV


y_slice = size(stack, 2) - 180;

foramina_filter = filter_sizes(1);

[~, foramina_median] = mean_median_filtering(stack(y_slice, :, :), foramina_filter);

figure;

subplot(1, 2, 1);
imshow(squeeze(stack(y_slice, :, :))');
title("Y-Slice ({\ity}=" + (size(stack, 2) - y_slice) + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(1, 2, 2);
imshow(edge(squeeze(stack(y_slice, :, :))', "Canny"));
title("Y-Slice ({\ity}=" + (size(stack, 2) - y_slice) + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");


figure;

subplot(4, 2, 1);
imshow(foramina_median');
title("Median filtering ({\itn}=" + foramina_filter + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 2);
imshow(edge(foramina_median', "Canny"));
title("Median filtering ({\itn}=" + foramina_filter + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 3);
imshow(cast(conv2(squeeze(stack(y_slice, :, :)), xz_psf_filter, "same"), class(stack))');
title("PSF System ({\it\delta x}=" + delta_x + ", {\it\delta z}=" + delta_z + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 4);
imshow(edge(cast(conv2(squeeze(stack(y_slice, :, :)), xz_psf_filter, "same"), class(stack))', "Canny"));
title("PSF System ({\it\delta x}=" + delta_x + ", {\it\delta z}=" + delta_z + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");

gaussian_image = imnoise(squeeze(stack(y_slice, :, :)), "gaussian", gaussian_noise_stds(1));

subplot(4, 2, 5);
imshow(gaussian_image');
title("Gaussian Noise ({\it\sigma}=" + gaussian_noise_stds(1) + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 6);
imshow(edge(gaussian_image', "Canny"));
title("Gaussian Noise ({\it\sigma}=" + gaussian_noise_stds(1) + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 7);
imshow(imgaussfilt(squeeze(stack(y_slice, :, :)), 3)');
title("Gaussian filtering ({\it\sigma}=" + 3 + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(4, 2, 8);
imshow(edge(imgaussfilt(squeeze(stack(y_slice, :, :)), 3)', "Canny"));
title("Gaussian filtering ({\it\sigma}=" + 3 + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");


figure;

subplot(1, 2, 1);
imshow(squeeze(stack(:, :, z_slice)));
title("Z-Slice ({\itz}=" + (size(stack, 3) - z_slice) + ")");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(1, 2, 2);
imshow(edge(squeeze(stack(:, :, z_slice)), "Canny"));
title("Z-Slice ({\itz}=" + (size(stack, 3) - z_slice) + ") Canny Edge Detector");
xlabel("{\itx}");
ylabel("{\itz}");