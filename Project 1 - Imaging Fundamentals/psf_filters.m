function [xyz_psf_filter, xy_psf_filter, xz_psf_filter, yz_psf_fitler] = psf_filters(delta_x, delta_z, pixel_spacing, slice_spacing)
%PSF_FILTERS Return PSF filters in 3D and 2D space
    FWHM_coef = sqrt(-2 * log(.5));

    std_x = delta_x / (2 * FWHM_coef);
    std_z = delta_z / (2 * FWHM_coef);

    xyz_psf = @(x, y, z) 1 / ((2*pi)^(3/2) * std_x * std_x * std_z) * exp(-1/2 * ((x / std_x).^2 + (y / std_x).^2 + (z / std_z) .^2));
    % Should be equal to one
    disp("Integral of 3D PSF: " + integral3(xyz_psf, -Inf, Inf, -Inf, Inf, -Inf, Inf));

    xy_psf = @(x, y) 1 / ((2*pi) * std_x * std_x) * exp(-1/2 * ((x / std_x).^2 + (y / std_x).^2));
    % Should be equal to one
    disp("Integral of 2D PSF ({x,y}): " + integral2(xy_psf, -Inf, Inf, -Inf, Inf));

    xz_yz_psf = @(x, z) 1 / ((2*pi) * std_x * std_z) * exp(-1/2 * ((x / std_x).^2 + (z / std_z) .^2));
    % Should be equal to one
    disp("Integral of 2D PSF ({x,z} & {y,z}): " + integral2(xz_yz_psf, -Inf, Inf, -Inf, Inf));

    xy_filt_size = ceil(2 * std_x / pixel_spacing);
    z_filt_size = ceil(2 * std_z / slice_spacing);
    
    [X, Y, Z] = meshgrid(-xy_filt_size * pixel_spacing:pixel_spacing:xy_filt_size * pixel_spacing, ...
                            -xy_filt_size * pixel_spacing:pixel_spacing:xy_filt_size * pixel_spacing, ...
                            -z_filt_size * slice_spacing:slice_spacing:z_filt_size * slice_spacing);
                        
    xyz_psf_filter = xyz_psf(X, Y, Z);
    xyz_psf_filter = xyz_psf_filter / sum(xyz_psf_filter, "all");
                
    [X, Y] = meshgrid(-xy_filt_size * pixel_spacing:pixel_spacing:xy_filt_size * pixel_spacing, ...
                        -xy_filt_size * pixel_spacing:pixel_spacing:xy_filt_size * pixel_spacing);
                    
    xy_psf_filter = xy_psf(X, Y);
    xy_psf_filter = xy_psf_filter / sum(xy_psf_filter, "all");

    [X, Z] = meshgrid(-xy_filt_size * pixel_spacing:pixel_spacing:xy_filt_size * pixel_spacing, ...
                        -z_filt_size * slice_spacing:slice_spacing:z_filt_size * slice_spacing);
                
    xz_psf_filter = xz_yz_psf(X, Z);
    xz_psf_filter = xz_psf_filter / sum(xz_psf_filter, "all");
    yz_psf_fitler = xz_psf_filter;
end