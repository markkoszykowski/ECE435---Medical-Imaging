function [snp_stack] = salt_and_pepper(image, noise_probs)
%SALT_AND_PEPPER Add salt and pepper noise to an image with given
%probabilities
    image = squeeze(image);
    
    salt = intmin(class(image));
    pepper = intmax(class(image));
    
    snp_stack = zeros([size(image) length(noise_probs)], class(image));
    for i = 1:length(noise_probs)
        p = noise_probs(i);
        
        rands = rand(size(image));

        image(rands <= p) = salt;
        image(rands >= (1-p)) = pepper;

        snp_stack(:, :, i) = image;
    end
end