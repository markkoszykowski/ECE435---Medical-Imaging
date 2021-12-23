function [value] = lasso(z, y, A, lambda, s, waveletDomain)
% LASSO Returns value of LASSO objective function
    value = waverec2(z, s, waveletDomain);
    value = A .* fftshift(fft2(value));
    value = (1/2) * norm((value - y), 2) ^ 2 + lambda * norm(z, 1);
end