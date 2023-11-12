function plotwavelet(array, level, waveletDomain)
% PLOTWAVELET Plots wavelet transform at two levels
    [c, s] = wavedec2(array, level, waveletDomain);

    waveletImage = [];

    figure;
    for l = level:-1:1
        [H, V, D] = detcoef2('all', c, s, l);
        expectedSize = size(array) / (2 ^ l);

        subplot(level, 4, ((l - 1)*4)+2);
        imagesc(abs(fftshift(fft2(H))));
        title(["|{\itF(u,v)}| Horizontal Detail Coefficients", "Level " + l + " (" + waveletDomain + ")"]);

        subplot(level, 4, ((l - 1)*4)+3);
        imagesc(abs(fftshift(fft2(V))));
        title(["|{\itF(u,v)}| Vertical Detail Coefficients", "Level " + l + " (" + waveletDomain + ")"]);

        subplot(level, 4, ((l - 1)*4)+4);
        imagesc(abs(fftshift(fft2(D))));
        title(["|{\itF(u,v)}| Diagonal Detail Coefficients", "Level " + l + " (" + waveletDomain + ")"]);

        offset = s(level-l+2, :) - expectedSize + 1;
        H = wcodemat(H(ceil(offset(1)/2):(end - floor(offset(2)/2)), ceil(offset(1)/2):(end - floor(offset(2)/2))), 255, "mat", 1);
        V = wcodemat(V(ceil(offset(1)/2):(end - floor(offset(2)/2)), ceil(offset(1)/2):(end - floor(offset(2)/2))), 255, "mat", 1);
        D = wcodemat(D(ceil(offset(1)/2):(end - floor(offset(2)/2)), ceil(offset(1)/2):(end - floor(offset(2)/2))), 255, "mat", 1);

        A = appcoef2(c, s, waveletDomain, l);
        
        subplot(level, 4, ((l - 1)*4)+1);
        imagesc(abs(fftshift(fft2(A))));
        title(["|{\itF(u,v)}| Approximation Coefficients", "Level " + l + " (" + waveletDomain + ")"]);

        A = wcodemat(A(ceil(offset(1)/2):(end - floor(offset(2)/2)), ceil(offset(1)/2):(end - floor(offset(2)/2))), 255, "mat", 1);
       
        if l == level
            waveletImage = [A H; V D];
        else
            waveletImage = [waveletImage H; V D];
        end
    end

    figure;
    imagesc(waveletImage);
    colormap bone(255);
    title("Wavelet Transform at " + level + " Levels (" + waveletDomain + ")");
end