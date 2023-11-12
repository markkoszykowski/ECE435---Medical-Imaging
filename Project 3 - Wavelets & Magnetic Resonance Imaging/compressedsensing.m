function [zs, s, history] = compressedsensing(groundTruth, lr, lambda, p, iterations, sensingMatrices, waveletDomain, levels)
% COMPRESSEDSENSING Performs proximal gradient descent to compute
% "sensingMatrices" estimations in domain of sparsity
    f = waitbar(0, "Performing Proximal Gradient Descent...");
    start = tic;
    As = (rand([size(groundTruth) sensingMatrices]) < p);

    groundTruth = fftshift(fft2(groundTruth));
    ys = repmat(groundTruth, [1 1 sensingMatrices]);
    ys = As .* ys;

    [forSize, s] = wavedec2(groundTruth, levels, waveletDomain);
    zs = randi([intmin("int16") intmax("int16")], [size(forSize) sensingMatrices]);
    
    history = zeros([1 iterations]);

    itTimes = zeros([1 iterations]);
    waitbar(0.01, f, "Preliminary Setup Done, Starting Algorithm...");
    for i = 1:iterations
        itStart = tic;
        errors = zeros([1 sensingMatrices]);
        for k = 1:sensingMatrices
            measurementMap = waverec2(zs(:, :, k), s, waveletDomain);
            measurementMap = fftshift(fft2(measurementMap));
            measurementMap = As(:, :, k) .* measurementMap;

            measurementMap = measurementMap - ys(:, :, k);

            measurementMap = As(:, :, k) .* measurementMap;
            measurementMap = ifft2(ifftshift(measurementMap));
            measurementMap = wavedec2(measurementMap, levels, waveletDomain);

            zHat = zs(:, :, k) - lr * measurementMap;

            zs(:, :, k) = (abs(zHat - (lr * lambda)) .* exp(1j * angle(zHat))) .* (abs(zHat) > (lr * lambda));

            errors(k) = lasso(zs(:, :, k), ys(:, :, k), As(:, :, k), lambda, s, waveletDomain);
        end
        history(i) = mean(errors);
        itTimes(i) = toc(itStart);
        waitbar(i/iterations, f, ["Iterations Completed: " + i + "/" + iterations, ...
            "Runtime: " + toc(start) + "s, ETA: " + ...
            ((iterations - i) * mean(itTimes(1:i))) + "s"]);
    end
    endTime = toc(start);
    close(f);
    disp("Total Runtime: " + endTime);
end