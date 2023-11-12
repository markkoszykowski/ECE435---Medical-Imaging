function [spatialScans] = tospatial(scans, averageBackground, smoothBackground, hammingWindow, L2K, subtract, deconvolve)
    disp("Running tospatial...");
    scanCount = size(scans, 2);
    
    backgroundRepeated = repmat(averageBackground, 1, scanCount);
    smoothBackgroundRepeated = repmat(smoothBackground, 1, scanCount);
    windowRepeated = repmat(hammingWindow, 1, scanCount);
    
    if subtract
        disp("    Subtracting Background...");
        tic
        scans = scans - backgroundRepeated;
        disp("    Done! Time elapsed: " + toc + " seconds");
    end
    
    disp("    Applying Window...");
    tic
    scans = scans .* windowRepeated;
    disp("    Done! Time elapsed: " + toc + " seconds");
    
    if deconvolve
        disp("    Applying Deconvolution...");
        tic
        scans = scans ./ smoothBackgroundRepeated;
        disp("    Done! Time elapsed: " + toc + " seconds");
    end
    
    disp("    Converting to k domain...");
    tic
    scans = L2K * scans;
    disp("    Done! Time elapsed: " + toc + " seconds");
    
    disp("    Taking FFT...");
    tic
    spatialScans = fft(scans);
    spatialScans = spatialScans(1:size(scans, 1)/2, :);
    disp("    Done! Time elapsed: " + toc + " seconds");
   
    disp(" ");
end