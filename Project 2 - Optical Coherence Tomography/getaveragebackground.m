function [averageBackgroundScan, smoothedBackground, nonBackgroundScans] = getaveragebackground(scans, numBackgrounds, smoothDegree)
    disp("Running getaveragebackground...");
    lineCameraDomain = 0:size(scans, 1)-1;
    
    disp("    Finding Average Background...");
    tic
    backgroundScans = scans(:, 1:numBackgrounds);
    
    averageBackgroundScan = sum(backgroundScans, 2) / numBackgrounds;
    disp("    Done! Time elapsed: " + toc + " seconds");
    
    disp("    Finding Smooth Background of Degree " + smoothDegree + "...");
    tic
    [p, ~, mu] = polyfit(lineCameraDomain, averageBackgroundScan, smoothDegree);
    
    smoothedBackground = polyval(p, lineCameraDomain, [], mu).';
    disp("    Done! Time elapsed: " + toc + " seconds");
    
    nonBackgroundScans = scans(:, numBackgrounds+1:end);
    disp(" ");
end