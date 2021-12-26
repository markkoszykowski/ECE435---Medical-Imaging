%% Project 3 ECE435 Mark Koszykowski

clc;
clear;
close all;
%% Prelim


CSI1Folder = "CSI1";
CSI2Folder = "CSI2";
CSI3Folder = "CSI3";
CSI4Folder = "CSI4";

CSI1Files = gunzipandget(CSI1Folder);
CSI2Files = gunzipandget(CSI2Folder);
CSI3Files = gunzipandget(CSI3Folder);
CSI4Files = gunzipandget(CSI4Folder);

CSI1T1Run1File = string(CSI1Files(contains(CSI1Files, "run-01") & contains(CSI1Files, "T1")));
CSI1T1Run2File = string(CSI1Files(contains(CSI1Files, "run-02") & contains(CSI1Files, "T1")));
CSI1T2File = string(CSI1Files(contains(CSI1Files, "T2")));

CSI2T1File = string(CSI2Files(contains(CSI2Files, "T1")));
CSI2T2File = string(CSI2Files(contains(CSI2Files, "T2")));

CSI3T1File = string(CSI3Files(contains(CSI3Files, "T1")));

CSI4T1File = string(CSI4Files(contains(CSI4Files, "T1")));


CSI1T1Run1 = niftiread(CSI1T1Run1File);
CSI1T1Run2 = niftiread(CSI1T1Run2File);
CSI1T2 = niftiread(CSI1T2File);

CSI2T1 = niftiread(CSI2T1File);
CSI2T2 = niftiread(CSI2T2File);

CSI3T1 = niftiread(CSI3T1File);

CSI4T1 = niftiread(CSI4T1File);

%% I


sagittalCrossSection = round(size(CSI1T1Run1, 1) / 2);
coronalCrossSection = round(size(CSI1T1Run1, 2) / 2);
axialCrossSection = round(size(CSI1T1Run1, 3) / 2);

plotcrosssections(CSI1T1Run1, sagittalCrossSection, coronalCrossSection, axialCrossSection, "T1-weighted Run 1");
plotcrosssections(CSI1T1Run2, sagittalCrossSection, coronalCrossSection, axialCrossSection, "T1-weighted Run 2");
plotcrosssections(CSI1T2, sagittalCrossSection, coronalCrossSection, axialCrossSection, "T2-weighted");

guidingImage = rot90(squeeze(CSI1T1Run1(sagittalCrossSection, :, :)));

%% III


guidingImageFFT = fftshift(fft2(guidingImage));

figure;
imagesc(abs(guidingImageFFT));
title("Guiding Image in Fourier Domain");
xlabel("{\itk_{y}}");
ylabel("{\itk_{z}}");

flattenedGuidingImageFFT = reshape(abs(guidingImageFFT), [1, numel(guidingImageFFT)]);

figure;
histogram(flattenedGuidingImageFFT, 256);
title("Guiding Image Fourier Domain Magnitude Histogram");
xlabel("|{\itF(u,v)}|");
ylabel("Count");

zeroValues = [30 50 90];

figure;
subplot(2, 2, 1);
imagesc(ifft2(ifftshift(guidingImageFFT)));
title("Guiding Image No Fourier Coefficients Zeroed Out");
xlabel("{\ity}");
ylabel("{\itz}");

for i = 1:length(zeroValues)
    subplot(2, 2, i+1);
    zeroedGuidingImageFFT = simplecompression(abs(guidingImageFFT), zeroValues(i));
    zeroedGuidingImageFFT = guidingImageFFT .* (zeroedGuidingImageFFT ~= 0);
    reconstructed = ifft2(ifftshift(zeroedGuidingImageFFT));
    imagesc(reconstructed);
    title("Guiding Image " + zeroValues(i) + "% of Fourier Coefficients Zeroed Out (MSE=" + ...
        meansquarederror(cast(guidingImage, class(reconstructed)), reconstructed) + ")");
    xlabel("{\ity}");
    ylabel("{\itz}");
end

%% IV


dwtmode("zpd", "nodisp");

load woman;
waveletBases = ["haar" "db4" "coif3"];

for i = 1:length(waveletBases)
    plotwavelet(X, 2, waveletBases(i))
end

for i = 1:length(waveletBases)
    plotwavelet(guidingImage, 2, waveletBases(i))
end

levels = 1:3;

figure;
for i = 1:length(waveletBases)
    for l = levels
        [c, ~] = wavedec2(guidingImage, l, waveletBases(i));

        subplot(length(waveletBases), length(levels), ((i - 1)*length(levels))+l)
        histogram(c, 256);
        title("Wavelet Domain Histogram (" + waveletBases(i) + " Level " + l + ")");
        xlabel("|{\itW(u,v)}|");
        ylabel("Count");
    end
end

bestDomain = "coif3";
zeroLowest = 10;

[c, s] = wavedec2(guidingImage, 3, bestDomain);
compressedC = simplecompression(c, zeroLowest);

figure;
subplot(1, 2, 1);
imagesc(waverec2(c, s, bestDomain));
title("Guiding Image No Wavelet Coefficients Zeroed Out");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(1, 2, 2);
reconstructed = waverec2(compressedC, s, bestDomain);
imagesc(reconstructed);
title(["Guiding Image " + zeroLowest + "% of Wavelet Coefficients Zeroed Out", "(" + ...
    bestDomain + ", MSE=" + meansquarederror(cast(guidingImage, class(reconstructed)), reconstructed) + ")"]);
xlabel("{\ity}");
ylabel("{\itz}");

%% V


zeroValues = 0:100;
mses = zeros([length(waveletBases) length(s)]);
bestLevel = 3;

figure;
for i = 1:length(waveletBases)
    [c, s] = wavedec2(guidingImage, bestLevel, waveletBases(i));
    for z = zeroValues
        compressedC = simplecompression(c, z);
        reconstructed = waverec2(compressedC, s, waveletBases(i));
        mses(i, z+1) = meansquarederror(cast(guidingImage, class(reconstructed)), reconstructed);
    end
    subplot(1, length(waveletBases), i);
    plot(zeroValues, mses(i, :));
    title("MSE vs. Percent of Coefficients Deleted (" + waveletBases(i) + " Level " + bestLevel + ")");
    xlabel("Percent of Coefficients Deleted");
    ylabel("MSE");
end


zeroValues = [70 80 90];
figure;
for i = 1:length(waveletBases)
    [c, s] = wavedec2(guidingImage, bestLevel, waveletBases(i));
    for z = 1:length(zeroValues)
        compressedC = simplecompression(c, zeroValues(z));
        reconstructed = waverec2(compressedC, s, waveletBases(i));
        mse = meansquarederror(cast(guidingImage, class(reconstructed)), reconstructed);

        subplot(length(zeroValues), length(waveletBases), ((z - 1)*length(waveletBases))+i)
        imagesc(reconstructed);
        title("Reconstruction of Guiding Image (" + waveletBases(i) + " Level " + bestLevel + ", s=" + zeroValues(z) + ")");
        xlabel("{\ity}");
        ylabel("{\itz}");
    end
end

figure;
for l = levels
    [c, s] = wavedec2(rot90(squeeze(CSI1T1Run2(:, coronalCrossSection, :))), l, bestDomain);
    for z = 1:length(zeroValues)
        compressedC = simplecompression(c, zeroValues(z));
        reconstructed = waverec2(compressedC, s, bestDomain);
        mse = meansquarederror(cast(rot90(squeeze(CSI1T1Run2(:, coronalCrossSection, :))), class(reconstructed)), reconstructed);

        subplot(length(levels), length(zeroValues), ((l - 1)*length(zeroValues))+z)
        imagesc(reconstructed);
        title(["Reconstruction of Coronal Slice " + coronalCrossSection  + ", T1-weighted Run 2", "(" + bestDomain + " Level " + l + ", s=" + zeroValues(z) + ")"]);
        xlabel("{\itx}");
        ylabel("{\itz}");
    end
end

bestZeroValue = 75;
figure;
for i = 1:length(waveletBases)
    for l = levels
        [c, s] = wavedec2(rot90(squeeze(CSI1T2(sagittalCrossSection, :, :))), l, waveletBases(i));
        compressedC = simplecompression(c, bestZeroValue);
        reconstructed = waverec2(compressedC, s, waveletBases(i));
        mse = meansquarederror(cast(rot90(squeeze(CSI1T2(sagittalCrossSection, :, :))), class(reconstructed)), reconstructed);

        subplot(length(levels), length(waveletBases), ((l - 1)*length(waveletBases))+i)
        imagesc(reconstructed);
        title(["Reconstruction of Sagittal Slice " + sagittalCrossSection  + ", T2-weighted", "(" + waveletBases(i) + " Level " + l + ", s=" + bestZeroValue + ")"]);
        xlabel("{\ity}");
        ylabel("{\itz}");
    end
end


zeroValues = 0:100;
maxError = zeros(size(zeroValues));
for z = zeroValues
    errors = zeros([1 7]);

    [c, s] = wavedec2(squeeze(CSI1T1Run1(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(1) = percenterror(cast(squeeze(CSI1T1Run1(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI1T1Run2(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(2) = percenterror(cast(squeeze(CSI1T1Run2(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI1T2(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(3) = percenterror(cast(squeeze(CSI1T2(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI2T1(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(4) = percenterror(cast(squeeze(CSI2T1(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI2T2(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(5) = percenterror(cast(squeeze(CSI2T2(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI3T1(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(6) = percenterror(cast(squeeze(CSI3T1(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    [c, s] = wavedec2(squeeze(CSI4T1(sagittalCrossSection, :, :)), bestLevel, bestDomain);
    compressedC = simplecompression(c, z);
    reconstructed = waverec2(compressedC, s, bestDomain);
    error(7) = percenterror(cast(squeeze(CSI4T1(sagittalCrossSection, :, :)), class(reconstructed)), reconstructed);

    maxError(z+1) = max(error);
end

figure;
plot(zeroValues, maxError, zeroValues, .1*ones(size(zeroValues)));
title("Maximum Average Absolute Error Across Dataset (" + bestDomain + " Level " + bestLevel + ")");
xlabel("Percent of Coefficients Deleted");
ylabel("Maximum Average Absolute Error");

%% VI


close all;

sensingMatrices = 1000;

lr = .1;
lambda = .001;
iterations = 100;

ps = .2:.1:.8;

allZs = zeros([size(c) sensingMatrices length(ps)]);

figure;
for p = 1:length(ps)
    [zs, s, history] = compressedsensing(guidingImage, lr, lambda, ps(p), iterations, sensingMatrices, bestDomain, bestLevel);

    allZs(:, :, :, p) = zs;

    subplot(2, 4, p);
    plot(1:iterations, history);
    title(["Average Object Function Value", ...
        "(" + bestDomain + " Level " + bestLevel + ", p=" + ps(p) + ")"]);
    xlabel("Iteration");
    ylabel("Average LASSO Value");
end

averageSamples = 5;
mses = zeros([1 length(ps)]);

figure;
subplot(2, 4, 1);
imagesc(guidingImage);
title("Groud Truth Guiding Image");
xlabel("{\ity}");
ylabel("{\itz}");

for p = 1:length(ps)
    estimations = zeros([size(guidingImage) averageSamples]);
    for i = 1:averageSamples
        estimations(:, :, i) = waverec2(allZs(:, :, randi(sensingMatrices), p), s, bestDomain);
    end

    estimation = squeeze(mean(estimations, 3));

    subplot(2, 4, p+1);
    imagesc(abs(estimation));
    title(["Compressed Sensing Reconstruction", ...
        "(" + bestDomain + " Level " + bestLevel + ", p=" + ps(p) + ")"]);
    xlabel("{\ity}");
    ylabel("{\itz}");

    mses(p) = meansquarederror(cast(guidingImage, class(real(estimation))), real(estimation));
end

figure;
plot(ps, mses);
title("MSE vs. {\itp}");
xlabel("{\itp}");
ylabel("Estimation MSE");
