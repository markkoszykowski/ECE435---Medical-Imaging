%% Project 2 ECE435 Mark Koszykowski

clc;
clear;
close all;
%% Prelim


n = 1;
lambda0 = 1310;
deltaLambda = 100;
NA = 0.055;
BScanWidth = 1e6;
pixelHeight = 2700;

lateralResolution = 0.37 * (lambda0 / NA);
axialResolution = ((2*log(2))/pi) * lambda0^2 / deltaLambda;

disp("Lateral Resolution: " + lateralResolution);
disp("Axial Resolution: " + axialResolution);
disp(" ");

zipFile = "SDPM Project.zip";
unzipDirectory = "SDPM Project";

L2KFile = "L2K.mat";
BScanFile = "BScan_Layers.raw";
MScan1File = "MScan1.raw";
MScan40File = "MScan40.raw";

lineCameraPixels = 2048;
smoothDegree = 12;

BScanBackgroundCount = 175;
MScanBackgroundCount = 320;

fs = 97656.25;

unzip(zipFile, ".");

L2K = load(fullfile(unzipDirectory, L2KFile), "-mat").L2K;
BScan = readscan(fullfile(unzipDirectory, BScanFile), lineCameraPixels);
MScan1 = readscan(fullfile(unzipDirectory, MScan1File), lineCameraPixels);
MScan40 = readscan(fullfile(unzipDirectory, MScan40File), lineCameraPixels);

hammingWindow = hamming(lineCameraPixels);

%% II


[BScanBackground, BScanSmoothBackground, nonBackgroundBScan] = getaveragebackground(BScan, BScanBackgroundCount, smoothDegree);
[MScan1Background, MScan1SmoothBackground, nonBackgroundMScan1] = getaveragebackground(MScan1, MScanBackgroundCount, smoothDegree);
[MScan40Background, MScan40SmoothBackground, nonBackgroundMScan40] = getaveragebackground(MScan40, MScanBackgroundCount, smoothDegree);

figure;
subplot(1, 3, 1);
plot(BScanBackground);
hold on;
plot(BScanSmoothBackground);
title("BScan Background");
legend(["Average", " Average Smoothed"]);
xlabel("\it{\lambda}");
ylabel("Electron Count");
xlim([0 length(BScanBackground)-1]);
ylim([0 1.25*max(BScanBackground)]);
xticks([]);

subplot(1, 3, 2);
plot(MScan1Background);
hold on;
plot(MScan1SmoothBackground);
title("MScan1 Background");
legend(["Average", "Average Smoothed"]);
xlabel("\it{\lambda}");
ylabel("Electron Count");
xlim([0 length(MScan1Background)-1]);
ylim([0 1.25*max(MScan1Background)]);
xticks([]);

subplot(1, 3, 3);
plot(MScan40Background);
hold on;
plot(MScan40SmoothBackground);
title("MScan40 Background");
legend(["Average", "Average Smoothed"]);
xlabel("\it{\lambda}");
ylabel("Electron Count");
xlim([0 length(MScan40Background)-1]);
ylim([0 1.25*max(MScan40Background)]);
xticks([]);


AScanIndices = [2500, 7500];
AScan1_40Index = 5000;

AScansSpatial = abs(tospatial(nonBackgroundBScan(:, AScanIndices), BScanBackground, ...
                                BScanSmoothBackground, hammingWindow, L2K, true, true));
                        
AScan1Spatial = abs(tospatial(nonBackgroundMScan1(:, AScan1_40Index), MScan1Background, ...
                                MScan1SmoothBackground, hammingWindow, L2K, true, true));

AScan40Spatial = abs(tospatial(nonBackgroundMScan40(:, AScan1_40Index), MScan40Background, ...
                                MScan40SmoothBackground, hammingWindow, L2K, true, true));

figure;
subplot(2, 2, 1);
plot(mag2db(AScansSpatial(:, 1)));
title("A Scan of BScan\_Layers at Index = " + AScanIndices(1));
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScansSpatial(:, 1))-1]);
xticks([]);

subplot(2, 2, 2);
plot(mag2db(AScansSpatial(:, 2)));
title("A Scan of BScan\_Layers at Index = " + AScanIndices(2));
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScansSpatial(:, 2))-1]);
xticks([]);

subplot(2, 2, 3);
plot(mag2db(AScan1Spatial));
title("A Scan of MScan1 at Index = " + AScan1_40Index);
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScan1Spatial)-1]);
xticks([]);

subplot(2, 2, 4);
plot(mag2db(AScan40Spatial));
title("A Scan of MScan40 at Index = " + AScan1_40Index);
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScan40Spatial)-1]);
xticks([]);

AScansSpatialNoSubNoDeconv = abs(tospatial(nonBackgroundBScan(:, AScanIndices(1)), BScanBackground, ...
                                            BScanSmoothBackground, hammingWindow, L2K, false, false));

AScansSpatialNoDeconv = abs(tospatial(nonBackgroundBScan(:, AScanIndices(1)), BScanBackground, ...
                                        BScanSmoothBackground, hammingWindow, L2K, true, false));

figure;
subplot(1, 2, 1);
plot(mag2db(AScansSpatialNoSubNoDeconv));
title("A Scan of BScan\_Layers at Index = " + AScanIndices(1) + " (No Subtraction, No Deconvolution)");
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScansSpatialNoSubNoDeconv)-1]);
xticks([]);

subplot(1, 2, 2);
plot(mag2db(AScansSpatialNoDeconv));
title("A Scan of BScan\_Layers at Index = " + AScanIndices(1) + " (No Deconvolution)");
xlabel("\it{z} ({\it\mum})");
ylabel("Power (\it{dB})");
xlim([0 length(AScansSpatialNoDeconv)-1]);
xticks([]);

%% III


BScanSpatial = abs(tospatial(nonBackgroundBScan, BScanBackground, ...
                                BScanSmoothBackground, hammingWindow, L2K, true, true));

axialAxis = pixelHeight/2:pixelHeight:size(BScanSpatial, 1)*pixelHeight - pixelHeight/2;

lateralAxis = (BScanWidth/size(nonBackgroundBScan, 2))/2: ...
                    BScanWidth/size(nonBackgroundBScan, 2): ...
                    size(nonBackgroundBScan, 2)*BScanWidth/size(nonBackgroundBScan, 2) - ...
                    (BScanWidth/size(nonBackgroundBScan, 2))/2;

[X, Y] = meshgrid(lateralAxis, axialAxis);

trueLateralAxis = pixelHeight/2:pixelHeight:pixelHeight*floor(BScanWidth/pixelHeight) - pixelHeight/2;

[Xq, Yq] = meshgrid(trueLateralAxis, axialAxis);

BScanTrueSize = interp2(X, Y, abs(BScanSpatial), Xq, Yq).';

BScanTrueSize = (BScanTrueSize - min(BScanTrueSize, [], "all")) ./ ...
                (max(BScanTrueSize, [], "all") - min(BScanTrueSize, [], "all"));

figure;
subplot(2, 1, 1);
imshow(BScanTrueSize);
title("Unprocessed Image");
xlabel("Axial Axis");
ylabel("Lateral Axis");

subplot(2, 1, 2);
imshow(imadjust(BScanTrueSize));
title("Processed Image");
xlabel("Axial Axis");
ylabel("Lateral Axis");

BScanSpatialNoConv = abs(tospatial(nonBackgroundBScan, BScanBackground, ...
                                    BScanSmoothBackground, hammingWindow, L2K, true, false));

BScanNoConvTrueSize = interp2(X, Y, abs(BScanSpatialNoConv), Xq, Yq).';

BScanNoConvTrueSize = (BScanNoConvTrueSize - min(BScanNoConvTrueSize, [], "all")) ./ ...
                        (max(BScanNoConvTrueSize, [], "all") - min(BScanNoConvTrueSize, [], "all"));

figure;
imshow(imadjust(BScanNoConvTrueSize));
title("Processed Image (No Convolution)");
xlabel("Axial Axis");
ylabel("Lateral Axis");

%% IV


skipFirstNZ = 10;
skipFirstNDC = 1000;

MScan1Spatial = tospatial(nonBackgroundMScan1, MScan1Background, ...
                            MScan1SmoothBackground, hammingWindow, L2K, true, true);

averageAScan = sum(abs(MScan1Spatial), 2) / size(MScan1Spatial, 2);

[~, speakerIndices] = maxk(averageAScan(skipFirstNZ:end), 2);
speakerIndices = sort(speakerIndices + (skipFirstNZ-1));
disp("Indices of Speaker: " + speakerIndices);

figure;
plot(mag2db(averageAScan));
hold on
plot(speakerIndices, mag2db(averageAScan(speakerIndices)), "r*");
title("Average A Scan of MScan1");
xlabel("\it{z} ({\it\mum})");
ylabel("Power ({\itdB})");
xlim([0 length(averageAScan)-1]);
xticks([]);

[MScan1Motion1, ~] = gettimeseries(MScan1Spatial, 1/fs, speakerIndices(1), lambda0, n);
[MScan1Motion2, time] = gettimeseries(MScan1Spatial, 1/fs, speakerIndices(2), lambda0, n);

f = fs*(skipFirstNDC:((size(MScan1Spatial, 2)/2)-1))/size(MScan1Spatial, 2);

MScan1Freq1 = fft(MScan1Motion1);
MScan1Freq2 = fft(MScan1Motion2);

MScan1Freq1 = mag2db(abs(MScan1Freq1(skipFirstNDC+1:length(f)+skipFirstNDC)));
MScan1Freq2 = mag2db(abs(MScan1Freq2(skipFirstNDC+1:length(f)+skipFirstNDC)));

figure;
subplot(2, 1, 1);
plot(time, MScan1Motion1);
title("Motion at pixel " + speakerIndices(1) + " of MScan1");
xlabel("\it{t} ({\its})");
ylabel("Displacement ({\itnm})");
xlim([0 time(end)]);

subplot(2, 1, 2);
plot(time, MScan1Motion2);
title("Motion at pixel " + speakerIndices(2) + " of MScan1");
xlabel("\it{t} ({\its})");
ylabel("Displacement ({\itnm})");
xlim([0 time(end)]);

[~, freq1] = max(MScan1Freq1);
[~, freq2] = max(MScan1Freq2);

disp("Frequency at pixel " + speakerIndices(1) + " of MScan1: " + f(freq1) + " Hz");

disp("Frequency at pixel " + speakerIndices(2) + " of MScan1: " + f(freq2) + " Hz");
disp(" ");

figure;
subplot(2, 1, 1);
plot(f, MScan1Freq1, f(freq1), MScan1Freq1(freq1), "r*");
title("Frequency at pixel " + speakerIndices(1) + " of MScan1");
xlabel("\it{f} ({\itHz})");
ylabel("|D(\it{f})| ({\itdB})");
xlim([f(skipFirstNDC) f(end)]);

subplot(2, 1, 2);
plot(f, abs(MScan1Freq2), f(freq2), MScan1Freq2(freq2), "r*");
title("Frequency at pixel " + speakerIndices(2) + " of MScan1");
xlabel("\it{f} ({\itHz})");
ylabel("|D(\it{f})| ({\itdB})");
xlim([f(skipFirstNDC) f(end)]);


MScan40Spatial = tospatial(nonBackgroundMScan40, MScan40Background, ...
                            MScan40SmoothBackground, hammingWindow, L2K, true, true);

averageAScan = sum(abs(MScan40Spatial), 2) / size(MScan40Spatial, 2);

[~, speakerIndices] = maxk(averageAScan(skipFirstNZ:end), 2);
speakerIndices = sort(speakerIndices + (skipFirstNZ-1));
disp("Indices of Speaker: " + speakerIndices);

figure;
plot(mag2db(averageAScan));
hold on
plot(speakerIndices, mag2db(averageAScan(speakerIndices)), "r*");
title("Average A Scan of MScan40");
xlabel("\it{z} ({\it\mum})");
ylabel("Power ({\itdB})");
xlim([0 length(averageAScan)-1]);
xticks([]);

[MScan40Motion1, ~] = gettimeseries(MScan40Spatial, 1/fs, speakerIndices(1), lambda0, n);
[MScan40Motion2, time] = gettimeseries(MScan40Spatial, 1/fs, speakerIndices(2), lambda0, n);

f = fs*(skipFirstNDC:((size(MScan40Spatial, 2)/2)-1))/size(MScan40Spatial, 2);

MScan40Freq1 = fft(MScan40Motion1);
MScan40Freq2 = fft(MScan40Motion2);

MScan40Freq1Mag = mag2db(abs(MScan40Freq1(skipFirstNDC+1:length(f)+skipFirstNDC)));
MScan40Freq2Mag = mag2db(abs(MScan40Freq2(skipFirstNDC+1:length(f)+skipFirstNDC)));

figure;
subplot(2, 1, 1);
plot(time, MScan40Motion1);
title("Motion at pixel " + speakerIndices(1) + " of MScan40");
xlabel("\it{t} ({\its})");
ylabel("Displacement ({\itnm})");
xlim([0 time(end)]);

subplot(2, 1, 2);
plot(time, MScan40Motion2);
title("Motion at pixel " + speakerIndices(2) + " of MScan40");
xlabel("\it{t} ({\its})");
ylabel("Displacement ({\itnm})");
xlim([0 time(end)]);

freqs1 = islocalmax(MScan40Freq1Mag, "MaxNumExtrema", 40);
freqs2 = islocalmax(MScan40Freq2Mag, "MaxNumExtrema", 40);

disp("Frequencies at pixel " + speakerIndices(1) + " of MScan40 (in Hz):")
disp(sort(f(freqs1)));

disp("Frequencies at pixel " + speakerIndices(2) + " of MScan40 (in Hz):")
disp(sort(f(freqs2)));

figure;
subplot(2, 1, 1);
plot(f, MScan40Freq1Mag, f(freqs1), MScan40Freq1Mag(freqs1), "r*");
title("Frequency at pixel " + speakerIndices(1) + " of MScan40");
xlabel("\it{f} ({\itHz})");
ylabel("|D(\it{f})| ({\itdB})");
xlim([f(skipFirstNDC) f(end)]);

subplot(2, 1, 2);
plot(f, MScan40Freq2Mag, f(freqs2), MScan40Freq2Mag(freqs2), "r*");
title("Frequency at pixel " + speakerIndices(2) + " of MScan40");
xlabel("\it{f} ({\itHz})");
ylabel("|D(\it{f})| ({\itdB})");
xlim([f(skipFirstNDC) f(end)]);

figure;
plot(f, unwrap(angle(MScan40Freq1(skipFirstNDC+1:length(f)+skipFirstNDC))));
title("Phase Response of Speaker at pixel " + speakerIndices(1));
xlabel("\it{f} ({\itHz})");
ylabel("\angleD(\it{f}) (radians)");
xlim([f(skipFirstNDC) f(end)]);