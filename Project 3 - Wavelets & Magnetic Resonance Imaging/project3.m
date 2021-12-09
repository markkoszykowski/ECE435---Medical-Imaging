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

%% I


CSI1T1Run1File = string(CSI1Files(contains(CSI1Files, "run-01") & contains(CSI1Files, "T1")));
CSI1T1Run2File = string(CSI1Files(contains(CSI1Files, "run-02") & contains(CSI1Files, "T1")));
CSI1T2File = string(CSI1Files(contains(CSI1Files, "T2")));

CSI1T1Run1 = niftiread(CSI1T1Run1File);
CSI1T1Run2 = niftiread(CSI1T1Run2File);
CSI1T2 = niftiread(CSI1T2File);

sagittalCrossSection = round(size(CSI1T1Run1, 1) / 2);
coronalCrossSection = round(size(CSI1T1Run1, 2) / 2);
axialCrossSection = round(size(CSI1T1Run1, 3) / 2);


figure;
subplot(2, 3, 1);
imshow(rot90(squeeze(CSI1T1Run1(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T1-weighted Run 1");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 2);
imshow(rot90(squeeze(CSI1T1Run1(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T1-weighted Run 1");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 3);
imshow(rot90(squeeze(CSI1T1Run1(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T1-weighted Run 1");
xlabel("{\itx}");
ylabel("{\ity}");

subplot(2, 3, 4);
imagesc(rot90(squeeze(CSI1T1Run1(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T1-weighted Run 1");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 5);
imagesc(rot90(squeeze(CSI1T1Run1(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T1-weighted Run 1");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 6);
imagesc(rot90(squeeze(CSI1T1Run1(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T1-weighted Run 1");
xlabel("{\itx}");
ylabel("{\ity}");


figure;
subplot(2, 3, 1);
imshow(rot90(squeeze(CSI1T1Run2(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T1-weighted Run 2");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 2);
imshow(rot90(squeeze(CSI1T1Run2(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T1-weighted Run 2");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 3);
imshow(rot90(squeeze(CSI1T1Run2(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T1-weighted Run 2");
xlabel("{\itx}");
ylabel("{\ity}");

subplot(2, 3, 4);
imagesc(rot90(squeeze(CSI1T1Run2(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T1-weighted Run 2");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 5);
imagesc(rot90(squeeze(CSI1T1Run2(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T1-weighted Run 2");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 6);
imagesc(rot90(squeeze(CSI1T1Run2(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T1-weighted Run 2");
xlabel("{\itx}");
ylabel("{\ity}");


figure;
subplot(2, 3, 1);
imshow(rot90(squeeze(CSI1T2(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T2-weighted");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 2);
imshow(rot90(squeeze(CSI1T2(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T2-weighted");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 3);
imshow(rot90(squeeze(CSI1T2(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T2-weighted");
xlabel("{\itx}");
ylabel("{\ity}");

subplot(2, 3, 4);
imagesc(rot90(squeeze(CSI1T2(sagittalCrossSection, :, :))));
title("Sagittal Slice " + sagittalCrossSection + ", T2-weighted");
xlabel("{\ity}");
ylabel("{\itz}");

subplot(2, 3, 5);
imagesc(rot90(squeeze(CSI1T2(:, coronalCrossSection, :))));
title("Coronal Slice " + coronalCrossSection + ", T2-weighted");
xlabel("{\itx}");
ylabel("{\itz}");

subplot(2, 3, 6);
imagesc(rot90(squeeze(CSI1T2(:, :, axialCrossSection))));
title("Axial Slice " + axialCrossSection + ", T2-weighted");
xlabel("{\itx}");
ylabel("{\ity}");

%% III


