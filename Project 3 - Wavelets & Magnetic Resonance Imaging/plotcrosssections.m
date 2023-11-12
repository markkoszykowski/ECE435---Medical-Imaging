function plotcrosssections(image, xSlice, ySlice, zSlice, name)
%PLOTCROSSSECTIONS Plots sagittal, coronal, and axial slices of image stack
    figure;
    subplot(2, 3, 1);
    imshow(rot90(squeeze(image(xSlice, :, :))));
    title("Sagittal Slice " + xSlice + ", " + name);
    xlabel("{\ity}");
    ylabel("{\itz}");

    subplot(2, 3, 2);
    imshow(rot90(squeeze(image(:, ySlice, :))));
    title("Coronal Slice " + ySlice + ", " + name);
    xlabel("{\itx}");
    ylabel("{\itz}");

    subplot(2, 3, 3);
    imshow(rot90(squeeze(image(:, :, zSlice))));
    title("Axial Slice " + zSlice + ", " + name);
    xlabel("{\itx}");
    ylabel("{\ity}");

    subplot(2, 3, 4);
    imagesc(rot90(squeeze(image(xSlice, :, :))));
    title("Sagittal Slice " + xSlice + ", " + name);
    xlabel("{\ity}");
    ylabel("{\itz}");

    subplot(2, 3, 5);
    imagesc(rot90(squeeze(image(:, ySlice, :))));
    title("Coronal Slice " + ySlice + ", " + name);
    xlabel("{\itx}");
    ylabel("{\itz}");

    subplot(2, 3, 6);
    imagesc(rot90(squeeze(image(:, :, zSlice))));
    title("Axial Slice " + zSlice + ", " + name);
    xlabel("{\itx}");
    ylabel("{\ity}");
end