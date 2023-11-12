function [AScan] = readscan(filePath, pixelsPerScan)
    fileID = fopen(filePath);
    
    AScan = fread(fileID, "uint16", "ieee-le");
    
    width = length(AScan) / pixelsPerScan;
    AScan = reshape(AScan, [pixelsPerScan, width]);
    
    fclose(fileID);
end