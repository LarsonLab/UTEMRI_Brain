function write_dicom(fit_maps, pfile_name)
% generate DICOM images from fit_maps


fields = fieldnames(fit_maps);

for i = 1:numel(fields)
    
    fieldName = fields{i};

    switch i
        case 1 %uT2_corrected
            scaleFactor = 1000;
        case 2 %uT2_fraction
            scaleFactor = 1000;
        case 3 %uT2_T2
            scaleFactor = 100;
        case 4 %uT2_df
            scaleFactor = 10;
        case 5 %lT2_T2
            scaleFactor = 100; 
        case 6 %fieldmap
            scaleFactor = 10;
        case 7 %AIC
            scaleFactor = 10;
    end
    
    seriesNumber = i * 1000;
    
    map = fit_maps.(fieldName);
    
%     x = map(find(map));
%     figure;histogram(x);

    % find minimum pixel value, and set backgound = px_min
    % otherwise background in DICOM will look gray
    px_min = min(map(map~=-inf));
    px_max = max(map(map~=+inf));
    map(map==-inf) = px_min;
    
    ute_dicom(map, pfile_name, fieldName, 0, scaleFactor, seriesNumber);

    
    

end


end
