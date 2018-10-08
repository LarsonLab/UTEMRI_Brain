function write_ute_brain_dicom(fit_maps, pfile_name, im_ute)
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
    
 
    ute_dicom(map, pfile_name, fieldName, 0, scaleFactor, seriesNumber);

   
end


im_ute = abs(im_ute); % calculate magnitude from complex number

% 
ute_dicom(im_ute, pfile_name, 'UTE', 0, 32767/max(abs(im_ute(:))), 0000);

end
