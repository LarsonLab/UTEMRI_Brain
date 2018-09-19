function write_dicom(fit_maps, pfile_name)
% generate DICOM images from fit_maps


fields = fieldnames(fit_maps);

for i = 1:numel(fields)
    
    fieldName = fields{i};

    switch i
        case 1 %uT2_corrected
            scaleFactor = 1000;%intensiteis inconsistent
        case 2 %uT2_fraction
            scaleFactor = 1000;%intensiteis inconsistent
        case 3 %uT2_T2
            scaleFactor = 1000;% bg gray
        case 4 %uT2_df
            scaleFactor = 1000;% bg gray
        case 5 %lT2_T2
            scaleFactor = 1000;% bg gray
        case 6 %fieldmap
            scaleFactor = 100;%!!!
        case 7 %AIC
            scaleFactor = 100;%!!!
    end
    
    ute_dicom(fit_maps.(fieldName), pfile_name, fieldName, 0, scaleFactor);
    

end


end
