function fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield)
% fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield)
%

% visualization
load(fit_plus_filename);
fit_result2_plus = fit_result2;  %fit_result3_plus = fit_result3;
AIC2_plus = AIC2; %AIC3_plus = AIC3;
load(fit_filename);

if nargin < 4 || isempty(im_biasfield)
    im_biasfield = ones(imsize);
end


%%

if nargin < 3 || isempty(AIC_thresh)
    datamask = 1;
    %    AIC_thresh = Inf;
else
    %datamask = medfilt3(AIC2 < -230);
    datamask = (AIC2_plus < AIC_thresh);
    %datamask = imgaussfilt3(AIC2_plus, 2) < -250;
end

for Imaps = 1:7
    
    clear dataplot
    for Ix = 1:length(I)
        switch Imaps
            case 1
                dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ im_biasfield(I(Ix)); sc = [.0 .04]; root_fname = 'uT2_corrected';
            case 2
                dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ fit_result2_plus(Ix,1).rho; sc = [.05 .25]; root_fname = 'uT2_fraction';
            case 3
                dataplot(Ix) = fit_result2_plus(Ix,2).T2; sc = [0 1.0]; root_fname = 'uT2_T2';
            case 4
                dataplot(Ix) = fit_result2_plus(Ix,2).df*1e3; sc = [-1100 -600]; root_fname = 'uT2_df';
            case 5
                dataplot(Ix) = fit_result2(Ix,1).T2; sc = [15 30]; root_fname = 'lT2_T2';
            case 6
                dataplot(Ix) = (fit_result1(Ix,1).df + fit_result2(Ix,1).df)*1e3; sc = [-100 100]; root_fname = 'fieldmap';% good g/w diffs
            case 7
                dataplot(Ix) = -AIC2_plus(Ix); sc = [150 350];  root_fname = 'AIC';
                
        end
        
        
    end
    
    switch Imaps
        case {6,7}
            datamask_temp = 1;
        otherwise
            datamask_temp = datamask;
            
    end
    
    fit_maps.(root_fname) = zeros(imsize) - Inf;
    fit_maps.(root_fname)(I) = dataplot .* datamask_temp;
    
    

    

end