function fit_maps = generate_fit_maps(fit_filename, fit_plus_filename, AIC_thresh, im_biasfield)

plot_maps = 1;
write_flag = 0;

if nargin < 3 || isempty(AIC_thresh)
    AIC_thresh = Inf;
end

% visualization
load(fit_plus_filename);
fit_result2_plus = fit_result2;  %fit_result3_plus = fit_result3;
AIC2_plus = AIC2; %AIC3_plus = AIC3;
load(fit_filename);

if nargin < 4 || isempty(im_biasfield)
    im_biasfield = ones(imsize);
end


%%

Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};
Iax = imsize(3)/2; Icor = imsize(1)/2; Isag = imsize(2)/2;
%
% switch B0
%     case 3
%         sc_ut2 = [.055 .145]; sc_lt2 = [.2 1.4];
%                  sc_ut2_t2 = [0 1.4]; sc_lt2_t2 = [30 55];
%          sc_ut2_df = [-450 -200]; sc_lt2_df = [-100 100];
% Iax = 37; Icor = 38; Isag = 48;
% %    Icrop = {[7:103], [1:imsize(2)], [1:83]};     Icor = Icor -6;
%          Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};
%
%
%     case 7
%          sc_ut2 = [.06 .16]; sc_lt2 = [.2 1.4];
%          sc_ut2_t2 = [0 0.6]; sc_lt2_t2 = [10 40];
%          sc_ut2_df = [-1200 -700]; sc_lt2_df = [-250 250];
%          Iax = 36; Icor = 39; Isag = 49;
%          Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};
% end
%
%          sc_AIC = [200 320];

%%

%datamask = medfilt3(AIC2 < -230);
datamask = (AIC2_plus < AIC_thresh);
%datamask = imgaussfilt3(AIC2_plus, 2) < -250;


for Imaps = 1:7
    
    clear dataplot
    for Ix = 1:length(I)
        switch Imaps
            case 1
                dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ im_biasfield(I(Ix)); sc = [.05 .2]; root_fname = 'uT2_corrected';
            case 2
                dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ fit_result2_plus(Ix,1).rho; sc = [.05 .2]; root_fname = 'uT2_fraction';
            case 3
                dataplot(Ix) = fit_result2_plus(Ix,2).T2; sc = [0 1.0]; root_fname = 'uT2_T2';
            case 4
                dataplot(Ix) = -fit_result2_plus(Ix,2).df*1e3; sc = [-450 -200]; root_fname = 'uT2_df';
            case 5
                dataplot(Ix) = fit_result2(Ix,1).T2; sc = [30 50]; root_fname = 'lT2_T2';
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
    
    alphamask = zeros(imsize);
    alphamask(I) = datamask_temp;
    
    fit_maps.(root_fname) = zeros(imsize) - Inf;
    fit_maps.(root_fname)(I) = dataplot .* datamask_temp;
    
    implot = fit_maps.(root_fname)(Icrop{1}, Icrop{2}, Icrop{3});
    alphamask = alphamask(Icrop{1}, Icrop{2}, Icrop{3});
    
    % figure(1)
    % disp3d(implot,sc(1), sc(2), [Icor, Isag, Iax])
    
    if plot_maps
        %% plot maps
        
        
        h = figure;
        subplot(131),sb1= imagesc(flipud((implot(:,:,Iax))), sc); axis off equal
        set(sb1, 'AlphaData', squeeze(flipud((alphamask(:,:,Iax)))));
        subplot(132), sb2= imagesc(squeeze((implot(:,Isag,:))).', sc); axis off equal
        set(sb2, 'AlphaData', squeeze((alphamask(:,Isag,:))).');
        subplot(133), sb3= imagesc(squeeze((implot(Icor,:,:))).', sc); axis off equal
        set(sb3, 'AlphaData', squeeze(alphamask(Icor,:,:)).');
        colormap(hot)
        colorbar
        set(gcf, 'color', [1 1 1]*.8)
        % apply mask
        
        if write_flag
            print([foldername{Ifolder} '/' root_fname], '-dpdf')
        end
    end
end