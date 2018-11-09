function map_scales = plot_fit_maps(fit_maps, im_ute, B0, Iplot, write_flag, write_root)
%plot_fit_maps(fit_maps, im_ute, [B0 or map_scales] ,write_flag)

%%
if nargin > 2 && isstruct(B0)
    map_scales = B0;
else
    map_scales = struct();
    if  nargin < 3 || isempty(B0)
        B0 = 3;
    end
        switch B0
    case 3
        map_scales.uT2_corrected = [0, .05];
        map_scales.uT2_fraction = [.1 .2];
        map_scales.uT2_T2 = [0 1.4];
        map_scales.uT2_df =[-450 -200];
        map_scales.lT2_T2 = [20 55];
        map_scales.fieldmap = [-100 100];
        map_scales.AIC = [100 300];
        
        
    case 7
        map_scales.uT2_corrected = [0, .04];
        map_scales.uT2_fraction = [.05 .25];
        map_scales.uT2_T2 = [0 0.6];
        map_scales.uT2_df =[-1200 -700];
        map_scales.lT2_T2 = [10 40];
        map_scales.fieldmap = [-250 250];
        map_scales.AIC = [100 300];
end
end

if nargin < 5 || isempty(write_flag)
    write_flag = 0;
end



root_fnames = {'uT2_corrected', 'uT2_fraction', 'uT2_T2', 'uT2_df', 'lT2_T2','fieldmap','AIC'};

imsize = size(fit_maps.(root_fnames{1}));
Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};
if nargin < 4 || isempty(Iplot)
    Iax = imsize(3)/2; Icor = imsize(1)/2; Isag = imsize(2)/2;
else
    Iax = Iplot(3); Icor = Iplot(1); Isag = Iplot(2);
end
    
for Imaps = 1:length(root_fnames)
    root_fname = root_fnames{Imaps};
    if isfield(fit_maps, root_fname)
    
    implot = fit_maps.(root_fname)(Icrop{1}, Icrop{2}, Icrop{3});
    
    alphamask = zeros(size(implot));
    alphamask(implot > -Inf) = 1;
    
    % figure(1)
    % disp3d(implot,sc(1), sc(2), [Icor, Isag, Iax])
    
    
    h = figure;
    subplot(131),sb1= imagesc(flipud((implot(:,:,Iax))), map_scales.(root_fname)); axis off equal
    set(sb1, 'AlphaData', squeeze(flipud((alphamask(:,:,Iax)))));
    subplot(132), sb2= imagesc(squeeze((implot(:,Isag,:))).', map_scales.(root_fname)); axis off equal
    set(sb2, 'AlphaData', squeeze((alphamask(:,Isag,:))).');
    subplot(133), sb3= imagesc(squeeze((implot(Icor,:,:))).', map_scales.(root_fname)); axis off equal
    set(sb3, 'AlphaData', squeeze(alphamask(Icor,:,:)).');
    colormap(hot)
    colorbar
    set(gcf, 'color', [1 1 1]*.8)
    % apply mask
    
    if write_flag
        print([write_root, root_fname], '-dpdf')
    end
    end
end

if ~isempty(im_ute)
    implot= im_ute(Icrop{1}, Icrop{2}, Icrop{3});
    h = figure;
    subplot(131),sb1= imagesc(flipud((implot(:,:,Iax)))); axis off equal
    subplot(132), sb2= imagesc(squeeze((implot(:,Isag,:))).'); axis off equal
    subplot(133), sb3= imagesc(squeeze((implot(Icor,:,:))).'); axis off equal
    colormap(gray)
    colorbar
    set(gcf, 'color', [1 1 1]*.8)
    % apply mask
    
    if write_flag
        print([write_root 'im_ute'], '-dpdf')
    end
end

end

