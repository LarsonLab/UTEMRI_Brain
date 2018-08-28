write_flag = 0;
AIC_thresh = 0;

%% visualization 
load([foldername{Ifolder} '/fit_results_plus.mat'])
fit_result2_plus = fit_result2;  fit_result3_plus = fit_result3; 
AIC2_plus = AIC2; AIC3_plus = AIC3;
load([foldername{Ifolder} '/fit_results.mat'])



%% estimate bias field, normalized & remove
imute = imall(:,:,:,1); 
imute_plus = imallplus(:,:,:,1); 
brainmask_open = imdilate(brainmask, strel('disk', 16));
brainmask_open = permute(imdilate(permute(brainmask_open, [3 2 1]), strel('disk', 6)), [3 2 1]);
imbias = imgaussfilt3(abs(imute).*brainmask_open,6);
Snormall = max(abs(imute(I)));

imbias = imbias / Snormall;

% figure(1)
% disp3d(imbias)

%%

switch B0
    case 3
        sc_ut2 = [.055 .145]; sc_lt2 = [.2 1.4];
                 sc_ut2_t2 = [0 1.4]; sc_lt2_t2 = [30 55];
         sc_ut2_df = [-450 -200]; sc_lt2_df = [-100 100]; 
Iax = 37; Icor = 38; Isag = 48;
%    Icrop = {[7:103], [1:imsize(2)], [1:83]};     Icor = Icor -6;
         Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};

    
    case 7
         sc_ut2 = [.06 .16]; sc_lt2 = [.2 1.4]; 
         sc_ut2_t2 = [0 0.6]; sc_lt2_t2 = [10 40];
         sc_ut2_df = [-1200 -700]; sc_lt2_df = [-250 250]; 
         Iax = 36; Icor = 39; Isag = 49;
         Icrop = {[1:imsize(1)], [1:imsize(2)], [1:imsize(3)]};
end

         sc_AIC = [200 320];

%%
for Imaps = 1:7

clear dataplot
for Ix = 1:NI
    switch Imaps
        case 1
            dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ imbias(I(Ix)); sc = sc_ut2; root_fname = 'uT2_corrected';
        case 2
            dataplot(Ix) = fit_result2_plus(Ix,2).rho ./ fit_result2_plus(Ix,1).rho; sc = sc_ut2; root_fname = 'uT2_fraction'; 
        case 3
     dataplot(Ix) = fit_result2_plus(Ix,2).T2; sc = sc_ut2_t2; root_fname = 'uT2_T2';
        case 4
     dataplot(Ix) = -fit_result2_plus(Ix,2).df*1e3; sc = sc_ut2_df; root_fname = 'uT2_df';
        case 5
     dataplot(Ix) = fit_result2(Ix,1).T2; sc = sc_lt2_t2; root_fname = 'lT2_T2';
        case 6
            dataplot(Ix) = (fit_result1(Ix,1).df + fit_result2(Ix,1).df)*1e3; sc = sc_lt2_df; root_fname = 'fieldmap';% good g/w diffs
        case 7
               dataplot(Ix) = -AIC2_plus(Ix); sc = sc_AIC;  root_fname = 'AIC';

    end
                

end

%datamask = medfilt3(AIC2 < -230);
datamask = (AIC2_plus < AIC_thresh(Ifolder));
%datamask = imgaussfilt3(AIC2_plus, 2) < -250;

switch Imaps
    case {6,7}
datamask = 1;
end

alphamask = zeros(imsize);
alphamask(I) = datamask;

implot = zeros(imsize) - Inf;
implot(I) = dataplot .* datamask;

implot = implot(Icrop{1}, Icrop{2}, Icrop{3});
alphamask = alphamask(Icrop{1}, Icrop{2}, Icrop{3});

% figure(1)
% disp3d(implot,sc(1), sc(2), [Icor, Isag, Iax])

%% generate maps


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