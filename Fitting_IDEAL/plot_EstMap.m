function [] = plot_EstMap(EstValue, brainmask_cut_test, idx)

EstMap.Rho1_map = flip(img_from_fit(EstValue.Rho(:,1), brainmask_cut_test));
EstMap.Rho2_map = flip(img_from_fit(EstValue.Rho(:,2), brainmask_cut_test));
EstMap.RhoF_map = flip(img_from_fit(EstValue.Rho_frac(:,1), brainmask_cut_test));
EstMap.T2s1_map = flip(img_from_fit(1000./EstValue.R2s(:,1), brainmask_cut_test));
EstMap.T2s2_map = flip(img_from_fit(1000./EstValue.R2s(:,2), brainmask_cut_test));
EstMap.T11_map = flip(img_from_fit(1./EstValue.R1(:,1), brainmask_cut_test));
EstMap.T12_map = flip(img_from_fit(1./EstValue.R1(:,2), brainmask_cut_test));
EstMap.Phi_map = flip(img_from_fit(EstValue.Phi(:,1), brainmask_cut_test));
EstMap.Rmse_map = flip(img_from_fit(EstValue.Rmse(:,1), brainmask_cut_test));
EstMap.Iter_map = flip(img_from_fit(EstValue.Iter(:,1), brainmask_cut_test));

figure('Position',[100 100 1600 600]);
subplot(241)
imagesc(EstMap.Rho1_map,[0 20]); axis off; axis tight; axis equal; title('Rho1'); colorbar;
subplot(242)
imagesc(EstMap.Rho2_map,[0 5]); axis off; axis tight; axis equal; title('Rho2'); colorbar;
subplot(243)
imagesc(EstMap.RhoF_map,[0 0.2]); axis off; axis tight; axis equal; title('Rho fraction'); colorbar;
subplot(244)
imagesc(EstMap.Phi_map,[-100 200]); axis off; axis tight; axis equal; title('delta B0'); colorbar;

subplot(247)
imagesc(EstMap.T2s1_map,[0 100]); axis off; axis tight; axis equal; title('T2s comp1'); colorbar;
subplot(248)
imagesc(EstMap.T2s2_map,[0 1]); axis off; axis tight; axis equal; title('T2s comp2'); colorbar;
subplot(245)
imagesc(EstMap.T11_map,[0 5]); axis off; axis tight; axis equal; title('T1 comp1'); colorbar;
subplot(246)
imagesc(EstMap.T12_map,[0 1]); axis off; axis tight; axis equal; title('T1 comp2'); colorbar;


figure('Position',[100 100 1200 600]);
subplot(141)
imagesc(EstMap.Rmse_map,[0 0.02]); axis off; axis tight; axis equal; title('RMSE'); colorbar;
subplot(142)
imagesc(EstMap.Iter_map); axis off; axis tight; axis equal; title('Iter num.'); colorbar;

if ~isempty(idx)
    subplot(143)
    imagesc(flip(img_from_fit(idx, brainmask_cut_test))); axis off; axis tight; axis equal;
    title('Clusters'); colorbar;
end
