% 7T opera ROI measurements

addpath(genpath('/Users/nikhi/Documents/opera'));
path = '/Users/nikhil/Documents/opera';

sub = {'subj001', 'subj002', 'subj003', 'subj004', 'subj005', 'subj006', ...
    'subj007', 'subj008', 'subj009', 'subj010'};

% import ROI nifti volumes

% subject test
for n = 2:2
    nLesions = 5; % this will change depending on the subject
    for i = 1:nLesions
        lesion_volumes{i} = niftiread([path filesep sub{n} '/lesion_' int2str(i) '.nii']);
        % NAWM
        nawm_volume = double(niftiread([path filesep sub{n} '/nawm.nii']));
        
        lesion_volumes{i} = double(lesion_volumes{i}); % convert volumes to type double
        lesion_mean(i) = mean(lesion_volumes{i}(:));
        lesion_std(i) = std(lesion_volumes{i}(:));

    end

    % boxplot preparation
    matrix = [lesion_volumes{1}(:) / mean(nawm_volume(:)); lesion_volumes{2}(:) / mean(nawm_volume(:)); lesion_volumes{3}(:) / mean(nawm_volume(:)); lesion_volumes{4}(:) / mean(nawm_volume(:)); ...
        lesion_volumes{5}(:) / mean(nawm_volume(:)); nawm_volume(:) / mean(nawm_volume(:))];

    bp_prepare = [ones(size(lesion_volumes{1}(:))); 2*ones(size(lesion_volumes{2}(:))); 3*ones(size(lesion_volumes{3}(:))); ...
        4*ones(size(lesion_volumes{4}(:))); 5*ones(size(lesion_volumes{5}(:))); ...
        6*ones(size(nawm_volume(:)))];
    
    subject_cell{4,1} = matrix;
    subject_cell{4,2} = bp_prepare;

    figure
    boxplot(subject_cell{4,1},subject_cell{4,2})
    set(gca,'xticklabel',{'lesion 1', 'lesion 2', 'lesion 3', 'lesion 4', 'lesion 5', 'NAWM'}, 'fontsize', 18);
    title(['Subject ' int2str(n) ' Lesion intensities'], 'fontsize', 21)
    xlabel('Lesion number', 'fontsize', 18)
    xtickangle(30)
    ylabel('Normalized Intensity', 'fontsize', 18)

end

%%
figure
subplot(221)
boxplot(subject_cell{1,1},subject_cell{1,2})
set(gca,'xticklabel',{'lesion 1', 'lesion 2', 'lesion 3', 'lesion 4', 'lesion 5', 'lesion 6','NAWM'}, 'fontsize', 18);
title(['Subject 1 Lesion intensities'], 'fontsize', 21)
xlabel('Lesion number', 'fontsize', 18)
xtickangle(30)
ylabel('Normalized Intensity', 'fontsize', 18)

subplot(222)
boxplot(subject_cell{2,1},subject_cell{2,2})
set(gca,'xticklabel',{'lesion 1', 'lesion 2', 'lesion 3', 'lesion 4', 'lesion 5', 'lesion 6','NAWM'}, 'fontsize', 18);
title(['Subject 2 Lesion intensities'], 'fontsize', 21)
xlabel('Lesion number', 'fontsize', 18)
xtickangle(30)
ylabel('Normalized Intensity', 'fontsize', 18)

subplot(223)
boxplot(subject_cell{3,1},subject_cell{3,2})
set(gca,'xticklabel',{'lesion 1', 'lesion 2', 'lesion 3', 'lesion 4', 'lesion 5', 'lesion 6','NAWM'}, 'fontsize', 18);
title(['Subject 3 Lesion intensities'], 'fontsize', 21)
xlabel('Lesion number', 'fontsize', 18)
xtickangle(30)
ylabel('Normalized Intensity', 'fontsize', 18)

subplot(224)
boxplot(subject_cell{4,1},subject_cell{4,2})
set(gca,'xticklabel',{'lesion 1', 'lesion 2', 'lesion 3', 'lesion 4', 'lesion 5','NAWM'}, 'fontsize', 18);
title(['Subject 4 Lesion intensities'], 'fontsize', 21)
xlabel('Lesion number', 'fontsize', 18)
xtickangle(30)
ylabel('Normalized Intensity', 'fontsize', 18)


