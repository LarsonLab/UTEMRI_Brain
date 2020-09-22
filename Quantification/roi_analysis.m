%% brain ultrashort-T2 analysis

addpath(genpath('/Users/nikhil/Documents/ute_t2'));

filepath = '/Users/nikhil/Documents/ute_t2';
scandate = {'2019-10-17' '2019-10-31' '2019-11-21' '2020-01-22' '2020-01-28' '2020-02-03', ...
    '2020-06-02', '2020-06-10', '2020-06-17', '2020-06-24', '2020-06-26'}; % list datasets
for x = 1: length(scandate)
    scandate{x} = [scandate{x} '_images'];
end
prefix = 'DICOMs_uT2_';
parameter = {'fraction' 'T2' 'df' 'T1_fat'}; % check if folders have 'T1-fat' in name...

brain_ROI_list = {'right_ct' 'left_ct' 'gcc' 'bcc' 'scc' 'forceps_minor' 'forceps_major' 'right_ptr' 'left_ptr' 'right_pic' 'left_pic' ...
    'right_aic' 'left_aic' 'right_ec' 'left_ec' 'right_pu' 'left_pu' 'right_cd' 'left_cd' 'right_thalamus' 'left_thalamus'}; % brain ROIs

%% read and prepare files

% read in raw ROIs into proper data structures; calculate mean/std for each
% dataset
for i = 1:length(scandate)
    for j = 1:length(parameter)
        for k = 1:length(brain_ROI_list)
            roi_volumes(j).(brain_ROI_list{k}) = niftiread([filepath filesep scandate{i} filesep prefix parameter{j} filesep brain_ROI_list{k} '.nii']);
            roi_mean(j).(brain_ROI_list{k}) = mean(roi_volumes(j).(brain_ROI_list{k})(:)); % mean intensity of volume  
            roi_std(j).(brain_ROI_list{k}) = std(double(roi_volumes(j).(brain_ROI_list{k})(:))); % standard deviation of given ROI intensity 
            switch j
                % scale each parameter for proper quantification
                case 1
                    roi_mean(1).(brain_ROI_list{k}) = roi_mean(1).(brain_ROI_list{k}) ./ 1000; % fractional component scaling
                    roi_std(1).(brain_ROI_list{k}) = roi_std(1).(brain_ROI_list{k}) ./ 1000; % fractional component scaling
                case 2
                    roi_mean(2).(brain_ROI_list{k}) = roi_mean(2).(brain_ROI_list{k}) ./ 10; % ultrashort T2* scaling (ms)
                    roi_std(2).(brain_ROI_list{k}) = roi_std(2).(brain_ROI_list{k}) ./ 10; % ultrashort T2* scaling (ms)
                case 3
                    roi_mean(3).(brain_ROI_list{k}) = roi_mean(3).(brain_ROI_list{k}) ./ 10; % frequency shift scaling (Hz)
                    roi_std(3).(brain_ROI_list{k}) = roi_std(3).(brain_ROI_list{k}) ./ 10; % frequency shift scaling (hz)
                case 4
                    roi_mean(4).(brain_ROI_list{k}) = roi_mean(4).(brain_ROI_list{k}) ./ 10; % T1 scaling (ms)
                    roi_std(4).(brain_ROI_list{k}) = roi_std(4).(brain_ROI_list{k}) ./ 10; % T1 scaling (ms)
            end
        end
    end
    
    all_scans_roi_volumes{i} = roi_volumes;
    all_scans_roi_mean{i} = roi_mean;
    all_scans_roi_std{i} = roi_std;
    
    
end

% generate "scandate x brian_ROI x parameter" data matrix
for p = 1:length(parameter)
    for q = 1:length(brain_ROI_list)
        for r = 1:length(scandate)
        mean_matrix(r,q,p) = all_scans_roi_mean{r}(p).(brain_ROI_list{q}); 
            if p == 1
                mean_matrix(mean_matrix < 0 ) = nan;
            end
        end
    end
end



% combine right and left ROIs
for z = 1:length(parameter)
     for x = 1:length(brain_ROI_list)-1
        if x == 1
            combined_matrix(:,x,z) = (mean_matrix(:,x,z) + mean_matrix(:,x+1,z))/2;
        elseif (x < 3 || x > 7) && (mod(x,2) == 0) && x ~= 2
            combined_matrix(:,x,z) = (mean_matrix(:,x,z) + mean_matrix(:,x+1,z))/2;
        elseif (3 <= x) && (x <= 7)
            combined_matrix(:,x,z) = mean_matrix(:,x,z);
        end
    end
end
combined_matrix = reshape(nonzeros(combined_matrix), length(scandate), 13, length(parameter));

%% fraction

combined_brain_ROI_list = {'Corticospinal Tract' 'Genu Corpus Callosum' 'Body Corpus Callosum' 'Splenium Corpus Callosum' ...
    'Forceps Minor' 'Forceps Major' 'Post Thalamic Radiation' 'External Capsule' 'Anterior Internal capsule'...
    'Posterior Internal capsule' 'Putamen' 'Caudate' 'Thalamus'};

figure
% uT2 fraction
boxplot(x(:,1:10,1))
set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
title('Ultrashort-T2* Fractional Component across 11 healthy volunteers', 'fontsize', 25)
xlabel('Region of Interest', 'fontsize', 22)
xtickangle(30)
ylim([.02 .15])
ylabel('Fractional component value', 'fontsize', 22)
%% plotting

combined_brain_ROI_list = {'1. corticospinal tract' '2. genu corpus callosum' '3. body corpus callosum' '4. splenium corpus callosum' ...
    '5. forceps minor' '6. forceps major' '7. post thalamic radiation' '8. posterior internal capsule' '9. anterior internal capsule'...
    '10. external capsule' 'putamen' 'caudate' 'thalamus'};

figure
% uT2 fraction
subplot(221)
boxplot(combined_matrix(:,:,1))
set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
title('fractional component', 'fontsize', 28)
xlabel('Region of Interest', 'fontsize', 22)
xtickangle(30)
ylabel('Fractional component value', 'fontsize', 22)


% uT2 T2
subplot(222)
boxplot(combined_matrix(:,:,2))
set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
title('T2* (ms)', 'fontsize', 28)
xlabel('Region of Interest', 'fontsize', 22)
xtickangle(30)
ylabel('relaxation time', 'fontsize', 22)


% uT2 df
subplot(223)
boxplot(combined_matrix(:,:,3))
set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
title('frequency shift (Hz)', 'fontsize', 28)
xlabel('Region of Interest', 'fontsize', 22)
xtickangle(30)
ylabel('frequency shift', 'fontsize', 22)

% uT2 T1
subplot(224)
boxplot(combined_matrix(:,:,4))
set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
title('T1 (s)', 'fontsize', 28)
xlabel('Region of Interest', 'fontsize', 22)
xtickangle(30)
ylabel('relaxation time (s)', 'fontsize', 22)

%% 

% uT2 T1 fraction

% subplot(235)
% boxplot(combined_matrix(:,:,5))
% set(gca,'xticklabel',combined_brain_ROI_list, 'fontsize', 18);
% title('ultrashort component T1 (ms)', 'fontsize', 28)
% xlabel('Region of Interest', 'fontsize', 22)
% xtickangle(30)
% ylabel('relaxation time (ms)', 'fontsize', 22)

% plot right and left structures separately

% % fractional component
% figure(1)
% boxplot(mean_matrix(:,:,1))
% set(gca,'xticklabel',brain_ROI_list, 'fontsize', 18);
% title('ultrashort T2 fractional component', 'fontsize', 28)
% xlabel('Region of Interest', 'fontsize', 22)
% xtickangle(30)
% ylabel('Fractional component value', 'fontsize', 22)
% 
% 
% % T2*
% figure(2)
% boxplot(mean_matrix(:,:,2))
% set(gca,'xticklabel',brain_ROI_list, 'fontsize', 18);
% title('ultrashort T2 T2*', 'fontsize', 28)
% xlabel('Region of Interest', 'fontsize', 22)
% xtickangle(30)
% ylabel('Fractional component value', 'fontsize', 22)
% 
% % frequency shift
% figure(3)
% boxplot(mean_matrix(:,:,3))
% set(gca,'xticklabel',brain_ROI_list, 'fontsize', 18);
% title('ultrashort T2 frequency shift', 'fontsize', 28)
% xlabel('Region of Interest', 'fontsize', 22)
% xtickangle(30)
% ylabel('Fractional component value', 'fontsize', 22)