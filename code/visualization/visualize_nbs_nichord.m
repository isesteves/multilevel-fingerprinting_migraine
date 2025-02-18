nbs_analysis_path = '/home/iesteves/FC/data/results/nbs/analysis';
nichordfiles_path = '/home/iesteves/FC/files/nichord';
atlas_path = '/home/iesteves/FC/files/atlases';
nr_areas = 130;

c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;


nComp = 68;

edge_files_nichord = ['/home/iesteves/FC/files/nichord/edges-', num2str(nComp),'PCs'];
edge_files = ['/home/iesteves/FC/files/brainnetviewer/edges-', num2str(nComp),'PCs'];

%%
t_labels_H = readtable([atlas_path, '/network_labels_hemisphere.txt'], 'HeaderLines', 0, 'ReadVariableNames', false);
mapping_H =t_labels_H{:,1};
labels_H = unique(mapping_H);

c_labels_H =  readtable([atlas_path, '/coord_SchaeferSubCRB7100_130.txt'], 'HeaderLines', 0, 'ReadVariableNames', false);
coords = c_labels_H{:,1:3};

%% create files 

%comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'preic-ic'};
contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};
thresholds = [3.1, 4];
method = 'Extent';
components = nComp;

for t = 1:length(thresholds)
    threshold = thresholds(t);
    for comp = 1:length(components)
        nComp = components(comp);

        load([nbs_analysis_path, '/analysis-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_sig-metrics_networks_method-',method, '_cluster-', num2str(round(threshold))])
        countH_sig_FC_aux = sigmetrics.countH_sig_FC;
        percH_sig_FC_aux = sigmetrics.percH_sig_FC;
        for c = 1:length(comparisons)
            comparison = comparisons{c};
            contrast_type = contrast_types{c};

            countH_sig_FC = countH_sig_FC_aux(:,:,c) + tril(countH_sig_FC_aux(:,:,c),-1)';

            managefolders([nichordfiles_path, '/edges-', num2str(nComp), 'PCs/'], 'create')
            dlmwrite([nichordfiles_path, '/edges-', num2str(nComp), 'PCs', ...
                '/countH_edges_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '.txt'], countH_sig_FC, 'delimiter', ' ');

            percH_sig_FC = percH_sig_FC_aux(:,:,c) + tril(percH_sig_FC_aux(:,:,c),-1)';
            managefolders([nichordfiles_path, '/edges-', num2str(nComp), 'PCs/'], 'create')
            dlmwrite([nichordfiles_path, '/edges-', num2str(nComp), 'PCs', ...
                '/percH_edges_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '.txt'], percH_sig_FC, 'delimiter', ' ');

        end % comparisons
    end % components
end % thresholds

%% coord
mapping_ind = zeros(length(mapping_H),1);
coords_H = zeros(length(labels_H), 3);
for i = 1:numel(labels_H)
    matching_ind = ismember(mapping_H, labels_H{i});
    mapping_ind(matching_ind,1) = i;
    coords_H(i,:) = mean(coords(matching_ind,:));
end
dlmwrite([nichordfiles_path,...
    '/coords_H.txt'], coords_H, 'delimiter', ',');

labels_H_t = cell2table(labels_H);
writetable(labels_H_t, [nichordfiles_path,...
    '/labels_H.txt'], 'delimiter', ' ', 'WriteVariableNames', false);

%% Joining significant edges - individual edges

% ict-pm + post-pm
ictpm = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ic-prem_4.txt']);
postpm = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postic-prem_4.txt']);

ictpm_m = ictpm{:,:};
postpm_m = postpm{:,:};
ictpm_postpm = ictpm_m + 2*postpm_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictpm-postpm_4.txt'], ictpm_postpm,'Delimiter', ' ');

cmap = [c_ictal; c_postictal; 0 0 0];
figure;
imAlpha=ones(size(ictpm_postpm));
imAlpha(ismember(ictpm_postpm, 0))=0;
imagesc(ictpm_postpm, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap)
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}post \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])

% ict-pm + post-pm + postov-pm
postovpm = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_prem-mid_4.txt']);
postovpm_m = postovpm{:,:};

ictpm_postpm_postovpm = postovpm_m + 2*ictpm_m + 4*postpm_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictpm-postpm-postovpm_4.txt'], ictpm_postpm_postovpm,'Delimiter', ' ');

cmap = [c_midcycle; c_ictal; 0.8 0.8 0.8; c_postictal; 0.7 0.7 0.7; 0.6 0.6 0.6; 0 0 0];
figure;
imAlpha=ones(size(ictpm_postpm_postovpm));
imAlpha(ismember(ictpm_postpm_postovpm, 0))=0;
imagesc(ictpm_postpm_postovpm, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap)
caxis([1 7])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}post \color{black}+ HC-\color[rgb]{',num2str(c_midcycle),'}postov \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])


% pre-inter + ict-inter + post-inter
preinter = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preic-interic_4.txt']);
ictinter = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ic-interic_4.txt']);
postinter = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postic-interic_4.txt']);

preinter_m = preinter{:,:};
ictinter_m = ictinter{:,:};
postinter_m = postinter{:,:};

preinter_ictinter_postinter = preinter_m + 2*ictinter_m + 4*postinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-ictinter-postinter_4.txt'], preinter_ictinter_postinter,'Delimiter', ' ');


cmap = [c_preictal; c_ictal; 0.8 0.8 0.8; c_postictal; 0.7 0.7 0.7; 0.6 0.6 0.6; 0 0 0];
figure;
imAlpha=ones(size(preinter_ictinter_postinter));
imAlpha(ismember(preinter_ictinter_postinter, 0))=0;
imagesc(preinter_ictinter_postinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap)
caxis([1 7])
colorbar
title('M-pre + M-ict + M-post (vs M-inter)')
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{',num2str(c_postictal),'}post \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

% pre-inter + pre-ict
preict = readtable([edge_files, '/edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preic-ic_4.txt']);
preict_m = preict{:,:};

preinter_preict = preinter_m + 2*preict_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-preict_4.txt'], preinter_preict,'Delimiter', ' ');


cmap = [c_interictal; c_ictal; 0 0 0];
figure;
imAlpha=ones(size(preinter_preict));
imAlpha(ismember(preinter_preict, 0))=0;
imagesc(preinter_preict, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap)
caxis([1 3])
colorbar
title('M-inter + M-ict (vs M-pre)')
title(['M-\color[rgb]{',num2str(c_interictal),'}inter \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_preictal),'}pre\color{black})'])

% [row, col] = find(preinter_preict==3);
% common = unique(sort([row, col], 2), 'rows');
% 
% a = categorical(common(:));
% tabulate(a)

% pre-inter + ict-inter
preinter_ictinter = preinter_m + 2*ictinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-ictinter_4.txt'], preinter_ictinter,'Delimiter', ' ');

% pre-inter + post-inter
preinter_postinter = preinter_m + 2*postinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-postinter_4.txt'], preinter_postinter,'Delimiter', ' ');

% ict-inter + post-inter
ictinter_postinter = ictinter_m + 2*postinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictinter-postinter_4.txt'], ictinter_postinter,'Delimiter', ' ');

cmap = [c_preictal; c_ictal; 0 0 0];
figure('pos', [50 50 1700 400]);
ax = subplot(1,3,1);
imAlpha=ones(size(preinter_ictinter));
imAlpha(ismember(preinter_ictinter, 0))=0;
imagesc(preinter_ictinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax, cmap)
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_preictal; c_postictal; 0 0 0];
ax1 = subplot(1,3,2);
imAlpha=ones(size(preinter_postinter));
imAlpha(ismember(preinter_postinter, 0))=0;
imagesc(preinter_postinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax1, cmap)
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_ictal; c_postictal; 0 0 0];
ax2 = subplot(1,3,3);
imAlpha=ones(size(ictinter_postinter));
imAlpha(ismember(ictinter_postinter, 0))=0;
imagesc(ictinter_postinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax2, cmap)
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

%% Mig diff vs HC diff

% pm-postov vs pre-inter 
postovpm_preinter = postovpm_m + 2*preinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-preinter_4.txt'], postovpm_preinter,'Delimiter', ' ');

% pm-postov + ict-inter
postovpm_ictinter = postovpm_m + 2*ictinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-ictinter_4.txt'], postovpm_ictinter,'Delimiter', ' ');

% pm-postov + post-inter
postovpm_postinter = postovpm_m + 2*postinter_m;
dlmwrite([edge_files_nichord, '/joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-postinter_4.txt'], postovpm_postinter,'Delimiter', ' ');

cmap = [c_premenstrual; c_preictal; 0 0 0];
figure('pos', [50 50 1700 400]);
ax = subplot(1,3,1);
imAlpha=ones(size(postovpm_preinter));
imAlpha(ismember(postovpm_preinter, 0))=0;
imagesc(postovpm_preinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax, cmap)
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_preictal),'}pre \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_premenstrual; c_ictal; 0 0 0];
ax1 = subplot(1,3,2);
imAlpha=ones(size(postovpm_ictinter));
imAlpha(ismember(postovpm_ictinter, 0))=0;
imagesc(postovpm_ictinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax1, cmap)
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_premenstrual; c_postictal; 0 0 0];
ax2 = subplot(1,3,3);
imAlpha=ones(size(postovpm_postinter));
imAlpha(ismember(postovpm_postinter, 0))=0;
imagesc(postovpm_postinter, 'AlphaData', imAlpha)
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax2, cmap)
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

%% COUNT - Joining significant edges - individual edges

% ict-pm + post-pm
ictpm = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ic-prem_4.txt']);
postpm = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postic-prem_4.txt']);

ictpm_m = ictpm{:,:};
postpm_m = postpm{:,:};
ictpm_postpm = (double(ictpm_m>1) + 2*(postpm_m>1)); % do not change the order, or the colorscale in nichord will start at 2, since there is no overlap
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictpm-postpm_4.txt'], ictpm_postpm,'Delimiter', ' ');

cmap = [c_ictal; c_postictal; 0 0 0];
figure;
imAlpha=ones(size(ictpm_postpm));
imAlpha(ismember(ictpm_postpm, 0))=0;
imagesc(ictpm_postpm, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap)
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}post \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])

% ict-pm + post-pm + postov-pm
postovpm = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_prem-mid_4.txt']);
postovpm_m = postovpm{:,:};

ictpm_postpm_postovpm = reverse_matrix(double(postovpm_m>1) + 2*(ictpm_m>1) + 4*(postpm_m>1),7);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictpm-postpm-postovpm_4.txt'], ictpm_postpm_postovpm,'Delimiter', ' ');

cmap = [c_midcycle; c_ictal; 0.8 0.8 0.8; c_postictal; 0.7 0.7 0.7; 0.6 0.6 0.6; 0 0 0];
figure;
imAlpha=ones(size(ictpm_postpm_postovpm));
imAlpha(ismember(ictpm_postpm_postovpm, 0))=0;
imagesc(ictpm_postpm_postovpm, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 7])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}post \color{black}+ HC-\color[rgb]{',num2str(c_midcycle),'}postov \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])

% ict-pm + postovpm
ictpm_postovpm = reverse_matrix(double(postovpm_m>1) + 2*(ictpm_m>1),3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictpm-postovpm_4.txt'], ictpm_postovpm,'Delimiter', ' ');

cmap = [c_midcycle; c_ictal; 0 0 0];
figure;
imAlpha=ones(size(ictpm_postovpm));
imAlpha(ismember(ictpm_postovpm, 0))=0;
imagesc(ictpm_postovpm, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ HC-\color[rgb]{',num2str(c_midcycle),...
    '}postov \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])

% postpm + postovpm
postpm_postovpm = reverse_matrix(double(postovpm_m>1) + 2*(postpm_m>1),3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postpm-postovpm_4.txt'], postpm_postovpm,'Delimiter', ' ');

cmap = [c_midcycle; c_postictal; 0 0 0];
figure;
imAlpha=ones(size(postpm_postovpm));
imAlpha(ismember(postpm_postovpm, 0))=0;
imagesc(postpm_postovpm, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_postictal),'}post \color{black}+ HC-\color[rgb]{',num2str(c_midcycle),...
    '}postov \color{black}(vs HC-\color[rgb]{', num2str(c_premenstrual),'}pm\color{black})'])

% pre-inter + ict-inter + post-inter
preinter = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preic-interic_4.txt']);
ictinter = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ic-interic_4.txt']);
postinter = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postic-interic_4.txt']);

preinter_m = preinter{:,:};
ictinter_m = ictinter{:,:};
postinter_m = postinter{:,:};

preinter_ictinter_postinter = reverse_matrix(double(preinter_m>1) + 2*(ictinter_m>1) + 4*(postinter_m>1), 7);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-ictinter-postinter_4.txt'], preinter_ictinter_postinter,'Delimiter', ' ');


cmap = [c_preictal; c_ictal; 0.8 0.8 0.8; c_postictal; 0.7 0.7 0.7; 0.6 0.6 0.6; 0 0 0];
figure;
imAlpha=ones(size(preinter_ictinter_postinter));
imAlpha(ismember(preinter_ictinter_postinter, 0))=0;
imagesc(preinter_ictinter_postinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 7])
colorbar
title('M-pre + M-ict + M-post (vs M-inter)')
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{',num2str(c_postictal),'}post \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

% pre-inter + pre-ict
preict = readtable([edge_files_nichord, '/countH_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preic-ic_4.txt']);
preict_m = preict{:,:};

preinter_preict = reverse_matrix(double(preinter_m>1) + 2*(preict_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-preict_4.txt'], preinter_preict,'Delimiter', ' ');


cmap = [c_interictal; c_ictal; 0 0 0];
figure;
imAlpha=ones(size(preinter_preict));
imAlpha(ismember(preinter_preict, 0))=0;
imagesc(preinter_preict, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 3])
colorbar
title('M-inter + M-ict (vs M-pre)')
title(['M-\color[rgb]{',num2str(c_interictal),'}inter \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_preictal),'}pre\color{black})'])

% [row, col] = find(preinter_preict==3);
% common = unique(sort([row, col], 2), 'rows');
% 
% a = categorical(common(:));
% tabulate(a)

% pre-inter + ict-inter
preinter_ictinter = reverse_matrix(2*(preinter_m>1) + double(ictinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-ictinter_4.txt'], preinter_ictinter,'Delimiter', ' ');

% pre-inter + post-inter
preinter_postinter = reverse_matrix(2*(preinter_m>1) + double(postinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_preinter-postinter_4.txt'], preinter_postinter,'Delimiter', ' ');

% ict-inter + post-inter
ictinter_postinter = reverse_matrix(double(ictinter_m>1) + 2*(postinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_ictinter-postinter_4.txt'], ictinter_postinter,'Delimiter', ' ');

cmap = [c_ictal; c_preictal; 0 0 0];
figure('pos', [50 50 1700 400]);
ax = subplot(1,3,1);
imAlpha=ones(size(preinter_ictinter));
imAlpha(ismember(preinter_ictinter, 0))=0;
imagesc(preinter_ictinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_postictal; c_preictal; 0 0 0];
ax1 = subplot(1,3,2);
imAlpha=ones(size(preinter_postinter));
imAlpha(ismember(preinter_postinter, 0))=0;
imagesc(preinter_postinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax1, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_preictal),'}pre \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_ictal; c_postictal; 0 0 0];
ax2 = subplot(1,3,3);
imAlpha=ones(size(ictinter_postinter));
imAlpha(ismember(ictinter_postinter, 0))=0;
imagesc(ictinter_postinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax2, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['M-\color[rgb]{',num2str(c_ictal),'}ict \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

%% Mig diff vs HC diff

% pm-postov vs pre-inter 
postovpm_preinter = reverse_matrix(2*(postovpm_m>1) + double(preinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-preinter_4.txt'], postovpm_preinter,'Delimiter', ' ');

% pm-postov + ict-inter
postovpm_ictinter = reverse_matrix(2*(postovpm_m>1) + double(ictinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-ictinter_4.txt'], postovpm_ictinter,'Delimiter', ' ');

% pm-postov + post-inter
postovpm_postinter = reverse_matrix(2*(postovpm_m>1) + double(postinter_m>1), 3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_postovpm-postinter_4.txt'], postovpm_postinter,'Delimiter', ' ');

cmap = [c_preictal; c_premenstrual; 0 0 0];
figure('pos', [50 50 1700 400]);
ax = subplot(1,3,1);
imAlpha=ones(size(postovpm_preinter));
imAlpha(ismember(postovpm_preinter, 0))=0;
imagesc(postovpm_preinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_preictal),'}pre \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_ictal; c_premenstrual; 0 0 0];
ax1 = subplot(1,3,2);
imAlpha=ones(size(postovpm_ictinter));
imAlpha(ismember(postovpm_ictinter, 0))=0;
imagesc(postovpm_ictinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax1, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_ictal),'}ict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

cmap = [c_postictal; c_premenstrual; 0 0 0];
ax2 = subplot(1,3,3);
imAlpha=ones(size(postovpm_postinter));
imAlpha(ismember(postovpm_postinter, 0))=0;
imagesc(postovpm_postinter, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(ax2, cmap(end:-1:1,:))
caxis([1 3])
colorbar
title(['HC-\color[rgb]{',num2str(c_premenstrual),'}pm \color{black}(vs HC-\color[rgb]{', num2str(c_midcycle),...
    '}postov\color{black}) \color{black}+ M-\color[rgb]{', ...
    num2str(c_postictal),'}postict \color{black}(vs M-\color[rgb]{', num2str(c_interictal),'}inter\color{black})'])

%%
% (preinter + ictinter + postinter) + postov-pm
mig_hc = reverse_matrix(double((preinter_m.*ictinter_m.*postinter_m)>1) + 2*double(postovpm_m>1),3);
dlmwrite([edge_files_nichord, '/countH_joint_edges_SchaeferSubCRB7100_130_nbs_', num2str(nComp),'PCs_migall-hc_4.txt'], mig_hc,'Delimiter', ' ');

cmap = [1 0 0; 0 1 0; 0 0 0;];
figure;
imAlpha=ones(size(mig_hc));
imAlpha(ismember(mig_hc, 0))=0;
imagesc(mig_hc, 'AlphaData', imAlpha)
hold on;
%plot_atlas_labels('SchaeferSubCRB7100', 130)
colormap(cmap(end:-1:1,:))
caxis([1 3])
colorbar

