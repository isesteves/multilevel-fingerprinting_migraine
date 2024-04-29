%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 4 - DI COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% define paths
main_path = '/home/iesteves/FC';

template_path = [main_path, '/files/template_matrices'];
sFC_path = [main_path, '/data/results/sFC_Pearson/MIGcomplete-HCcomplete_SchaeferSubCRB7100_nonzero-50']; 
DI_path = [main_path, '/data/results/DI'];
pipelinelog_path = [main_path, '/logs/pipeline_logs'];

% create relevant folders if they don't already exist
managefolders([DI_path, '/DI_values'], 'create');
managefolders([DI_path, '/DI_corr'], 'create');
managefolders([DI_path, '/DI_FC_PCA'], 'create');

% variables
pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';
networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'Subcortical', 'Cerebellum', 'WB'};

% number of networks
nr_networks = length(networks);

% load Pearson correlation variables (iFC_all - areas x areas x subject_session, ordered by session type)
FC_file = load([sFC_path, '/Pearson-results-', atlas,'-', pipeline,'-bysubject.mat']);
nr_areas = FC_file.N_areas;
nr_subjects = FC_file.N_subjects;
filenames = FC_file.filenames_bysubject;
FC_all = FC_file.FC_all_bysubject;

% log file name
DIcomputation_log_filename = [pipelinelog_path,'/step4_DI_computation-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%% FC vectorization
% get the vectorized lower triangular matrix
FCvec =  conversion_vecmat(FC_all, nr_areas, 'mat2vec');

save([DI_path, '/DI_FC_PCA/FCvec.mat'], 'FCvec');

% demean data
mean_data = mean(FCvec);
group_data = FCvec-mean_data;

%% Networks mask
matrix_all = zeros(nr_areas, nr_areas, 9);
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 130 9];
for g = 1:size(areas,1)
    matrix_all(areas(g,1):areas(g,2),areas(g,1):areas(g,2),g) = g;
    matrix_all(:,:,g) = matrix_all(:,:,g) - diag(diag(matrix_all(:,:,g)));
end

matrix_FC_all = sum(matrix_all,3);

% figure; 
% imagesc(matrix_FC_all);
% hold on;
% plot_atlas_labels('SchaeferSubCRB7100', 130)

% get the vectorized lower triangular matrix of the networks template
network_mask_vec = conversion_vecmat(matrix_FC_all, nr_areas, 'mat2vec')';

%% PCA
% number of PCs, from 1 to the number of data samples
components = 1:nr_subjects;
nr_components = length(components);

% pre-allocate for concatenated FC matrices after PCA reconstruction 
FCvec_PCA_recon = zeros(((nr_areas-1)*nr_areas)/2, nr_subjects, nr_components);

% pre-allocate for nr_areas x nr_areas FC matrices for each subject after PCA reconstruction 
FC_PCA_recon = zeros(nr_areas, nr_areas, nr_subjects, nr_components);

% pre-allocate for explained variance as a function of number of components
PCA_VE = zeros(nr_components, 1);

% pre-allocate for identifiability matrix from PCA reconstructed matrix
corr_PCA_recon_networks = zeros(size(group_data, 2), size(group_data, 2), nr_networks+1, nr_components);

% PCA
[coeffs, score, latent,~,explained,~] = pca(group_data, 'NumComponents', components(end));
for comp = 1:nr_components
      
    % multiply each PC by the corresponding weight
    PCA_matrix = score(:,1:components(comp)) * coeffs(:,1:components(comp))';
    
    % sum data original mean
    group_PCA_recon = bsxfun(@plus, PCA_matrix, mean_data);
    
    % store concatenated vectorized FC for each #PCs 
    FCvec_PCA_recon(:,:,comp) = group_PCA_recon;
    
    % store nr_areas x nr_areas FC matrix for each subject after PCA reconstruction 
    FC_PCA_recon(:,:,:,comp) = conversion_vecmat(group_PCA_recon, nr_areas, 'vec2mat');
    
    % explained variance
    exp_var = cumsum(explained(1:comp));
    
    % store explained variance
    PCA_VE(comp, 1) = exp_var(end);
    
    % for each network/WB
    for n = 1:nr_networks
        if n <= nr_networks-1
            network_data = group_PCA_recon(network_mask_vec==n,:); 
        else
            network_data = group_PCA_recon;
        end
        
        % compute Pearson correlation between vectorized FC matrices
        corr_val = corr(network_data, network_data, 'Type', 'Pearson');
        
        % store identifiability matrix for each netwrok/WB and #PCs
        corr_PCA_recon_networks(:,:,n,comp) = corr_val;

    end   
end

save([DI_path, '/DI_FC_PCA/PCA_coeffs.mat'], 'coeffs');
save([DI_path, '/DI_FC_PCA/PCA_score.mat'], 'score');
save([DI_path, '/DI_FC_PCA/FCvec_PCA_recon.mat'], 'FCvec_PCA_recon');
save([DI_path, '/DI_FC_PCA/FC_PCA_recon.mat'], 'FC_PCA_recon');
save([DI_path, '/DI_FC_PCA/corr_PCA_recon_networks.mat'], 'corr_PCA_recon_networks');
save([DI_path, '/DI_FC_PCA/PCA_VE.mat'], 'PCA_VE');

% Log file
DIcomputation_log_msg = ['------ ', 'Finished and stored PCA results for: ', pipeline, '-', atlas, ', using ', num2str(nr_networks), ' networks and ', num2str(nr_components), ' PCs'];
disp(DIcomputation_log_msg);
logCustom(DIcomputation_log_msg, DIcomputation_log_filename)

%% Load template matrices
% general (for Idiff)
load([template_path, '/within_group'])
load([template_path, '/within_subject'])
load([template_path, '/within_session'])
load([template_path, '/within_menstrual_session'])
load([template_path, '/between_group_within_session'])

% subdivisions to compute average correlation (to explore more specific aspects, for example, for specific sessions)
load([template_path, '/wgroup_corr'])
load([template_path, '/wgroup_bsession_corr'])
load([template_path, '/wsubject_bsession_corr'])
load([template_path, '/wsession_corr'])
load([template_path, '/bgroup_corr'])

%% Idiff
% For each level, computation of Idiff, Iself/Iwithin and Iothers for each network and
% number of PCs used for reconstruction + storage

% pre-allocate space to store DI values
within_session_DI = zeros(3, nr_networks+1, nr_components);
within_subject_DI = zeros(3, nr_networks+1, nr_components);
within_group_DI = zeros(3, nr_networks+1, nr_components);
within_menstrual_session_DI = zeros(3, nr_networks+1, nr_components);
between_group_within_session_DI = zeros(3, nr_networks+1, nr_components);

% for each network + WB
for n = 1:nr_networks
    % for each number of PCs used for reconstruction
    for comp = 1:nr_components
        corr_matrix = corr_PCA_recon_networks(:,:,n,comp); % get the corresponding multilevel identifiability matrix

        [within_session_Idiff, within_session_Iself, within_session_Iothers] = diff_identifiability(corr_matrix, within_session); % within session
        [within_subject_Idiff, within_subject_Iself, within_subject_Iothers] = diff_identifiability(corr_matrix, within_subject); % within subject
        [within_group_Idiff, within_group_Iself, within_group_Iothers] = diff_identifiability(corr_matrix, within_group); % within group
        [within_menstrual_session_Idiff, within_menstrual_session_Iself, within_menstrual_session_Iothers] = diff_identifiability(corr_matrix, within_menstrual_session); % within menstrual session
        [between_group_within_session_Idiff, between_group_within_session_Iself, between_group_within_session_Iothers] = diff_identifiability(corr_matrix, between_group_within_session); % between groups (within menstrual session)
        
        within_session_DI(:,n,comp) = [within_session_Idiff, within_session_Iself, within_session_Iothers]'; % within session
        within_subject_DI(:,n,comp) = [within_subject_Idiff, within_subject_Iself, within_subject_Iothers]; % within subject
        within_group_DI(:,n,comp) = [within_group_Idiff, within_group_Iself, within_group_Iothers]; % within group
        within_menstrual_session_DI(:,n,comp) = [within_menstrual_session_Idiff, within_menstrual_session_Iself, within_menstrual_session_Iothers]; % within menstrual session
        between_group_within_session_DI(:,n,comp) = [between_group_within_session_Idiff, between_group_within_session_Iself, between_group_within_session_Iothers]; % between groups (within menstrual session)

    end % components
end % networks

% store DI values
save([DI_path, '/DI_values/within_session_DI.mat'], 'within_session_DI');
save([DI_path, '/DI_values/within_subject_DI.mat'], 'within_subject_DI');
save([DI_path, '/DI_values/within_group_DI.mat'], 'within_group_DI');
save([DI_path, '/DI_values/within_menstrual_session_DI.mat'], 'within_menstrual_session_DI');
save([DI_path, '/DI_values/between_group_within_session_DI.mat'], 'between_group_within_session_DI');

% Log file
DIcomputation_log_msg = ['------ ', 'Finished and stored DI values for: ', pipeline, '-', atlas, ', using ', num2str(nr_networks), ' networks and ', num2str(nr_components), ' PCs'];
disp(DIcomputation_log_msg);
logCustom(DIcomputation_log_msg, DIcomputation_log_filename)

%% Correlation
% For each level, computation of Idiff, Iself/Iwithin and Iothers for each network and
% number of PCs used for reconstruction + storage

% pre-allocate space to store correlation values
corr_val_wgroup = cell(nr_networks+1, nr_components);
corr_val_wgroup_bsession = cell(nr_networks+1, nr_components);
corr_val_wsubject_bsession = cell(nr_networks+1, nr_components);
corr_val_wsession = cell(nr_networks+1, nr_components);
corr_val_bgroup = cell(nr_networks+1, nr_components);

% for each network + WB
for n = 1:nr_networks
    % for each number of PCs used for reconstruction
    for comp = 1:nr_components
        corr_matrix = corr_PCA_recon_networks(:,:,n,comp); % get the corresponding multilevel identifiability matrix
        
        corr_val_wgroup{n,comp} = diff_identifiability_corr(corr_matrix, wgroup_corr); % within group (ex. patients)
        corr_val_wgroup_bsession{n,comp} = diff_identifiability_corr(corr_matrix, wgroup_bsession_corr); % within group between sessions (ex. preictal vs ictal)
        corr_val_wsubject_bsession{n,comp} = diff_identifiability_corr(corr_matrix, wsubject_bsession_corr); % within subject between sessions (ex. preictal vs ictal same subject)
        corr_val_wsession{n,comp} = diff_identifiability_corr(corr_matrix, wsession_corr); % within session (ex. preictal)
        corr_val_bgroup{n,comp} = diff_identifiability_corr(corr_matrix, bgroup_corr); % between groups (ex. preictal vs premenstrual)
    end
end

save([DI_path, '/DI_corr/corr_val_wgroup.mat'], 'corr_val_wgroup');
save([DI_path, '/DI_corr/corr_val_wgroup_bsession.mat'], 'corr_val_wgroup_bsession');
save([DI_path, '/DI_corr/corr_val_wsubject_bsession.mat'], 'corr_val_wsubject_bsession');
save([DI_path, '/DI_corr/corr_val_wsession.mat'], 'corr_val_wsession');
save([DI_path, '/DI_corr/corr_val_bgroup.mat'], 'corr_val_bgroup');

% Log file
DIcomputation_log_msg = ['------ ', 'Finished and stored DI correlation for: ', pipeline, '-', atlas, ', using ', num2str(nr_networks), ' networks and ', num2str(nr_components), ' PCs'];
disp(DIcomputation_log_msg);
logCustom(DIcomputation_log_msg, DIcomputation_log_filename)
