% add multiple comparisons correction toolbox
addpath('/home/iesteves/stats/MultipleTestingToolbox/MultipleTestingToolbox')

main_path = '/home/iesteves/FC';

DI_path = [main_path, '/data/results/DI'];
nbs_stats_path = [main_path, 'data/results/nbs/stats'];
fig_path = [main_path, 'figures/clinical'];

managefolders(fig_path, 'create');

% load patient/control clinical information
load('/home/iesteves/mign2treat_sample/MIGN2TREATcontrolsample.mat')
load('/home/iesteves/mign2treat_sample/MIGN2TREATpatientsample.mat')

% clinical variables
clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'yrs', 'yrs', '#/month', 'hrs', 'score'};

% patients
patients = MIGN2TREATpatientsample;
complete_patients_ind = strcmp(patients.status, 'complete');
complete_patients = patients(complete_patients_ind, :);
preictal = complete_patients(strcmp(complete_patients.session, 'ses-preictal'),:);
postictal = complete_patients(strcmp(complete_patients.session, 'ses-postictal'),:);
ictal = complete_patients(strcmp(complete_patients.session, 'ses-ictal'),:);
interictal = complete_patients(strcmp(complete_patients.session, 'ses-interictal'),:);
nr_patients = length(unique(complete_patients.subject));

% controls
controls = MIGN2TREATcontrolsample;
complete_controls_ind = strcmp(controls.status, 'complete');
complete_controls = controls(complete_controls_ind, :);
premenstrual = complete_controls(strcmp(complete_controls.session, 'ses-premenstrual'),:);
midcycle = complete_controls(strcmp(complete_controls.session, 'ses-midcycle'),:);

nr_controls = length(unique(complete_controls.subject))-1; % 025 was excluded

% load Idiff information
load([DI_path, '/DI_corr/corr_val_wsubject_bsession.mat']) % within subject, between sessions (corr_wsubject_bsession)
load([DI_path, '/DI_corr/corr_val_bgroup.mat']) % between groups (corr_bgroup)

% load FC data
%load([DI_path, '/DI_FC_PCA/FCvec_PCA_recon.mat']);
load([DI_path, '/DI_FC_PCA/FC_PCA_recon.mat']);

% FC variables
n = 10; % whole brain
comp = 19; % component

corr_type = 'Spearman';

% figure variables
font = 10;
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;


%% clinical features
% > Idiff
% - within subject, between sessions avg - 1 value per subject
% - Mig vs HC - 4 values per subject (each session)
% - Mig vs HC - 1 value per subject (avg)
%
% > FC
% - mask ict vs pm - 1 value per subject
% - mask post vs pm - 1 value per subject


%% Within-subject
nr_val_wpatient = 6; % 4 sessions, pairwise

% specify network and number of PCs
corr_wsubject = corr_val_wsubject_bsession{n, comp};

% Idiff values for patients and controls
m_val_patient = mean(corr_wsubject(1:nr_patients,1:nr_val_wpatient),2);
m_val_control = corr_wsubject(:,end);

c_patient = [1 0 0];
c_control = [0 0 1];

d = 1;
figure('pos', [50 50 1450 800])
for c = 1:length(clinical_names)

    clinical_name = clinical_names{c};

    clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
    subplot(3,5,d)

    scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
    lsline
    hold on;
    scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
    [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
    text(mean(clinical_data),mean(m_val_patient), ...
        ['r_s = ', num2str(round(r,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
    
    ht = title(clinical_name);
    ht.Interpreter = 'none';

    if c == 1
        ylabel({'Avg correlation', '(within-subject)'})
    end
    d = d +1;

    xlabel(clinical_units{c})

end

%print([fig_path, '/general-clinical_ses-', session_patient,'_network-global_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

%[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(4,:,:), 0.05, false);  


%% between groups (more similar to Iclinical)
sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

mean_patient_allsessions = zeros(10,4);

pvec = zeros(length(sessions), length(clinical_names));

d = 1;
figure('pos', [50 50 1450 800])

for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    corr_bgroup = corr_val_bgroup{n, comp};
    corr_session_vec = corr_bgroup(:,s);
    corr_session = reshape(corr_session_vec, nr_controls, nr_patients);

    m_val_patient = mean(corr_session)';
    mean_patient_allsessions(:,s) = m_val_patient;

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
        subplot(4,5,d)

        scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
        lsline
        hold on;
        %scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
        [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
        pvec(s,c) = p;
        text(mean(clinical_data),mean(m_val_patient), ...
            ['r_s = ', num2str(round(r,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)

        ht = title(clinical_name);
        ht.Interpreter = 'none';

        if c == 1
            ylabel({'Avg correlation', '(between-group)'})
        end
        d = d +1;

        xlabel(clinical_units{c})
        ylim([0.4 0.9])

    end
end

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(3,:), 0.05, false);

%% avg between
m_val_patient = mean(mean_patient_allsessions,2);
m_val_control = corr_wsubject(:,end);

c_patient = [1 0 0];
c_control = [0 0 1];
font = 10;

d = 1;
figure('pos', [50 50 1450 800])

for c = 1:length(clinical_names)

    clinical_name = clinical_names{c};

    clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
    subplot(3,5,d)

    scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
    lsline
    hold on;
    %scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
    [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
    text(mean(clinical_data),mean(m_val_patient), ...
        ['r_s = ', num2str(round(r,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
    
    ht = title(clinical_name);
    ht.Interpreter = 'none';

    if c == 1
        ylabel({'Avg correlation', '(between-group avg)'})
    end
    d = d +1;

    xlabel(clinical_units{c})

end

%% Avg sig edges

contrast_type = 1;
method = 'Extent';
threshold = 4;

sessions = {'ictal', 'postictal'};

% stats
hvec = zeros(length(sessions), length(clinical_names));
pvec = zeros(length(sessions), length(clinical_names));

figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Avg sig FC edges vs General Clinical - No Thresh - ', num2str(comp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;
for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'ic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'postic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    load([nbs_stats_path, '/stats-', num2str(comp), 'PCs/nbs_', num2str(comp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
    adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
    
%     figure('pos', [50, 50, 300 450]);
%     subplot(3,1,[1 2])
%     imAlpha=ones(size(adj_matrix));
%     imAlpha(adj_matrix==0)=0;
%     imagesc(tril(adj_matrix,0),'AlphaData', imAlpha)
%     hold on;
%     plot_atlas_labels(atlas, nr_areas);
%     cmap = gray(256);
%     colormap(flipud(cmap))
%     set(gca, 'XTickLabel', '')
%     set(gca, 'YTickLabel', '')
%     title(comparison)

    all_data = FC_PCA_recon(:,:,:,comp);
    %all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient = zeros(nr_patients, 1);
    area = all_data(:,:,ind_patient:4:40);
    for k = 1:nr_patients
        area_patient = area(:,:,k);
        area_sig = area_patient(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_patient(k) = val;      
    end
    
    m_val_control = zeros(nr_controls, 1);
    area = all_data(:,:,ind_control:2:68);
    for k = 1:nr_controls
        area_control = area(:,:,k);
        area_sig = area_control(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_control(k) = val;      
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};
       
        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
        subplot(2,5,d)

        scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
        lsline
        hold on;
        scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
        [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
         % stats
        pvec(s, c) = p;
        text(mean(clinical_data),mean(m_val_patient), ...
            ['r_s = ', num2str(round(r,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
        ylim([0.15 0.75])
        if s == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
        end
        if c == 1
            ylabel(['\bf \color[rgb]{', num2str(c_patient),'}', session_patient]);
        end

        if s == 2
            xlabel(clinical_units{c})
        end
        
        d = d +1;
    end
end

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(2, :), 0.05, true); 

                                  

