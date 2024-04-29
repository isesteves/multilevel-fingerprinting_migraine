%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 5 - NETWORK-BASED STATISTIC (NBS) STATS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

addpath(genpath('/home/iesteves/toolboxes/NBS1.2/NBS1.2'))

main_path = '/home/iesteves/FC';

DI_path = [main_path, '/data/results/DI'];
nbs_files_path = [main_path, '/files/nbs'];
nbs_stats_path = [main_path, '/data/results/nbs/stats'];
atlas_path = [main_path, '/files/atlases'];
pipelinelog_path = [main_path, '/logs/pipeline_logs'];

nr_patients = 10;
nr_controls = 14;
nr_sessions_patients = 4;
nr_sessions_controls = 2;

% Load FC matrix with data for all subjects and number of PCs used for
% reconstruction
load([DI_path, '/DI_FC_PCA/FC_PCA_recon.mat']);

% log file name
NBSstats_log_filename = [pipelinelog_path,'/step5_NBSstats-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%% Design
% Create NBS design matrices for one-sample t-test for patients and
% controls and two-sample t-test for between group comparison

managefolders([nbs_files_path, '/design'], 'create');

% patients
design = zeros(nr_patients*2, nr_patients+1);
design(:,nr_patients+1) = [ones(nr_patients,1); -1*ones(nr_patients,1)];
design(1:nr_patients,1:nr_patients) = eye(nr_patients);
design((nr_patients+1):(nr_patients*2),1:nr_patients) = eye(nr_patients);
dlmwrite([nbs_files_path, '/design/design_patients.txt'],design, 'delimiter', ',');

% figure; 
% imagesc(design)
% title('design patients')

% controls
design = zeros(nr_controls*2,nr_controls+1);
design(:,nr_controls+1) = [ones(nr_controls,1); -1*ones(nr_controls,1)];
design(1:nr_controls,1:nr_controls) = eye(nr_controls);
design((nr_controls+1):(nr_controls*2),1:nr_controls) = eye(nr_controls);
dlmwrite([nbs_files_path, '/design/design_controls.txt'],design, 'delimiter', ',');

% figure; 
% imagesc(design)
% title('design controls')

% patients vs controls
design = zeros(nr_patients+nr_controls,2);
design(1:nr_patients,1) = ones(nr_patients,1);
design((nr_patients+1):(nr_patients+nr_controls), 2) = ones(nr_controls, 1);
dlmwrite([nbs_files_path, '/design/design_patients-controls.txt'],design, 'delimiter', ',');

% figure; 
% imagesc(design)
% title('design patients vs controls')

% Log file
NBSstats_log_msg = ['------ ', 'Stored NBS design matrices'];
disp(NBSstats_log_msg);
logCustom(NBSstats_log_msg, NBSstats_log_filename)


%% Nodes - coordinates and labels
% Create .txt file with node coordinates and labels
managefolders([nbs_files_path, '/nodes'], 'create');

% read file with node information
nodes = readtable([atlas_path, '/coord_SchaeferSubCRB7100_130.txt']);

% coordinates 
dlmwrite([nbs_files_path, '/nodes/coord.txt'], nodes{:,1:3})

% labels
writetable(table(nodes{:,6}), [nbs_files_path, '/nodes/labels.txt'], 'WriteVariableNames', false);

% Log file
NBSstats_log_msg = ['------ ', 'Stored node .txt files'];
disp(NBSstats_log_msg);
logCustom(NBSstats_log_msg, NBSstats_log_filename)

%% 
components = [19, 68];

session1 = [1, 1, 1, 3, 2, 3, 41, 1, 2, 3, 4]; % 1-preictal, 2-ictal, 3-postictal, 4-interictal, 41-premenstrual, 42-midcycle
session2 = [2, 3, 4, 2, 4, 4, 42, 41, 41, 41, 42]; % 1-preictal, 2-ictal, 3-postictal, 4-interictal, 41-premenstrual, 42-midcycle
comparisons = {'preic-ic', 'preic-postic', 'preic-interic', 'postic-ic', 'ic-interic', 'postic-interic', 'prem-mid', 'preic-prem', 'ic-prem', 'postic-prem', 'interic-mid'};

methods = {'Extent','Intensity'};
contrasts = [1, -1];
thresholds = [2, 3.1, 4, 5, 6];

exchange_block_controls = '[1 2 3 4 5 6 7 8 9 10 11 12 13 14 1 2 3 4 5 6 7 8 9 10 11 12 13 14]'; 
exchange_block_patients = '[1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10]'; 
exchange_block_patients_controls = ' '; 

for comp = 1:length(components)
    
    nComp = components(comp);

    group_data = FC_PCA_recon(:,:,:,nComp); % get FC matrix for all subjects for this number of PCs 

    managefolders([nbs_stats_path, '/stats-', num2str(nComp), 'PCs'], 'create');

    %% Build matrices with concatenated FC data for each comparison

    for c = 1:length(comparisons)

        if c==7
           data1 = group_data(:,:,session1(c):nr_sessions_controls:(nr_sessions_patients*nr_patients+nr_sessions_controls*nr_controls)); 
           data2 = group_data(:,:,session2(c):nr_sessions_controls:(nr_sessions_patients*nr_patients+nr_sessions_controls*nr_controls)); 

        elseif c < 7
           data1 = group_data(:,:,session1(c):nr_sessions_patients:(nr_sessions_patients*nr_patients));
           data2 = group_data(:,:,session2(c):nr_sessions_patients:(nr_sessions_patients*nr_patients));

        elseif c > 7 
           data1 = group_data(:,:,session1(c):nr_sessions_patients:(nr_sessions_patients*nr_patients));
           data2 = group_data(:,:,session2(c):nr_sessions_controls:(nr_sessions_patients*nr_patients+nr_sessions_controls*nr_controls));
        end
        data = cat(3, data1, data2);

        managefolders([nbs_files_path, '/FCdata-',num2str(nComp), 'PCs'], 'create');
        save([nbs_files_path, '/FCdata-',num2str(nComp),'PCs/FCdata-',num2str(nComp), 'PCs_' comparisons{c}], 'data');
    end

    % Run NBS for each comparison
    for m = 1:length(methods)
        method = methods{m};
        for t = 1:length(thresholds)
            cluster_threshold = num2str(thresholds(t));
            for a = 1:length(contrasts)
                contrast_type = contrasts(a);
                
                contrast_controls =['[0 0 0 0 0 0 0 0 0 0 0 0 0 0 ', num2str(contrast_type),']'];               
                contrast_patients = ['[0 0 0 0 0 0 0 0 0 0 ', num2str(contrast_type),']'];               
                contrast_patients_controls = ['[',num2str(contrast_type),',', num2str(-1*contrast_type),']'];

                for c = 1:length(comparisons)
                    comparison = comparisons{c};

                    if c < 7
                        group = 'patients';
                        contrast_val = contrast_patients;
                        exchange_block = exchange_block_patients;
                    elseif c == 7
                        group = 'controls';
                        contrast_val = contrast_controls;
                        exchange_block = exchange_block_controls;
                    elseif c > 7
                        group = 'patients-controls';
                        contrast_val = contrast_patients_controls;
                        exchange_block = exchange_block_patients_controls;
                     end

                    UI.method.ui='Run NBS'; 
                    UI.test.ui='t-test';
                    UI.size.ui= method; % Extent
                    UI.thresh.ui= cluster_threshold;
                    UI.perms.ui='5000';
                    UI.alpha.ui='0.05';
                    UI.contrast.ui= contrast_val;
                    UI.design.ui= [nbs_files_path, '/design/design_', group,'.txt'];
                    UI.exchange.ui= exchange_block;
                    UI.matrices.ui=[nbs_files_path, '/FCdata-', num2str(nComp), 'PCs/FCdata-',num2str(nComp), 'PCs_', comparison,'.mat'];
                    UI.node_coor.ui= [nbs_files_path, '/nodes/coord.txt'];                         
                    UI.node_label.ui= [nbs_files_path, '/nodes/labels.txt'];   

                    NBSrun(UI, 'a')

                    global nbs;
                    if ~isempty(nbs.NBS.con_mat)
                        save([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(thresholds(t)))], 'nbs');
                    end

                    close all
                    
                    % Log file
                    NBSstats_log_msg = ['------ ', 'Computed NBS for: ', num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(thresholds(t)))];
                    disp(NBSstats_log_msg);
                    logCustom(NBSstats_log_msg, NBSstats_log_filename)
                   
                end % comparisons
                
            end % contrasts 
            
        end % thresholds
        
    end % methods
    
end % components