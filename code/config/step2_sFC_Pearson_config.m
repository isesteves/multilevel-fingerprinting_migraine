%% Paths and variables
config_path = '/home/iesteves/FC/code/config/';
pipeline_step = 'step2_sFC_Pearson';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define subjects and sessions - consider only complete subjects
patient_preictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

patient_ictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

patient_postictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

patient_interictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

control_premenstrual = {'sub-control019','sub-control020','sub-control026','sub-control027','sub-control028','sub-control029',...
'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

control_midcycle = {'sub-control019','sub-control020', 'sub-control026','sub-control027','sub-control028','sub-control029',...
'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

%% Combine control and patient groups for all sessions
all_patient_subjects = {patient_preictal; patient_ictal; patient_postictal; patient_interictal};
all_control_subjects = {control_premenstrual; control_midcycle};
all_subjects = {all_patient_subjects{:}, all_control_subjects{:}};

all_patient_sessions = {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-interictal'};
all_control_sessions = {'ses-premenstrual', 'ses-midcycle'};
all_sessions = [all_patient_sessions, all_control_sessions];

%% Creat configuration file
config = struct();
% config.project_path = '/home/iesteves/FC';
config.in_path = '/home/iesteves/FC/data/parcellated_data';
config.out_path = '/home/iesteves/FC/data/results/sFC_Pearson/MIGcomplete-HCcomplete_SchaeferSubCRB7100_nonzero-50';
% config.atlas_path = '/home/iesteves/FC/files/atlases';
config.pipelinelog_path = '/home/iesteves/FC/logs/pipeline_logs';
config.subjects = all_subjects;
config.Nsubjects = length([all_subjects{:}]);
config.sessions = all_sessions;
config.pipelines = {'preprocessed_rp_mo_csf_wm_nui'};
config.atlases = {'SchaeferSubCRB7100'};
config.atlases_Nareas = 138;
config.TR = 1.26;
config.nr_windows = 1;
config.exclude_areas = load('/home/iesteves/FC/files/atlases/parcellation_nullROI_parcel2exclude_SchaeferSubCRB7100_idx'); %[115, 116, 123, 124, 125, 126, 127, 128];
config.include_areas = setdiff(1:config.atlases_Nareas, config.exclude_areas.parcel2exclude);
config.include_Nareas = length(config.include_areas);
config.order_bysubject = [1:10:31 2:10:32 3:10:33 4:10:34 5:10:35 6:10:36 7:10:37 8:10:38 9:10:39 10:10:40 ...
    41 55 42 56 43 57 44 58 45 59 46 60 47 61 48 62 49 63 50 64 51 65 52 66 53 67 54 68];%  
% config.plt = 1;
% config.fig_path = '';

%% Save configuration file and update archive if needed
% Save configuration file
save([config_path, '/', pipeline_step,'_config.mat'], 'config');

% Check if the config archive needs to be updated
config_tracker(config_path, pipeline_step, config);

