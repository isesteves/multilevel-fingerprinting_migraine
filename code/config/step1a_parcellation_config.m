%% Paths and variables
config_path = '/home/iesteves/FC/code/config/';
pipeline_step = 'step1a_parcellation';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define subjects and sessions
control_premenstrual = {'sub-control019','sub-control020','sub-control025','sub-control026','sub-control027','sub-control028','sub-control029',...
'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

control_midcycle = {'sub-control019','sub-control020','sub-control025','sub-control026','sub-control027','sub-control028','sub-control029',...
'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

patient_interictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006',  'sub-patient007', 'sub-patient008',...
'sub-patient009', 'sub-patient012', 'sub-patient013', 'sub-patient034', 'sub-patient038', 'sub-patient041',...
'sub-patient043', 'sub-patient045'};

patient_ictal = {'sub-patient001', 'sub-patient002', 'sub-patient003', 'sub-patient004', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient013', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045', 'sub-patient052'};

patient_preictal = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient013', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient043', 'sub-patient045'};

patient_postictal = {'sub-patient001', 'sub-patient002' 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
'sub-patient009', 'sub-patient012', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient043', 'sub-patient045'};

%% Combine control and patient groups for all sessions
all_control_subjects = {control_premenstrual; control_midcycle};
all_patient_subjects = {patient_preictal; patient_ictal; patient_postictal; patient_interictal};
all_subjects = {all_control_subjects{:}, all_patient_subjects{:}};

all_control_sessions = {'ses-premenstrual', 'ses-midcycle'};
all_patient_sessions = {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-interictal'};
all_sessions = [all_control_sessions, all_patient_sessions];

%% Creat configuration file
config = struct();
config.project_path = '/home/iesteves/FC';
config.in_path = '/home/iesteves/FC/data/parcellated_data';
config.atlas_path = '/home/iesteves/FC/files/atlases';
config.pipelinelog_path = '/home/iesteves/FC/logs/pipeline_logs';
config.subjects = all_subjects;
config.sessions = all_sessions;
config.pipelines = {'nocleanup_preprocessed'};%{'preprocessed_rp_mo_csf_wm_nui'};
config.atlases = {'SchaeferSubCRB7100'};
config.atlases_Nareas = 138;
config.plt = 1;
config.fig_path = '/home/iesteves/FC/figures/parcellation_nullROI';

%% Save configuration file and update archive if needed
% Save configuration file
save([config_path, '/', pipeline_step,'_config.mat'], 'config');

% Check if the config archive needs to be updated
config_tracker(config_path, pipeline_step, config);

