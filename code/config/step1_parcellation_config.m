%% Paths and variables
config_path = '/home/iesteves/FC/code/config/';
pipeline_step = 'step1_parcellation';


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
config.in_path = '/home/mig_n2treatdata/derivatives/func-preproc/func-task_rest';
config.out_path = '/home/iesteves/FC/data/parcellated_data';
config.atlas_path = '/home/iesteves/FC/files/atlases';
config.pipelinelog_path = '/home/iesteves/FC/logs/pipeline_logs';
config.subjects = all_subjects;
config.sessions = all_sessions;
config.pipelines = {'nocleanup_preprocessed', 'preprocessed_rp_mo_csf_wm_nui'};
config.atlases =  {'AAL90'};%{'AAL116woSchaefer', 'SchaeferH7100', 'HOcort', 'Desikan'};
config.plt = 1;

config.visualize_carpetplot.plt = 1;
config.visualize_carpetplot.error_log_path = '/home/iesteves/FC/logs/error_logs';
config.visualize_carpetplot.resp_path = '/home/iesteves/Physio/derivatives/resp-preproc';
config.visualize_carpetplot.fig_path = '/home/iesteves/FC/figures/qc_carpetplot';
config.visualize_carpetplot.rois_path = '/home/mig_n2treatdata/derivatives/func-preproc/func-task_rest';
config.visualize_carpetplot.task = 'task-rest';
config.visualize_carpetplot.mask_types = {'GM', 'WM', 'CSF'};
config.visualize_carpetplot.TR = 1.26;

%% Save configuration file and update archive if needed
% Save configuration file
save([config_path, '/', pipeline_step, '_config.mat'], 'config');

% Check if the config archive needs to be updated
config_tracker(config_path, pipeline_step, config);

