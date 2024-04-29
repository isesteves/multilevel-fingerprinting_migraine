%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 1 - PARCELLATION OF PREPROCESSED FUNCTIONAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

disp(['>>>>>>> Step 1 - Parcellation of preprocessed functional data']);

config_path = '/home/iesteves/FC/code/config';

config_detailts = dir('/home/iesteves/FC/code/config/step1_parcellation_config.mat');
disp(['--- Step 1 - Config file last modified on ', datestr(config_detailts.date)]);

% Load configuration settings
config = load([config_path, '/step1_parcellation_config.mat']);
config = config.config;  % Extract the structure from the loaded file

% Access script-specific parameters
parcellation_inpath = config.in_path;
parcellation_outpath = config.out_path;
parcellation_atlas_path = config.atlas_path;
pipelinelog_path = config.pipelinelog_path;

subjects = config.subjects;
sessions = config.sessions;
cleanup_list = config.pipelines;
atlas_list = config.atlases;

plt = config.plt;

parcellation_log_filename = [pipelinelog_path,'/step1_parcellation-log_', datestr(now, 'yyyymmdd-HHMMSS'),'.txt'];

for ses = 1:length(sessions)
    session = sessions{ses};
    session_subjects = subjects{1, ses};
    
    for s = 1:length(session_subjects)
        
        subject = session_subjects{s};   

        for p = 1:length(cleanup_list)
            pipeline = cleanup_list{p};

            % Specify input path
            img_path = [parcellation_inpath,'/', subject, '/', session, '/preprocessed/filtered_func_data_', pipeline, '2standard.nii.gz'];

            % Load functional image in nii.gz format
            gunzip(img_path);
            data_load = load_untouch_nii(img_path(1:end-3));
            delete(img_path(1:end-3)); % delete unzipped image after assigning it to a variable (saves space)

            funcimg = data_load.img;

            for a = 1:length(atlas_list)
                atlas = atlas_list{a};
                outpath_av_nonzero = [parcellation_outpath, '/', subject, '/', session, '/BOLD', pipeline, '-parcellated', atlas, '.mat'];
                outpath_av_all = [parcellation_outpath, '/', subject, '/', session, '/BOLDall-', pipeline, '-parcellated', atlas, '.mat'];
                outpath_percent = [parcellation_outpath, '/', subject, '/', session, '/percent-', pipeline, '-parcellated', atlas, '.mat'];

                parcellation_log_msg = ['------ ', atlas,' Parcellation for: ', subject, '_', session,' with pipeline ', pipeline];
                disp(parcellation_log_msg);
                logCustom(parcellation_log_msg, parcellation_log_filename)

                [all_av_allvoxels, all_av, percentage_included] = parcellation(parcellation_atlas_path, atlas, funcimg, plt);

                out_folder = [parcellation_outpath, '/', subject, '/', session];
                managefolders(out_folder, 'create')
                save(outpath_av_nonzero,'all_av','-mat');
                save(outpath_av_all,'all_av_allvoxels','-mat');
                save(outpath_percent,'percentage_included','-mat');
            end
        end
    end
end
