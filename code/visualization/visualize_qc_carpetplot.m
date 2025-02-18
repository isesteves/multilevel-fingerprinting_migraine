% for matlab version 2016; 
% it needs this toolbox: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

addpath('/home/iesteves/MIG_N2Treat_fMRI/NIfTI_20140122')

config_path = '/home/iesteves/FC/code/config';

config_detailts = dir('/home/iesteves/FC/code/config/step1_parcellation_config.mat');
disp(['--- Step 1 - Config file last modified on ', datestr(config_detailts.date)]);

% Load configuration settings
config = load([config_path, '/step1_parcellation_config.mat']);
config = config.config;  % Extract the structure from the loaded file

% Access script-specific parameters

preprocessed_func_path = config.in_path;
subjects = config.subjects;
sessions = config.sessions;
cleanup_list = config.pipelines;
atlas_list = config.atlases;

plt = config.visualize_carpetplot.plt;

error_log_path = config.visualize_carpetplot.error_log_path;
resp_path = config.visualize_carpetplot.resp_path;
fig_path = config.visualize_carpetplot.fig_path;
rois_path = config.visualize_carpetplot.rois_path;
task = config.visualize_carpetplot.task;
mask_types = config.visualize_carpetplot.mask_types;
TR = config.visualize_carpetplot.TR;

errorLogFileName = [error_log_path, '/visualize_carpetplot-errorlog_', datestr(now, 'yyyymmdd-HHMMSS')];

for ses = 1:length(sessions)
    session = sessions{ses};
    session_subjects = subjects{1, ses};
    
    for s = 1:length(session_subjects)
        
        subject = session_subjects{s};   

        for p = 1:length(cleanup_list)
            pipeline = cleanup_list{p};
   
            try 
            funcpath = [preprocessed_func_path, '/', subject, '/', session,'/preprocessed/filtered_func_data_', pipeline,'.nii.gz'];
            brainmaskpath = [preprocessed_func_path, '/', subject, '/', session,'/mask.nii.gz'];

            %functional image
            gunzip(funcpath)
            file = load_untouch_nii(funcpath(1:end-3));
            delete([funcpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;

            %brain mask
            gunzip(brainmaskpath)
            file_brainmask = load_untouch_nii(brainmaskpath(1:end-3));
            delete([brainmaskpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;

            %%
            func_img = file.img;
            n_x = size(func_img, 1); 
            n_y = size(func_img, 2); 
            n_z = size(func_img, 3); 
            nr_vol = size(func_img, 4);

            brainmask = file_brainmask.img;

            vol_brainmask = repmat(brainmask, [1 1 1 nr_vol]);

            bmfunc_img = func_img;
            bmfunc_img(~vol_brainmask) = 0; 

            masks_data = cell(1,3); 
            for m = 1:length(mask_types)

                mask_type = mask_types{m};

                if any(strcmp(mask_type, {'GM', 'WM'}))
                    maskpath =  [rois_path, '/', subject, '/', session,'/rois/EF_',mask_type,'_ero.nii.gz'];
                elseif strcmp(mask_type, 'CSF')
                    maskpath = [rois_path, '/', subject, '/', session,'/rois/EF_CSF_Ventricle.nii.gz'];
                end

                gunzip(maskpath)
                file_mask = load_untouch_nii(maskpath(1:end-3));
                delete([maskpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;

                mask_img = file_mask.img;

                mask = mask_img;

                % GM/WM/CSF mask: #voxels x time (TRs)
                mask_tc = reshape(mask, [n_x*n_y*n_z 1]);

                % functional image: #voxels x time (TRs)
                func_tc = reshape(func_img, [n_x*n_y*n_z nr_vol]);

                % functional image with GM/WM/CSF mask
                mfunc_tc = func_tc(logical(mask_tc),:);

                masks_data{1,m} = mfunc_tc;

            end

            %% Voxels for each mask
            gm_data = masks_data{1,1};
            wm_data = masks_data{1,2};
            csf_data = masks_data{1,3};

            %% Load Motion outliers
            mofile = [preprocessed_func_path, '/', subject, '/', session,'/mo_confound_dvars.txt'];
            
            try
            mo_data = load(mofile, '-ascii');
            catch
                mo_data = [];
               
            end
            %% Load MO metric
            dvarsfile = [preprocessed_func_path, '/', subject, '/', session, '/dvars.txt'];
            
            dvarsdata = load(dvarsfile, '-ascii');

            mo_metric = 'dvars';
            mo_metric_data = dvarsdata;
          
            %% Load respiration data
            respfile = [resp_path, '/', task, '/', session, '/', subject, '_', session, '_', task, '_resp-preproc.mat'];
            try
                respdata_struct = load(respfile);

                resp_data = respdata_struct.resp.data;
                resp_times = respdata_struct.resp.times;
            catch
                resp_data = [];
                resp_times = [];
            end
            %% Plot carpet
            figure('pos', [20 20 1000 900]);
            h = suptitle({[subject,' - ', session, ' - ', task, ' - ', pipeline], '  ', ' '});
            h.Interpreter = 'none';

            plot_carpet_all(gm_data, wm_data, csf_data, resp_data, resp_times, mo_metric, mo_metric_data, mo_data, TR)
            
            managefolders([fig_path, '/', subject, '/'], 'create');
            print([fig_path, '/', subject, '/carpetplot_', subject, '_', session, '_', task], '-dpng')
            
            catch ME
                % Handle errors and log them
                errorMessage = [subject, '_', session, '_', task, ' --- ', ME.identifier, ' - ', ME.message]; % Removes ANSI color codes
                logError(errorMessage, errorLogFileName);
                disp('An error occurred. Please check the error log for details.');
            end
        end
    end
end
