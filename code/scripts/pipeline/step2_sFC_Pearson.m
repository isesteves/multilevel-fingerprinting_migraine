%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 2 - COMPUTATION OF sFC USING PEARSON CORRELATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

disp(['>>>>>>> Step 2 - Computation of sFC using Pearson correlation']);

config_path = '/home/iesteves/FC/code/config';

config_details = dir('/home/iesteves/FC/code/config/step2_sFC_Pearson_config.mat');
disp(['--- Step 2 - Config file last modified on ', datestr(config_details.date)]);

% Load configuration settings
config = load([config_path, '/step2_sFC_Pearson_config.mat']);
config = config.config;  % Extract the structure from the loaded file

% Access script-specific parameters
sFC_inpath = config.in_path;
sFC_outpath = config.out_path;
pipelinelog_path = config.pipelinelog_path;

subjects = config.subjects;
sessions = config.sessions;
cleanup_list = config.pipelines;
atlas_list = config.atlases;
include_ind = config.include_areas;
N_areas_list = config.include_Nareas;
N_subjects = config.Nsubjects;
TR = config.TR;
N_windows = config.nr_windows;
order_bysubject = config.order_bysubject; 
sFC_log_filename = [pipelinelog_path,'/step2_sFC_Pearson-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%plt = config.plt;

% Bandpass filter settings
fnq = 1/(2*TR);               % Nyquist frequency
flp = .01;                    % Lowpass frequency of filter (Hz)
fhi = .1;                     % Highest BOLD frequency considered (Hz)

Wn = [flp/fnq fhi/fnq];       % Butterworth bandpass non-dimensional frequency
k = 2;                        % 2nd order butterworth filter
[bfilt,afilt] = butter(k,Wn); % Construct the filter

%%
% Loop over atlases
for a = 1:length(atlas_list)
    atlas = atlas_list{a};
    N_areas = N_areas_list(a);
    
    % Loop over cleanup pipelines
    for p = 1:length(cleanup_list)        
        pipeline = cleanup_list{p};
         
        % List of filenames for all subjects and sessions
        auxfilenames = cell(length(sessions),1);
        for k = 1:length(sessions)
            kfiles = subjects{k};
            files = cellfun(@(x) [sFC_inpath,'/',x,'/',sessions{k},'/BOLD', pipeline, '-parcellated', atlas, '.mat'],kfiles,'UniformOutput',false);
            auxfilenames{k,1} = files; 
        end
        filenames = [auxfilenames{:}];
        
               
        %% Obtain FC matrices using Pearson correlation
        % Preallocate variables to save FC patterns and associated information
        V1_all = zeros(N_windows*N_subjects,N_areas); % All leading Eigenvectors
        Var_Eig = zeros(N_subjects,N_windows); % Saves the variance explained by the leading Eigenvectors
        FC_all = zeros(N_areas,N_areas,N_subjects);
    
        % Loop over subjects and sessions
        for s = 1:length(filenames)
            
            % Specify input data
            data_in = filenames{s};

            % Load the BOLD matrix (NxT) from each subject
            load_BOLD = load(data_in);
            BOLD_aux = load_BOLD.all_av;
            BOLD = BOLD_aux(include_ind, :);
            
            %De-mean BOLD
            demeanedBOLD = BOLD-mean(BOLD, 2);
            
            % Filter the demeaned BOLD data for each seed (area)   
            BOLD_filt = zeros(size(BOLD)); 
            for seed = 1:N_areas
                signal_filt = filtfilt(bfilt,afilt,demeanedBOLD(seed,:));
                BOLD_filt(seed,:) = signal_filt;
            end
            
            % Compute Pearson correlation (FC)
            FC = zeros(N_areas);
            for n = 1:N_areas
                for k = 1:N_areas
                    FC(n,k) = corr(BOLD_filt(n, :)', BOLD_filt(k, :)', 'Type', 'Pearson');     
                end    
            end
            FC_all(:, :, s) = FC;   
            
            %%% From LEiDA pipeline, to match dFC
            % Get the leading Eigenvector
            [eVec,eigVal] = eigs(FC);
            eVal = diag(eigVal);
            [val1, i_vec_1] = max(eVal);
            V1 = eVec(:,i_vec_1);
            % Make sure the largest component is negative 
            % This step is important because the same eigenvector can 
            % be returned either as V or its symmetric -V and we need
            % to make sure it is always the same (so we choose always
            % the most negative one)
            if sum(V1)<0
                V1 = -V1;
            end
            V1_all(s,:) = V1;

            % Compute the variance explained by the leading Eigenvector
            Var_Eig(s) = val1/sum(eVal);
            
            % Log file
            sFC_log_msg = ['------ ', 'sFC Pearson for: ', data_in];
            disp(sFC_log_msg);
            logCustom(sFC_log_msg, sFC_log_filename)

          
        end % filenames: subjects and sessions
        
        % Create output path if not existent and define output file name
        managefolders(sFC_outpath, 'create');
        
        % Save FC data for the whole sample (organized by session type)
        save([sFC_outpath, '/Pearson-results-',atlas,'-', pipeline, '.mat'], 'V1_all','Var_Eig','N_subjects', 'filenames', 'N_areas','FC_all');
        
        % Save FC data for the whole sample (organized by subject)
        filenames_bysubject = filenames(order_bysubject);
        FC_all_bysubject = FC_all(:,:,order_bysubject);
        V1_all_bysubject = V1_all(order_bysubject);
        Var_Eig_bysubject = Var_Eig(order_bysubject); 
        save([sFC_outpath, '/Pearson-results-',atlas,'-', pipeline, '-bysubject.mat'], 'V1_all_bysubject', 'Var_Eig_bysubject', 'N_subjects', 'filenames_bysubject', 'N_areas', 'FC_all_bysubject');

        
    end % cleanup pipelines
end % atlases


