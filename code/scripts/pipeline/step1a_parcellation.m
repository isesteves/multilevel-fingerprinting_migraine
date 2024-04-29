%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 1A - REORGANIZATION OF PARCELS TO GET ONLY NON-ZERO SchaeferSubCRB7100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
disp('>>>>>>> [Optional] Step 1a - Reorganization of parcels to get SchaeferSubCRB7100');

config_path = '/home/iesteves/FC/code/config';

config_details = dir('/home/iesteves/FC/code/config/step1a_parcellation_config.mat');
disp(['--- Step 1a - Config file last modified on ', datestr(config_details.date)]);

% Load configuration settings
config = load([config_path, '/step1a_parcellation_config.mat']);
config = config.config;  % Extract the structure from the loaded file

% Access script-specific parameters
parcellation_inpath = config.in_path;
atlas_path = config.atlas_path;
pipelinelog_path = config.pipelinelog_path;
project_path = config.project_path;
atlas = config.atlases{1};

parcellation_log_filename = [pipelinelog_path,'/step1a_parcellation-log_', datestr(now, 'yyyymmdd-HHMMSS'),'.txt'];

ind_H2N = load([atlas_path, '/Schaefer/Schaefer_H2N_ind.mat']);
ind = ind_H2N.ind;

% AAL116 indices for subcortical and cerebellum parcels
subcortical_ind_aal = [37, 38, 41, 42, 71:78];
cerebellum_ind_aal = 91:116;

nr_areas = 138; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if ~exist([atlas_path, '/Schaefer/Schaefer2018_100Parcels_7Networks_order_new.txt'], 'file')
    run([project_path, '/code/scripts/utils/reorganize_Schaefer_pernetwork.m'])
end


% load Schaefer (reorganized by network)
atlas_Schaefer = readtable([atlas_path, '/Schaefer/Schaefer2018_100Parcels_7Networks_order_new.txt']);

% load original Schaefer centroid MNI coordinates and reorganize them by
% network
centroids = readtable([atlas_path, '/Schaefer/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv']);
order = atlas_Schaefer.Var1;
new_centroids = centroids(order, :);

% load AAL116 coordinates
atlas_aal = readtable([atlas_path,'/AAL/AAL116_coordinates_brainGraph.xlsx']); 
SubCRB_areas = [37 38 41 42 71:78 91:116];

atlas_aal.name  = cellfun(@(x) replace(x, '.', '_'), atlas_aal.name, 'UniformOutput', false);

labels = atlas_Schaefer.Var5;
labels_all = cell(nr_areas,1);
labels_all(1:100) = labels;
labels_all(101:end,:) = atlas_aal{SubCRB_areas, 2};

coordinates_all = zeros(nr_areas,3);
coordinates_all(1:100,:) = new_centroids{:,3:5};
coordinates_all(101:end,:) = atlas_aal{SubCRB_areas, 3:5};

areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 138 9];

colors_all = zeros(nr_areas, 1);
for k = 1:size(areas, 1)
    colors_all(areas(k,1):areas(k,2)) = areas(k,3);
end

size_all = 3*ones(nr_areas, 1);

T = table(coordinates_all, colors_all, size_all, labels_all);
writetable(T,[atlas_path, '/coord_SchaeferSubCRB7100_', num2str(nr_areas), ''],'Delimiter',' ', 'WriteVariableNames', 0);


%% Load Schaefer parcellated files (hemisphere), reorganize by network; add AAL116 subcortical + cerebellum
% Load all Schaefer files - organized by hemisphere 
files = dir([parcellation_inpath, '/sub-*/ses-*/BOLDpreprocessed_*SchaeferH7100*']);
filenames = {files(:).name};
filefolders = {files(:).folder};

for f = 1:length(filenames)
    
    % load parcellated BOLD data with Schaefer organized by hemisphere
    filename = filenames{f};
    load([filefolders{f}, '/', filename]);
    
    % reorganize parcellated BOLD data by network using the pre-computed
    % order
    newBOLD = all_av(ind,:);
    
    % load percentage included of each parcel - Schaefer organized by
    % hemisphere
    load([filefolders{f}, '/percent-', filename(5:end)]);
    
    % reorganized percentage included of each parcel - Schaefer by
    % network
    newpercent = percentage_included(ind, :);
    
    % load correspondind parcellated BOLD using AAL116
    file_aux = strsplit(filename, '-');
    pipeline = file_aux{1,1};   
    aal = load([filefolders{f}, '/', pipeline, '-parcellatedAAL116woSchaefer.mat']);
    aal_data = aal.all_av;
    
    % get subcortical and cerebellum parcels
    subcortical = aal_data(subcortical_ind_aal,:);
    cerebellum = aal_data(cerebellum_ind_aal,:);
    
    % join BOLD parcellated data using Schaefer reorganized by network,
    % AAL116 subcortical and AAL116 cerebellum
    all_av = [newBOLD; subcortical; cerebellum];
    
    aal_percent = load([filefolders{f}, '/percent-', pipeline(5:end), '-parcellatedAAL116woSchaefer.mat']);
    aal_percentage_included = aal_percent.percentage_included;
    
    % percentage included of each parcel of AAL116 subcortical and
    % cerebellum
    subcortical_percent = aal_percentage_included(subcortical_ind_aal,1);
    cerebellum_percent = aal_percentage_included(cerebellum_ind_aal,:);
    
    % join parcel percentage included using Schaefer reorganized by network,
    % AAL116 subcortical and AAL116 cerebellum
    percentage_included = [newpercent; subcortical_percent; cerebellum_percent];
    
    % save parcellated BOLD and percentage included of each parcel for
    % SchaeferSubCRB7100 combination
    save([filefolders{f}, '/', filename(1:end-9), 'SubCRB', filename(end-7:end)], 'all_av')
    save([filefolders{f}, '/percent-', filename(5:end-9), 'SubCRB', filename(end-7:end)], 'percentage_included')
 
    filedetails_aux = strsplit(filefolders{f}, '/');
    subject = filedetails_aux{end-1};
    session = filedetails_aux{end};
    parcellation_log_msg = ['------ [Step 1a] ', atlas,' Parcellation for: ', subject, '_', session,' with pipeline ', pipeline];
    disp(parcellation_log_msg);
    logCustom(parcellation_log_msg, parcellation_log_filename)

end


%% 
disp('=== [Visualization] - Parcels to exclude ');
run([project_path, '/code/visualization/visualize_parcellation_nullROI.m']);

exclude_areas = load([atlas_path, '/parcellation_nullROI_parcel2exclude_SchaeferSubCRB7100_idx']); %[115, 116, 123, 124, 125, 126, 127, 128];
include_areas = setdiff(1:nr_areas, exclude_areas.parcel2exclude);
include_Nareas = length(include_areas);

T_included = table(coordinates_all(include_areas,:), colors_all(include_areas,:), size_all(include_areas,:), labels_all(include_areas,:));
writetable(T_included,[atlas_path, '/coord_SchaeferSubCRB7100_', num2str(include_Nareas)],'Delimiter',' ', 'WriteVariableNames', 0)  
