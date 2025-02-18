% for matlab version 2016; 
% it needs this toolbox: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

addpath('/home/iesteves/MIG_N2Treat_fMRI/NIfTI_20140122')
path = '/home/iesteves/eeg-fmri';
resppath = '/home/iesteves/Physio/derivatives/resp-complete';
figurepath = '/home/iesteves/MIG_N2Treat_fMRI/Figures/'; 


dataset = 'MIGN2TREAT';
task = 'task-rest';
session = 'ses-midcycle';
mask_types = {'GM', 'WM', 'CSF'};
funcpipeline = 'nocleanup_preprocessed';
TR = 1.26;

subject = 'sub-control048';


funcpath = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/filtered_func_data_', funcpipeline,'.nii.gz'];
brainmaskpath = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mask.nii.gz'];

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
        maskpath =  [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/masks/EF_',mask_type,'_ero.nii.gz'];
    elseif strcmp(mask_type, 'CSF')
        maskpath =  [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/masks/EF_CSF_Ventricle.nii.gz'];
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
mofile = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mo_confound_dvars.txt'];
mo_data = load(mofile, '-ascii');

%% Load MO metric
dvarsfile = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mc/dvars.txt'];
dvarsdata = load(dvarsfile, '-ascii');

mo_metric = 'dvars';
mo_metric_data = dvarsdata;

%% Load respiration data
respfile = [resppath, '/', task, '/', session, '/', subject, '_', session, '_', task, '_resp-complete.mat'];
respdata_struct = load(respfile);

resp_data = respdata_struct.resp.valid.data;
resp_times = respdata_struct.resp.valid.times;

%% Plot carpet
figure('pos', [20 20 1000 900]);
h = suptitle({[subject,' - ', session, ' - ', task, ' - ', funcpipeline], '  ', ' '});
h.Interpreter = 'none';

plot_carpet_all(gm_data, wm_data, csf_data, resp_data, resp_times, mo_metric, mo_metric_data, mo_data, TR)



