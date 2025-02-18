function [all_av_allvoxels, all_av_nonzero, percentage_included] = parcellation(atlasdir, atlas_name, funcimg, plt)
%PARCELLATION Parcellate functional image based on specified atlas.
%
%   [all_av_allvoxels, all_av_nonzero, percentage_included] = PARCELLATION(atlasdir, atlas_name, funcimg, plt)
%
%   INPUTS:
%   - atlasdir: Directory containing atlas files.
%   - atlas_name: Name of the atlas to use ('AAL90', 'AAL116', 'Desikan', 'HO', 'HOcort', 'HOsubcort', 'SchaeferH7100', 'SchaeferH17100', 'CDK50', 'CDK100', 'CDK500').
%   - funcimg: Functional image data.
%   - plt: Boolean flag indicating whether to plot the atlas (true or false).
%
%   OUTPUTS:
%   - all_av_allvoxels: Matrix containing mean intensity values for all
%   voxels in each region (dim = #parcels x #volumes).
%   - all_av_nonzero: Matrix containing mean intensity values for non-zero
%   voxels in each region (dim = #parcels x #volumes).
%   - percentage_included: Percentage of included voxels for each region (dim = #parcels x 1).
%
%   EXAMPLE USAGE:
%   [all_av_allvoxels, all_av_nonzero, percentage_included] = parcellation('/path/to/atlas/directory', 'AAL90', funcimg_data, true);
%
%   NOTES:
%   - Make sure to have the required atlas files in the specified directory.
%   - It relies on "Jimmy Shen (2023). Tools for NIfTI and ANALYZE image 
% (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), 
% MATLAB Central File Exchange. Retrieved September 29, 2023.
%   - Set plt to true if you want to visualize the atlas.
%
%  Note: All atlases and preprocessed data are registered in the MNI152 space
%  (91x109x91, 2x2x2 mm^3).
%
%   [Ines Esteves]

    cd(atlasdir)
    
    % load the selected atlas
    switch atlas_name
        case 'AAL90'
            atlas_load = load('AAL/AAL90_2mm_woSchaefer.mat'); 
            atlas = atlas_load.AAL90; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'AAL116'
            gunzip('AAL/AAL116_MNI152-2mm.nii.gz')
            atlas_load = load_nii('AAL/AAL116_MNI152-2mm.nii'); % extracted as 'AAL_space-MNI152NLin6_res-2x2x2.nii.gz'
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'AAL116woSchaefer'
            gunzip('AAL/AAL116_2mm_woSchaefer.nii.gz')
            atlas_load = load_nii('AAL/AAL116_2mm_woSchaefer.nii'); % extracted as 'AAL_space-MNI152NLin6_res-2x2x2.nii.gz'
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0   
                       
        case 'Desikan'
            atlas_load = load('Desikan/adjustedDesikan.mat'); 
            atlas = atlas_load.adjustedDesikan; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'HO'           
            atlas_load = load('Harvard-Oxford/HO.mat'); % HO cort+subcort recoded atlas
            atlas = atlas_load.newHO; % extract array from structure
            areas = 63;
            
        case 'HOcort' % Harvard-Oxford cortical atlas
            gunzip('Harvard-Oxford/HO-cortical_MNI152-2mm.nii.gz');
            atlas_load = load_untouch_nii('Harvard-Oxford/HO-cortical_MNI152-2mm.nii'); % extracted as 'HarvardOxfordcort-maxprob-thr25_space-MNI152NLin6_res-2x2x2.nii.gz' 
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'HOsubcort'
            gunzip('Harvard-Oxford/HO-subcortical_MNI152-2mm.nii.gz')
            atlas_load = load_untouch_nii('Harvard-Oxford/HO-subcortical_MNI152-2mm.nii'); % extracted as 'HarvardOxfordsub-maxprob-thr25_space-MNI152NLin6_label_all_res-2x2x2.nii.gz'
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
         case 'SchaeferH7100'
            atlas_load = load_nii('Schaefer/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii');
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'SchaeferH17100'
            gunzip('Schaefer/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz')
            atlas_load = load_nii('Schaefer/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii');
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'CDK50'
            gunzip('Craddock/Craddock-K50_MNI152-2mm.nii.gz')
            atlas_load = load_nii('Craddock/Craddock-K50_MNI152-2mm.nii');
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
            
        case 'CDK100'
            gunzip('Craddock/Craddock-K100_MNI152-2mm.nii.gz')
            atlas_load = load_nii('Craddock/Craddock-K100_MNI152-2mm.nii');
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
        
        case 'CDK500'
            gunzip('Craddock/Craddock-K500_MNI152-2mm.nii.gz')
            atlas_load = load_nii('Craddock/Craddock-K500_MNI152-2mm.nii');
            atlas = atlas_load.img; % extract array from structure
            areas = length(unique(atlas))-1; % excludes the intensity equal to 0
        
        otherwise
            disp('The specified atlas is not available.')
            
    end 
    
    % Plot atlas if plt is true
    if plt
        slice = 45;
        imshow(atlas(:,:,slice), [], 'InitialMagnification', 500);
        title([atlas_name, ' atlas']);
    end

    % Number of regions in the atlas
    disp(['This atlas parcellates the brain into ', num2str(areas), ' regions.']);

    % Create matrices to store the parcellated image
    timepoints = size(funcimg, 4);
    all_ind = unique(atlas); % all indexes or intensities
    all_ind(all_ind==0) = []; % excludes the intensity equal to 0
    all_ind2 = all_ind(1:areas,:);
    all_av_allvoxels = zeros(length(all_ind2),timepoints);
    all_av_nonzero = zeros(length(all_ind2),timepoints);
    percentage_included = zeros(areas, 1);
    
    % Parcellate functional image
    aux = funcimg;
    BOLD2D = reshape(aux,size(aux,1)*size(aux,2)*size(aux,3),size(aux,4)); % changes the matrix dimensions, lin
    for i = 1:length(all_ind2)
        ind = find(atlas == all_ind2(i));
        all_av_allvoxels(i,:) = mean(BOLD2D(ind,:));
        
        zero_points = BOLD2D(ind,:)== 0;
        vol_zero = sum(zero_points,2);
        vol_nonzero_voxel = find(vol_zero < 0.01*timepoints);
        all_av_nonzero(i,:) = mean(BOLD2D(ind(vol_nonzero_voxel),:));
        
        percentage_included(i) = 100*length(vol_nonzero_voxel)/length(ind);
        
        % Plot mean intensity for all and only non-zero voxels figure
%         figure; 
%         plot(mean(BOLD2D(ind,:))-mean(mean(BOLD2D(ind,:))))
%         hold on;
%         plot(mean(BOLD2D(ind(vol_nonzero_voxel),:))-mean(mean(BOLD2D(ind(vol_nonzero_voxel),:))))
%         xlabel('TR')
%         ylabel('Demeaned intensity')
%         legend('All voxels', 'Non-zero voxels')
%         title(['Region:', num2str(i), ' Percentage included: ', num2str(round(100*length(vol_nonzero_voxel)/length(ind))), ' %'])
    end
    
end

