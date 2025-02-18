addpath(genpath('/home/iesteves/toolboxes/NIfTI_20140122'))

img_Schaefer = load_untouch_nii('/home/iesteves/FC/files/atlases/Schaefer/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii');
img_AAL = load_untouch_nii('/home/iesteves/FC/files/atlases/AAL/AAL116_2mm_woSchaefer.nii');
%gunzip('/home/mig_n2treatdata/derivatives/func-preproc/func-task_rest/sub-control027/ses-midcycle/preprocessed/filtered_func_data_preprocessed_rp_mo_csf_wm_nui2standard.nii.gz')
%img_fmri_file = load_untouch_nii('/home/mig_n2treatdata/derivatives/func-preproc/func-task_rest/sub-control027/ses-midcycle/preprocessed/filtered_func_data_preprocessed_rp_mo_csf_wm_nui2standard.nii');


% preprocessed fMRI image
img_fmri = img_fmri_file.img;

%% Schaefer and AAL
slice_Schaefer = 40;
atlas_Schaefer = img_Schaefer.img;
slice_AAL = 37;
atlas_AAL = img_AAL.img;

mask_AAL = ismember(atlas_AAL, [37, 38, 41, 42, 71:78, 91, 92, 95:100, 107:116]);
im_AAL = zeros(size(atlas_AAL));
im_AAL(mask_AAL) = atlas_AAL(mask_AAL);
im_AAL_black = ones(size(atlas_AAL));
zeros_AAL = zeros(size(atlas_AAL));
im_AAL_black(mask_AAL) = zeros_AAL(mask_AAL);

% figure('pos', [50 50 1100 400]);
% subplot(1,2,1)
% imshow(imrotate(atlas_Schaefer(:,:,slice_Schaefer), 90, 'bilinear', 'crop'), parula)
% set(gca, 'XTick', [], 'YTick', [])
% 
% subplot(1,2,2)
% imshow(imrotate(im_AAL(:,:,slice_AAL), 90, 'bilinear', 'crop'), [0 0 1; 0 1 0; 1 0 0; 0.5 0.5 0.5; 0.3 0.4 0.6], 'InitialMagnification', 'fit')
% colormap(parula)
% set(gca, 'XTick', [], 'YTick', [])
% hold on;
% cmap = [0 0 0; 0 0 0];
% h = imshow(ind2rgb(imrotate(im_AAL_black(:,:,slice_AAL),90, 'bilinear', 'crop')+1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
% set(h, 'AlphaData', imrotate(im_AAL_black(:,:,slice_AAL),90, 'bilinear', 'crop')); % Set transparency based on mask values
% hold off;
% set(gca, 'FontSize',14)
% 
% uniqueValues = unique(im_AAL(im_AAL~=0));
% [valueMap, ~, remappedTensor] = unique(im_AAL, 'rows');
% remappedTensor(remappedTensor ~= 0) = 1:numel(uniqueValues);

% imagesc
figure('pos', [50 50 1100 400]);
subplot(1,2,1)
imagesc(imrotate(atlas_Schaefer(:,:,slice_Schaefer), 90, 'bilinear', 'crop'))
set(gca, 'XTick', [], 'YTick', [])

ax = subplot(1,2,2);
imagesc(imrotate(im_AAL(:,:,slice_AAL), 90, 'bilinear', 'crop'))
set(gca, 'XTick', [], 'YTick', [])
print('/home/iesteves/FC/figures/atlases/atlases', '-dpng')

% grouped
figure('pos', [50 50 1000 850]);
imagesc(imrotate(atlas_Schaefer(:,:,slice_Schaefer)+im_AAL(:,:,slice_Schaefer), 90, 'bilinear', 'crop'))
set(gca, 'XTick', [], 'YTick', [])
% print('/home/iesteves/FC/figures/atlases/atlases_grouped', '-dpng')

figure('pos', [50 50 1000 850]);
imagesc(imrotate(img_fmri(:,:,slice_Schaefer), 90, 'bilinear', 'crop'))
set(gca, 'XTick', [], 'YTick', [])
colormap(gray)
% print('/home/iesteves/FC/figures/atlases/preprocessed_fMRI', '-dpng')

%%
im_CRB = atlas_AAL>91;

figure; 
imagesc(img_fmri(:,:,10))
colormap(gray)
hold on;
imAlpha=ones(91,109);
imAlpha(im_CRB(:,:,10)==0)=0;
imagesc(double(im_CRB(:,:,10)),'AlphaData', imAlpha)
colormap(parula)

slice_CRB = 15;
figure;
imshow(img_fmri(:,:,slice_CRB,60), []);
colormap(gray);
% Create a colormap with transparency and red
cmap = [0 0 0; 1 0 0]; % Transparent for 0, Red for 1

% Overlay colored mask with transparency
hold on;
h = imshow(ind2rgb(im_CRB(:,:,slice_CRB) + 1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
set(h, 'AlphaData', im_CRB(:,:,slice_CRB)); % Set transparency based on mask values

im_subcortical = (atlas_AAL==37)+(atlas_AAL==38)+(atlas_AAL==41)+(atlas_AAL==42)+(atlas_AAL==71)+(atlas_AAL==72)+(atlas_AAL==73)+(atlas_AAL==74)+(atlas_AAL==75)+(atlas_AAL==76);
slice_subcortical = 35;
figure;
imshow(img_fmri(:,:,slice_subcortical,60), []);
colormap(gray);
% Create a colormap with transparency and red
cmap = [0 0 0; 0.8500 0.3250 0.0980]; % Transparent for 0, Red for 1

% Overlay colored mask with transparency
hold on;
h = imshow(ind2rgb(im_subcortical(:,:,slice_subcortical) + 1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
set(h, 'AlphaData', im_subcortical(:,:,slice_subcortical)); % Set transparency based on mask values

t1 = readtable('/home/iesteves/FC/files/atlases/Schaefer/Schaefer2018_100Parcels_7Networks_order_new.txt');
t2 = readtable('/home/iesteves/FC/files/atlases/Schaefer/coord_SchaeferSubCRB7100_130.txt');
ind = t1.Var1;
net=t2.Var4(1:100);
m_ind_net = [ind,net];

%%
% FPN, DMN, DAN, LN, VAN, SMN, VN, SUB, CRB
networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN'};
slice_cerebellar = 12;
slice_subcortical = 35;
slices = [61, 55, 62, 26, 51, 62, 30];

c_SUB = [0.8500 0.3250 0.0980];
c_CRB = [0.25 0.25 0.25];

colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330];

figure('pos', [50 20 1000 900]);
for net_ind = 1:7  
    im_net = ismember(atlas_Schaefer, m_ind_net(m_ind_net(:,2)==net_ind,1));
    slice = slices(net_ind);
     
    subplot(2, 5, net_ind)
    imshow(imrotate(img_fmri(:,:,slice,60), 90, 'bilinear', 'crop'), []);
    colormap(gray);
    cmap = [0 0 0; colors(net_ind,:)]; % Transparent for 0, Red for 1
    hold on;
    h = imshow(ind2rgb(imrotate(im_net(:,:,slice),90, 'bilinear', 'crop') + 1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
    set(h, 'AlphaData', imrotate(im_net(:,:,slice),90, 'bilinear', 'crop')); % Set transparency based on mask values
    title(networks{net_ind})
    hold off;
    set(gca, 'FontSize',14)
end
subplot(2, 5, 8)
subcortical_regions = [37, 38, 41, 42, 71:78];
im_subcortical = ismember(atlas_AAL, subcortical_regions);
imshow(imrotate(img_fmri(:,:,slice_subcortical,60), 90, 'bilinear', 'crop'), []);
colormap(gray);
cmap = [0 0 0; c_SUB]; % Transparent for 0, Red for 1
hold on;
h = imshow(ind2rgb(imrotate(im_subcortical(:,:,slice_subcortical), 90, 'bilinear', 'crop') + 1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
set(h, 'AlphaData', imrotate(im_subcortical(:,:,slice_subcortical), 90, 'bilinear', 'crop')); % Set transparency based on mask values
title('SUB')
set(gca, 'FontSize',14)

subplot(2, 5, 9)
cerebellar_regions = 91:116;
im_cerebellar = ismember(atlas_AAL, cerebellar_regions);
imshow(imrotate(img_fmri(:,:,slice_cerebellar,60), 90, 'bilinear', 'crop'), []);
colormap(gray);
cmap = [0 0 0; c_CRB]; % Transparent for 0, Red for 1
hold on;
h = imshow(ind2rgb(imrotate(im_cerebellar(:,:,slice_cerebellar), 90, 'bilinear', 'crop') + 1, cmap), 'InitialMagnification', 'fit'); % Adding 1 to make 0 values transparent
set(h, 'AlphaData', imrotate(im_cerebellar(:,:,slice_cerebellar), 90, 'bilinear', 'crop')); % Set transparency based on mask values
title('CRB')
set(gca, 'FontSize',14)
%print('/home/iesteves/FC/figures/atlases/networks', '-dpng')


