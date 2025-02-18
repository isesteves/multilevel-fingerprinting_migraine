%%% STILL NEED TO CONVERT TO FUNCTION AND DEFINE OUTPUTS AND INPUTS, USE AS
%%% SCRIPT FOR NOW

% tasks = {'task-breathhold', 'task-deepbreathing_run-02'};
% sessions = {'ses-midcycle'};
% subjects = {'sub-control019', 'sub-control020', 'sub-control025', 'sub-control026', 'sub-control029', 'sub-control030', 'sub-control031'};

addpath('/home/iesteves/MIG_N2Treat_fMRI/NIfTI_20140122')
path = '/home/iesteves/eeg-fmri';
figurepath = '/home/iesteves/MIG_N2Treat_fMRI/Figures/'; 

dataset = 'MIGN2TREAT';
subjects = {'sub-control019', 'sub-control020', 'sub-control025','sub-control026', 'sub-control027', 'sub-control029', 'sub-control030', 'sub-control031',...
    'sub-control033', 'sub-control044', 'sub-control046'};
%subjects = {'sub-control019', 'sub-control020', 'sub-control026', 'sub-control027', 'sub-control028','sub-control029', 'sub-control030', 'sub-control031',...
%    'sub-control033', 'sub-control046'};
%subject = 'sub-control025';
task = 'task-rest';
session = 'ses-midcycle';
mask_types = {'GM', 'WM', 'CSF'};
nr_vol = 333; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE

matlab_vers  = '2016';

exist_dvars = 1;
exist_refrms = 1;

%pipelines = {'', '_preprocessed_','_preprocessed_ica_mo_reg', '_preprocessed_ica_mo_csf_reg', '_preprocessed_ica_mo_csf_wm_reg', '_preprocessed_ica_mo_csf_wm_gs_reg'};
pipelines = {'icafixrp', 'icafixrp_mo_nui',  'icafixrp_mo_csf_wm_nui',  'ica_rp_mo_csf_wm_nui', ...
   'icafix', 'ica_rp_nui', 'ica_rp_mo_nui', 'rp_mo_csf_wm_nui'};

%% STILL NEED TO USE GENERIC OUTSIDE PATH
masks_pipelines_gm = zeros(length(subjects), nr_vol, length(pipelines));
masks_pipelines_wm = zeros(length(subjects), nr_vol, length(pipelines));
masks_pipelines_csf = zeros(length(subjects), nr_vol, length(pipelines));
for s = 1:length(subjects)
    subject = subjects{s};
    for p = 1:length(pipelines)

        funcpipeline = pipelines{p}; 

        funcpath = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/filtered_func_data_preprocessed_', funcpipeline,'.nii.gz'];
        brainmaskpath = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mask.nii.gz'];

        if strcmp(matlab_vers, '2016')
            %functional image
            gunzip(funcpath)
            file = load_untouch_nii(funcpath(1:end-3));
            delete([funcpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;

            %brain mask
            gunzip(brainmaskpath)
            file_brainmask = load_untouch_nii(brainmaskpath(1:end-3));
            delete([brainmaskpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;

        else
            % functional image
            file = load_nii(funcpath(1:end-3)); %%%%% CONFIRM

            % brain mask
            file_brainmask = load_nii(brainmaskpath(1:end-3)); %%%%% CONFIRM
        end

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

            if length(unique(mask_img(:)))>2
                thr = 0.8; % in case the mask was not binarized; %%%%%%% CHECK THIS

                mask(mask_img >= thr) = 1;
                mask(mask_img < thr) = 0;
            end

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

        gm_mean_vol = mean(gm_data);
        gm_psc = ((gm_mean_vol-mean(gm_mean_vol))/mean(gm_mean_vol))*100;

        wm_mean_vol = mean(wm_data);
        wm_psc = ((wm_mean_vol-mean(wm_mean_vol))/mean(wm_mean_vol))*100;

        csf_mean_vol = mean(csf_data);
        csf_psc = ((csf_mean_vol-mean(csf_mean_vol))/mean(csf_mean_vol))*100;

        masks_pipelines_gm(s, :, p) = gm_psc;
        masks_pipelines_wm(s, :, p) = wm_psc;
        masks_pipelines_csf(s, :, p) = csf_psc;

        %% Load Motion outliers
        mofile = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mo_confound_dvars.txt'];
        mo_data = load(mofile, '-ascii');

        %% Load MO metric
        dvarsfile = [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mc/dvars.txt'];
        dvarsdata = load(dvarsfile, '-ascii');

        mo_metric = 'dvars';
        mo_metric_data = dvarsdata;

        %% Plot carpet

        figure('pos', [20 20 1000 900]);
        h = suptitle({[subject,' - ', session, ' - ', task, ' - ', funcpipeline], '  ', ' '});
        h.Interpreter = 'none';

        plot_carpet(gm_data, wm_data, csf_data, mo_metric, mo_metric_data, mo_data)

        print([figurepath, 'carpet_',subject, '_', session, '_', task, funcpipeline], '-dpng')

    end
end

%% 
save([session, '_masks_pipelines_gm.mat'], 'masks_pipelines_gm');
save([session, '_masks_pipelines_wm.mat'], 'masks_pipelines_wm');
save([session, '_masks_pipelines_csf.mat'], 'masks_pipelines_csf');

%ind_pipelines = [1 ,2, 3, 5,8];
%ind_pipelines = [4, 6, 7];
ind_pipelines = 1:8;
for s = 1:length(subjects)
    subject = subjects{s};
    
    figure('pos', [50 50 1200 600]);
    h = suptitle({[subject,' - ', session, ' - ', task], '  ', ' '});  
    h.Interpreter = 'none';

    for p = ind_pipelines

        funcpipeline = pipelines{p}; 
        subplot(3,1,1)
        plot(masks_pipelines_gm(s,:,p))
        hold on
        title('GM')
        ylabel('PSC')
        xlim([1 nr_vol])
        %ylim([-2 2])

        subplot(3,1,2)
        plot(masks_pipelines_wm(s,:,p))
        hold on
        title('WM')
        ylabel('PSC')
        xlim([1 nr_vol])
        %ylim([-2 2])

        subplot(3,1,3)
        plot(masks_pipelines_csf(s,:,p))
        hold on
        title('CSF')
        ylabel('PSC')
        xlim([1 nr_vol])
        xlabel('Volume Number')
        %ylim([-2 2])
 
    end
    l = legend(pipelines(ind_pipelines));
    l.Interpreter = 'none';
    l.Position = [0.5, 0.01, 0.03, 0.03];
    l.Orientation = 'Horizontal';
    
    print([figurepath, 'pipeplines_avgpsc_',subject, '_', session, '_', task], '-dpng')
end

% %% Desikan parcellation
% gunzip('/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/task-rest/sub-control019/ses-midcycle/minusy/filtered_func_data_preprocessed2standard.nii.gz');
% stdimg = load_untouch_nii('/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/task-rest/sub-control019/ses-midcycle/minusy/filtered_func_data_preprocessed2standard.nii');
% delete('/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/task-rest/sub-control019/ses-midcycle/minusy/filtered_func_data_preprocessed2standard.nii')
% 
% func_img = stdimg.img;
% n_x = size(func_img, 1); 
% n_y = size(func_img, 2); 
% n_z = size(func_img, 3); 
% nr_vol = size(func_img, 4);
% 
% dskfile = load_untouch_nii('/home/iesteves/MIG_N2Treat_fMRI/Desikan_MNI152-2mm.nii');
% idx = setdiff(1:70, [1 2 36]);
% parceldata = cell(1,length(idx));
% for p = 1:length(idx)          
%    
%             % parcel mask: #voxels x time (TRs)
%             mask = dskfile.img==idx(p);
%             mask_tc = reshape(mask, [n_x*n_y*n_z 1]);
% 
%             % functional image: #voxels x time (TRs)
%             func_tc = reshape(func_img, [n_x*n_y*n_z nr_vol]);
% 
%             % functional image with GM/WM/CSF mask
%             mfunc_tc = func_tc(logical(mask_tc),:);
% 
%             parceldata{1,p} = mfunc_tc;  
% end
%%
% masks_data = cell(1,3); 
% for m = 1:length(mask_types)
%     
%     mask_type = mask_types{m};
%     maskpath =  [path, '/', dataset, '/PREPROCESS/', task,'/', subject, '/', session,'/minusy/masks/EF_',mask_type,'_ero.nii.gz'];
% 
%     gunzip(maskpath)
%     file_mask = load_untouch_nii(maskpath(1:end-3));
%     delete([maskpath(1:end-3)]) % no need to have unzipped file stored + creates conflict when using FEAT if both .nii and .nii.gz exist;
% 
%     %%
%     func_img = file.img;
%     n_x = size(func_img, 1); 
%     n_y = size(func_img, 2); 
%     n_z = size(func_img, 3); 
%     nr_vol = size(func_img, 4);
%     
%     brainmask = file_brainmask.img;
%     
%     vol_brainmask = repmat(brainmask, [1 1 1 nr_vol]);
% 
%     bmfunc_img = func_img;
%     bmfunc_img(~vol_brainmask) = 0; 
%     
%     mask_img = file_mask.img;
% 
%     mask_val = mask_img(:);
%     thr = 0.8; % in case the mask was not binarized; %%%%%%% CHECK THIS
% 
%     mask = mask_img;
%     mask(mask_img >= thr) = 1;
%     mask(mask_img < thr) = 0;
%     
% 
%     vol_mask = repmat(mask, [1 1 1 nr_vol]);
% 
%     mfunc_img = func_img;
%     mfunc_img(~vol_mask) = 0; 
% 
%     v = floor(nr_vol/2);
%     slices = [20, 30, 40];
%     c = 1;
%     fig = figure('pos', [20 20 1100 600]);
%     h = suptitle({[subject,' - ', session, ' - ', task, ' -  vol = ', num2str(v)], '  '});
%     h.Interpreter = 'none';
%     for sk = 1:length(slices)
%         s = slices(sk);
% 
%         subplot(3,3,c)
%         imshow(mask(:,:,s), [], 'InitialMagnification', 800)
%         title([mask_type, ' mask'])
%         ylabel(['slice = ', num2str(s), ])
% 
%         subplot(3,3,c+1)
%         imshow(func_img(:,:,s,v), [], 'InitialMagnification', 800)
%         title('Func')
% 
%         subplot(3,3,c+2)
%         imshow(mfunc_img(:,:,s, v), [], 'InitialMagnification', 800)
%         title(['Func with ', mask_type,' mask'])
% 
%         c = c + 3;
% 
%     end
%     
%    
%     % GM/WM/CSF mask: #voxels x time (TRs)
%     mask_tc = reshape(mask, [n_x*n_y*n_z 1]);
%     
%     % functional image: #voxels x time (TRs)
%     func_tc = reshape(func_img, [n_x*n_y*n_z nr_vol]);
%     
%     % functional image with GM/WM/CSF mask
%     mfunc_tc = func_tc(logical(mask_tc),:);
% 
%     masks_data{1,m} = mfunc_tc;
%     
% end
% 
% %% Carpet plot
% 
% gm_data = masks_data{1,1};
% wm_data = masks_data{1,2};
% csf_data = masks_data{1,3};
% 
% figure;
% imagesc([gm_data; wm_data; csf_data]);
% set(gca,'ytick',[]);
% colormap(gray);
% hgm=hline(size(gm_data,1),'g');
% hwm=hline(size([gm_data; wm_data],1),'b');
% hcsf=hline(size([gm_data; wm_data; csf_data],1),'r');
% xlim([-3 nr_vol]);
% title('Carpet plot');
% ylabel('Voxel Intensity');
% xlabel('Volume Number');
% 
% x1=-3; x2=0;
% y1=0; y2=size(gm_data,1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g1 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.9 0.3]);
% 
% x1=-3; x2=0;
% y1=size(gm_data,1); y2=size([gm_data; wm_data],1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g2 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.1 1]);
% 
% x1=-3; x2=0;
% y1=size([gm_data; wm_data],1); y2=size([gm_data; wm_data; csf_data],1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g3 = patch('Faces',f,'Vertices',v,'FaceColor',[1 0.1 0.1]);
% 
% g = [g1(1), g2(1), g3(1)];
% legend(g, 'GM', 'WM', 'CSF', 'Orientation', 'Horizontal', 'Location', 'southoutside')
% 
% 
% %% Load Motion parameters
% mcfile = ['/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/', task,'/', subject, '/', session,'/minusy/prefiltered_func_data_mcf.txt'];
% mcdata = load(mcfile, '-ascii');
% 
% figure;
% subplot(2,1,1)
% plot(mcdata(:,1:3))
% legend('rx', 'ry', 'rz')
% title('MCFLIRT estimated rotations')
% ylabel('mm')
% xlabel('Volume number')
% xlim([1 nr_vol])
% 
% subplot(2,1,2)
% plot(mcdata(:,4:6))
% legend('tx', 'ty', 'tz')
% title('MCFLIRT estimated translations')
% ylabel('mm')
% xlabel('Volume number')
% xlim([1 nr_vol])
% 
% %% Load Motion outliers
% 
% mofile = ['/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mo_confound.txt'];
% modata = load(mofile, '-ascii');
% 
% [irow,jcol,v] = find(modata);
% 
% %% Load FSL motion outliers metrics for comparison if available
% 
% if exist_refrms
%     refrmsfile = ['/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mo_refrms.txt'];
%     refrmsdata = load(refrmsfile, '-ascii');
% else
%     refrmsdata = [];
% end
% 
% if exist_dvars
%     dvarsfile = ['/home/iesteves/MIG_N2Treat_fMRI/PREPROCESS/', task,'/', subject, '/', session,'/minusy/mo_dvars.txt'];
%     dvarsdata = load(dvarsfile, '-ascii');
% else
%     dvarsdata = [];
% end
% 
% %% MOTION OUTLIER MEASURES
% % masks: It is possible to specify a mask to restrict the spatial region over
% % which the calculations are done (for the top three metrics). This, when used, 
% % is usually a brain mask to avoid calculations based on non-brain regions. 
% % It is not necessary to use and results are usually very good without needing any masking.
% 
% 
% % Brain mask for refrms and DVARS
% 
% perc_mask= sum(brainmask(:)==1)/numel(brainmask); % ==mean(brainmask(:));
% 
% % brain mask: #voxels x time (TRs)
% brainmask_tc = reshape(brainmask, [n_x*n_y*n_z 1]);
% 
% % functional image with brain mask: #voxels x time (TRs)
% bmfunc_tc = func_tc(logical(mask_tc),:);
% 
% %% Framewise displacement (FD)
% dermcdata = mcdata-circshift(mcdata, 1);
% dermcdata(1,:) = 0;
% fd = sum(abs(dermcdata),2);
% 
% figure; 
% plot(fd)
% xlabel('Volume number')
% title('Framewise Displacement')
% 
% 
% %% refrms: RMS intensity difference of volume N to the reference volume
% ref_vol_ind = floor(nr_vol/2);
% ref_vol = bmfunc_tc(:,ref_vol_ind);
% 
% ref_vol_rep = repmat(ref_vol(:), [1 nr_vol]);
% diff_vol_matrix = (bmfunc_tc - ref_vol_rep)/median(bmfunc_tc(:));
% refrms = sqrt(mean(diff_vol_matrix.^2, 1));
% 
% figure; plot(refrms)
% hold on;
% plot(refrmsdata)
% xlabel('Volume number')
% title('refrms')
% legend('refrms', 'refrms FSL')
% 
% % TRY TO MATCH FSL VERSION - NOT TESTED FULLY 
% diff_filter = refrms - circshift(refrms,1);
% 
% %% DVARS: RMS intensity difference of volume N to volume N+1 (see Power et al, NeuroImage, 59(3), 2012)
% 
% % PREVIOUS VERSION - CHANGED TO MATCH FSL SCALING
% % diff_vol_matrix = 1000*((bmfunc_tc -circshift(bmfunc_tc,1, 2))/median(bmfunc_tc(:)));
% % dvars = sqrt(mean(diff_vol_matrix.^2, 1));
% % dvars(1,1) = 0;
% 
% diff_vol_matrix = bmfunc_tc -circshift(bmfunc_tc,1, 2);
% dvars_unscaled = sqrt(mean(diff_vol_matrix.^2, 1)/perc_mask);
% dvars = 1000*(sqrt(mean(diff_vol_matrix.^2, 1)/perc_mask))/median(bmfunc_tc(:));
% dvars(1,1) = 0;
% 
% thr = prctile(dvars, 75)+1.75*iqr(dvars); % THRESHOLD THAT FSL USES
% 
% figure; 
% yrange = [min(dvars)-0.5*std(dvars) max(dvars)+0.5*std(dvars)];
% plot(dvars, '-*')
% hold on;
% plot(dvarsdata)
% line([1 nr_vol], [thr thr])
% line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r')
% xlabel('Volume number')
% title('DVARS')
% xlim([1 nr_vol])
% ylim(yrange)
% legend('dvars', 'dvars FSL')
% 
% %% Motion parameters + motion outliers metrics 
% 
% figure('pos', [20 20 900 900]);
% h = suptitle({[subject,' - ', session, ' - ', task], '  ', ' '})
% h.Interpreter = 'none';
% subplot(5,1,1)
% plot(mcdata(:,1:3))
% legend('rx', 'ry', 'rz')
% title('MCFLIRT estimated rotations')
% ylabel('mm')
% xlim([1 nr_vol])
% %xlabel('Volume number')
% 
% subplot(5,1,2)
% plot(mcdata(:,4:6))
% legend('tx', 'ty', 'tz')
% title('MCFLIRT estimated translations')
% ylabel('mm')
% xlim([1 nr_vol])
% %xlabel('Volume number')
% 
% % Framewise displacement (FD)
% yrange = [0 max(fd)+0.2*max(fd)];
% subplot(5,1,3)
% plot(fd, '-*')
% hold on
% line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r')
% title('Framewise Displacement (FD)')
% xlim([1 nr_vol])
% ylim(yrange)
% %xlabel('Volume number')
% 
% % refrms 
% yrange = [min(refrms)-0.5*std(refrms) max(refrms)+0.5*std(refrms)];
% subplot(5,1,4)
% plot(refrms, '-*')
% hold on
% line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r')
% title('refrms')
% xlim([1 nr_vol])
% ylim(yrange)
% %xlabel('Volume number')
% 
% % dvars
% subplot(5,1,5)
% yrange = [min(dvars)-0.5*std(dvars) max(dvars)+0.5*std(dvars)];
% plot(dvars, '-*')
% hold on
% go = line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r');
% title('dvars')
% xlim([1 nr_vol])
% ylim(yrange)
% xlabel('Volume number')
% 
% hl= legend(go(1), 'FSL Motion outliers');
% hl.Position = [0.5 0.03 0.01 0.01];
% 
% %% Carpet plot with dvars
% 
% figure('pos', [20 20 1000 900]);
% h = suptitle({[subject,' - ', session, ' - ', task], '  ', ' '});
% h.Interpreter = 'none';
% 
% % dvars 
% subplot(4,1,1)
% yrange = [min(dvars)-0.5*std(dvars) max(dvars)+0.5*std(dvars)];
% plot(dvars, '-*')
% hold on
% go = line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r');
% title('dvars')
% xlim([1 nr_vol])
% ylim(yrange)
% xlabel('Volume number')
% 
% subplot(4,1,[2 3 4])
% imagesc([gm_data; wm_data; csf_data]);
% set(gca,'ytick',[]);
% colormap(gray);
% hgm=hline(size(gm_data,1),'g');
% hwm=hline(size([gm_data; wm_data],1),'b');
% hcsf=hline(size([gm_data; wm_data; csf_data],1),'r');
% xlim([-3 nr_vol]);
% title('Carpet plot');
% ylabel('Voxel Intensity');
% xlabel('Volume Number');
% 
% x1=-3; x2=0;
% y1=0; y2=size(gm_data,1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g1 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.9 0.3]);
% 
% x1=-3; x2=0;
% y1=size(gm_data,1); y2=size([gm_data; wm_data],1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g2 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.1 1]);
% 
% x1=-3; x2=0;
% y1=size([gm_data; wm_data],1); y2=size([gm_data; wm_data; csf_data],1);
% v=[x1 y1; x1 y2; x2 y2; x2 y1];
% f=[1 2 3 4];
% g3 = patch('Faces',f,'Vertices',v,'FaceColor',[1 0.1 0.1]);
% 
% g = [g1(1), g2(1), g3(1)];
% legend(g, 'GM', 'WM', 'CSF', 'Orientation', 'Horizontal', 'Location', 'southoutside')
% 
% 
% %%
% 
