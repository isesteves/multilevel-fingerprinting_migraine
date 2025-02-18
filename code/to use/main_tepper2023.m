% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;


%% Networks matrix

matrix_all = zeros(128, 128, 9);
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 128 9];
for g = 1:size(areas,1)
    matrix_all(areas(g,1):areas(g,2),areas(g,1):areas(g,2),g) = g;
    matrix_all(:,:,g) = matrix_all(:,:,g) - diag(diag(matrix_all(:,:,g)));
end

matrix_FC_all = sum(matrix_all,3);

figure; 
imagesc(matrix_FC_all);
hold on;
plot_atlas_labels('SchaeferSubCRB7100', 128)

% get the vectorized upper triangular matrix of the networks template
tril_m  = tril(true(size(matrix_FC_all)), -1);
v  = matrix_FC_all(tril_m).';

matrix_FC_all(~tril_m)=NaN;

figure;
imagesc(matrix_FC_all)

figure;
imagesc(matrix_FC_all(tril_m))

figure;
imagesc(v)

%%
% load Pearson correlation variables (iFC_all - areas x areas x subject_session, ordered by session type)
pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';

%load(['/home/iesteves/dFC/sFC/MIGcomplete_nonzero-parcels-50/Pearson-results-patient_control-', atlas,'-', pipeline,'.mat'])
load(['/home/iesteves/dFC/sFC/MIGcomplete_nonzero-parcels-50/Pearson-results-patient_control-', atlas,'-', pipeline,'.mat'])
nr_areas = size(iFC_all, 1);
nr_subjects = size(iFC_all,3);


% get the vectorized lower triangular matrix
group_data_persession = zeros((N_areas*N_areas-N_areas)/2, n_Subjects);
for k = 1:n_Subjects
    %group_data_ltm = tril(iFC_all(:,:,k), -1)';
    subject_FC =  iFC_all(:,:,k);
    tril_m  = tril(true(size(subject_FC)), -1);
    group_data_vec  = subject_FC(tril_m).';
    %group_data_vec = group_data_ltm(group_data_ltm~=0);
    group_data_persession(:,k) = group_data_vec;
end

% reorganize correlation matrices by subject 
order = [1:10:31 2:10:32 3:10:33 4:10:34 5:10:35 6:10:36 7:10:37 8:10:38 9:10:39 10:10:40 55 41 56 42 57  43 58 44  59  45  60 46 61 47 62  48 63  49 64  50 65  51 66  52 67  53 68  54];
filenames_persubject = filenames(order);
data = group_data_persession(:, order);

group_data = (data);

% FIGURE -
figure('pos', [200 1000 800 1200]);
imagesc(group_data);
ylabel('FC edges')
xlabel('Data samples (All sessions for all subjects)')
    
%% Generate IA and IB template matrices 
within_session_aux = zeros(nr_subjects, nr_subjects);
within_subject_aux = zeros(nr_subjects, nr_subjects);
within_group_aux = zeros(nr_subjects, nr_subjects);
within_menstrual_session_aux = zeros(nr_subjects, nr_subjects);
between_group_within_session_aux = zeros(nr_subjects, nr_subjects);

for j = 1:nr_subjects
    
    for t = 1:nr_subjects
        
        % get file 1 and split to get subject and session
        filename1 = filenames_persubject{j};
        filename1_aux = strsplit(filename1, '/');
        subject1 = filename1_aux{6};
        session1 = filename1_aux{7};
        
        % get file 2 and split to get subject and session
        filename2 = filenames_persubject{t};
        filename2_aux = strsplit(filename2, '/');
        subject2 = filename2_aux{6};
        session2 = filename2_aux{7};
        
        % same session
        if strcmp(session1, session2)
            within_session_aux(j,t) = 1;
        end
        
        % same subject
        if strcmp(subject1, subject2)
            within_subject_aux(j,t) = 1;
        end
        
        % same group
        if strcmp(subject1(1:end-3), subject2(1:end-3))
            within_group_aux(j,t) = 1;
        end
        
        % different groups, same menstrual cycle phase
        cond1= any(strcmp(session1, {'ses-preictal', 'ses-ictal', 'ses-postictal'})) & any(strcmp(session2, {'ses-premenstrual'}));
        cond2= any(strcmp(session2, {'ses-preictal', 'ses-ictal', 'ses-postictal'})) & any(strcmp(session1, {'ses-premenstrual'}));
        cond3= any(strcmp(session1, {'ses-interictal'})) & any(strcmp(session2, {'ses-midcycle'}));
        cond4= any(strcmp(session1, {'ses-midcycle'})) & any(strcmp(session2, {'ses-interictal'}));
        if cond1 || cond2 || cond3 || cond4
            between_group_within_session_aux(j,t) = 1;
        end
        
        % same menstrual cycle phase
        cond1= any(strcmp(session1, {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-premenstrual'})) & any(strcmp(session2,  {'ses-preictal', 'ses-ictal', 'ses-postictal', 'ses-premenstrual'}));
        cond2= any(strcmp(session1, {'ses-midcycle', 'ses-interictal'})) & any(strcmp(session2, {'ses-midcycle', 'ses-interictal'}));
        if cond1 || cond2 || cond3 || cond4
            within_menstrual_session_aux(j,t) = 1;
        end
        
    end
    
end


figure('pos', [200 1000 1800 300]);
subplot(1,5,1)
imagesc(within_group_aux);
title('Within group')

subplot(1,5,2)
imagesc(within_subject_aux);
title('Within subject')

subplot(1,5,3)
imagesc(within_session_aux);
title('Within session')

subplot(1,5,4)
imagesc(within_menstrual_session_aux);
title('W menstrual ses')

subplot(1,5,5)
imagesc(between_group_within_session_aux);
title('B group w menstrual ses')

%% code diagonal and upper triangular part of the template matrices as 0.5 to ignore these values
triu_ind = ones(nr_subjects, nr_subjects) ;
triu_ind = triu(triu_ind) ; 

within_group = within_group_aux;
within_group(triu_ind==1) = 0.5 ;

within_subject = within_subject_aux;
within_subject(triu_ind==1) = 0.5 ; 

within_session = within_session_aux;
within_session(triu_ind==1) = 0.5 ; 

within_menstrual_session = within_menstrual_session_aux;
within_menstrual_session(triu_ind==1) = 0.5 ; 

between_group_within_session = between_group_within_session_aux;
between_group_within_session(triu_ind==1) = 0.5 ; 

figure('pos', [200 1000 1800 300]);
subplot(1,5,1)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_group==0.5)=0;
imagesc(within_group,'AlphaData',imAlpha);
title('Within group')

subplot(1,5,2)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_subject==0.5)=0;
imagesc(within_subject,'AlphaData',imAlpha);
title('Within subject')

subplot(1,5,3)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_session==0.5)=0;
imagesc(within_session,'AlphaData',imAlpha);
title('Within session')

subplot(1,5,4)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_menstrual_session==0.5)=0;
imagesc(within_menstrual_session,'AlphaData',imAlpha);
title('W menstrual ses')

subplot(1,5,5)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(between_group_within_session==0.5)=0;
imagesc(between_group_within_session,'AlphaData',imAlpha);
title('B group w menstrual ses')

%% Generate IA and IB template matrices - more specific
wgroup_corr_aux = zeros(nr_subjects, nr_subjects);
wgroup_bsession_corr_aux =  zeros(nr_subjects, nr_subjects);
wsession_corr_aux =  zeros(nr_subjects, nr_subjects);
bgroup_corr_aux = zeros(nr_subjects, nr_subjects);

for j = 1:nr_subjects
    
    for t = 1:nr_subjects
        
        % get file 1 and split to get subject and session
        filename1 = filenames_persubject{j};
        filename1_aux = strsplit(filename1, '/');
        subject1 = filename1_aux{6};
        session1 = filename1_aux{7};
        
        % get file 2 and split to get subject and session
        filename2 = filenames_persubject{t};
        filename2_aux = strsplit(filename2, '/');
        subject2 = filename2_aux{6};
        session2 = filename2_aux{7};
        
       
        if strcmp(session1,'ses-preictal')
            if strcmp(session2, 'ses-preictal')
                wsession_corr_aux(j, t) = 1;
            elseif strcmp(session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) =  1;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) =  2;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  3;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  1;
            end
            
        elseif strcmp(session1,'ses-ictal')
            if strcmp(session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 1;
            elseif strcmp(session2, 'ses-ictal')
                wsession_corr_aux(j, t) = 2;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) =  4;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  5;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  2;
            end
        
        elseif strcmp(session1,'ses-postictal')
            if strcmp (session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 2;
            elseif strcmp (session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) = 4;
            elseif strcmp(session2, 'ses-postictal')
                wsession_corr_aux(j, t) = 3;
            elseif strcmp(session2, 'ses-interictal')
                wgroup_bsession_corr_aux(j, t) =  6;
            elseif strcmp(session2, 'ses-premenstrual')
                bgroup_corr_aux(j, t) =  3;
            end
        
        elseif strcmp(session1,'ses-interictal') 
            if strcmp (session2, 'ses-preictal')
                wgroup_bsession_corr_aux(j, t) = 3;
            elseif strcmp(session2, 'ses-ictal')
                wgroup_bsession_corr_aux(j, t) = 5;
            elseif strcmp(session2, 'ses-postictal')
                wgroup_bsession_corr_aux(j, t) = 6;
            elseif strcmp(session2, 'ses-interictal')
                wsession_corr_aux(j, t) = 4;
            elseif strcmp(session2, 'ses-midcycle')
                bgroup_corr_aux(j, t) =  4;
            end
        
        elseif strcmp(session1,'ses-premenstrual') 
            if strcmp(session2, 'ses-premenstrual')
                wsession_corr_aux(j, t) = 5;
            elseif strcmp(session2, 'ses-midcycle')
                wgroup_bsession_corr_aux(j, t) =  7;
            elseif strcmp(session2, 'ses-preictal')
                bgroup_corr_aux(j, t) =  1;
            elseif strcmp(session2, 'ses-ictal')
                bgroup_corr_aux(j, t) =  2;
            elseif strcmp(session2, 'ses-postictal')
                bgroup_corr_aux(j, t) =  3;
            end
        
        elseif strcmp(session1,'ses-midcycle')
            if strcmp(session2, 'ses-premenstrual')
                wgroup_bsession_corr_aux(j, t) =  7;
            elseif strcmp(session2, 'ses-midcycle')
                wsession_corr_aux(j, t) = 6;
            elseif strcmp(session2, 'ses-interictal')
                bgroup_corr_aux(j, t) =  4;
            end
        end

        % same group
        if strcmp(subject1(1:end-3), subject2(1:end-3)) & strcmp(subject1(1:end-3), 'sub-patient')
            wgroup_corr_aux(j,t) = 1;
        elseif  strcmp(subject1(1:end-3), subject2(1:end-3)) & strcmp(subject1(1:end-3), 'sub-control')
            wgroup_corr_aux(j, t) = 2;
        end
        
        
    end
    
end

%% code diagonal and upper triangular part of the template matrices as 0.5
triu_ind = ones(nr_subjects, nr_subjects);
triu_ind = triu(triu_ind) ; 

wgroup_corr = wgroup_corr_aux;
wgroup_corr(triu_ind==1)= 0.5;

wgroup_bsession_corr =wgroup_bsession_corr_aux;
wgroup_bsession_corr(triu_ind==1) = 0.5 ;

wsubject_bsession_corr_aux = within_subject.*wgroup_bsession_corr;
wsubject_bsession_corr = wsubject_bsession_corr_aux;
wsubject_bsession_corr(triu_ind==1) = 0.5 ;

wsession_corr = wsession_corr_aux;
wsession_corr(triu_ind==1) = 0.5;

bgroup_corr = bgroup_corr_aux;
bgroup_corr(triu_ind==1) = 0.5;

%% PCA
nr_networks = 9;
corr_networks = zeros(size(group_data, 2), size(group_data, 2), nr_networks+1, length(components));
networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'Subcortical', 'Cerebellum', 'WB'};

components = [1, 8, 12, 35, 40, 68];
mean_data = mean(data);
PCA_group_data = data-mean_data;

% FIGURE - corr matrix for different nb of PCs and networks
for comp = 1:length(components)
    nComp = components(comp);
    [coeffs, score, ~] = pca(PCA_group_data, 'NumComponents', nComp);
    PCA_matrix = score * coeffs';
    group_PCA_recon = bsxfun(@plus, PCA_matrix, mean_data);

    figure('pos', [100, 200, 1900 400]);
    for n = 1:nr_networks+1
        if n <= nr_networks
            network_data_PCA = group_PCA_recon(v==n,:); 
        else
            network_data_PCA = group_PCA_recon;
        end

        corr_val = corr(network_data, network_data, 'Type', 'Pearson');
        corr_networks(:,:,n,comp) = corr_val;

        subplot(length(components), nr_networks+1, n)
        imagesc(corr_val)
        title(networks{n})
        if n == 1
            ylabel([num2str(nComp), ' PCs'])
        end

    end   
end
%%

corr_matrix = corr_val;

bdata_wgroup= NaN*ones(size(corr_matrix(:),1),2, nr_subjects);
bdata_wgroup_bsession= NaN*ones(size(corr_matrix(:),1),7, nr_subjects);
bdata_wsubject_bsession= NaN*ones(size(corr_matrix(:),1),7, nr_subjects);
bdata_wsession= NaN*ones(size(corr_matrix(:),1),6, nr_subjects);
bdata_bgroup= NaN*ones(size(corr_matrix(:),1),4, nr_subjects);
bdata_wgroup_PCA= NaN*ones(size(corr_matrix(:),1),2, nr_subjects);
bdata_wgroup_bsession_PCA= NaN*ones(size(corr_matrix(:),1),7, nr_subjects);
bdata_wsubject_bsession_PCA= NaN*ones(size(corr_matrix(:),1),7, nr_subjects);
bdata_wsession_PCA= NaN*ones(size(corr_matrix(:),1),6, nr_subjects);
bdata_bgroup_PCA= NaN*ones(size(corr_matrix(:),1),4, nr_subjects);
for n = 1:nr_networks+1
 
    corr_matrix = corr_networks(:,:,n);
    corr_matrix_PCA = corr_networks_PCA(:,:,n);

    bdata_wgroup(1:sum(sum(wgroup_corr == 1)),1, n) = corr_matrix(wgroup_corr == 1);
    bdata_wgroup(1:sum(sum(wgroup_corr == 2)),2, n) = corr_matrix(wgroup_corr == 2);
    
    bdata_wgroup_PCA(1:sum(sum(wgroup_corr == 1)),1, n) = corr_matrix_PCA(wgroup_corr == 1);
    bdata_wgroup_PCA(1:sum(sum(wgroup_corr == 2)),2, n) = corr_matrix_PCA(wgroup_corr == 2);


    for k = 1:7
        bdata_wgroup_bsession(1:sum(sum(wgroup_bsession_corr == k)),k, n) = corr_matrix(wgroup_bsession_corr == k);
        bdata_wsubject_bsession(1:sum(sum(wsubject_bsession_corr == k)),k, n) = corr_matrix(wsubject_bsession_corr == k);
        
        bdata_wgroup_bsession_PCA(1:sum(sum(wgroup_bsession_corr == k)),k, n) = corr_matrix_PCA(wgroup_bsession_corr == k);
        bdata_wsubject_bsession_PCA(1:sum(sum(wsubject_bsession_corr == k)),k, n) = corr_matrix_PCA(wsubject_bsession_corr == k);
    end
    
    for k = 1:6
        bdata_wsession(1:sum(sum(wsession_corr == k)),k, n) = corr_matrix(wsession_corr == k);
        
        bdata_wsession_PCA(1:sum(sum(wsession_corr == k)),k, n) = corr_matrix_PCA(wsession_corr == k);
    end
    
    for k = 1:4
        bdata_bgroup(1:sum(sum(bgroup_corr == k)),k, n) = corr_matrix(bgroup_corr == k);
        
        bdata_bgroup_PCA(1:sum(sum(bgroup_corr == k)),k, n) = corr_matrix_PCA(bgroup_corr == k);
    end
    
end


colors_aux = [c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle];
colors = colors_aux(end:-1:1,:);


for n = 1:nr_networks+1
    
    fig = figure('pos', [-400, 1000, 2000, 600]);
    subplot(2,6,1)
    corr_matrix = corr_networks(:,:,n);
    imagesc(corr_matrix);
    title(networks{n})
   
    subplot(2,6,7)
    corr_matrix = corr_networks_PCA(:,:,n);
    imagesc(corr_matrix);
    title(networks{n})

    subplot(2,6,2)
    boxplot(bdata_wgroup(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wgroup(:,:,n), 5);
    hold off;
    title('within group')
    
    subplot(2,6,8)
    boxplot(bdata_wgroup_PCA(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wgroup_PCA(:,:,n), 5);
    hold off;
    title('within group')
     
    subplot(2,6,3)
    y1 = 0.2;
    y2 = 0.6;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    hold on;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
    boxplot(bdata_wgroup_bsession(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wgroup_bsession(:,:,n), 5)
    ylim([0.2 1]);
    hold off;
    title('w/ group, b/ sessions')
    
    subplot(2,6,9)
    y1 = 0.2;
    y2 = 0.6;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    hold on;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
    boxplot(bdata_wgroup_bsession_PCA(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wgroup_bsession_PCA(:,:,n), 5)
    ylim([0.2 1]);
    hold off;
    title('w/ group, b/ sessions')

    subplot(2,6,4)  
    y1 = 0.2;
    y2 = 0.6;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    hold on;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
    boxplot(bdata_wsubject_bsession(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wsubject_bsession(:,:,n), 5)
    ylim([0.2 1])
    hold off;
    title('w/ subject, b/ sessions')

    subplot(2,6,10)  
    y1 = 0.2;
    y2 = 0.6;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    hold on;
    patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
    patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
    boxplot(bdata_wsubject_bsession_PCA(:,:,n))
    hold on;
    plot_bp_pnts(bdata_wsubject_bsession_PCA(:,:,n), 5)
    ylim([0.2 1])
    hold off;
    title('w/ subject, b/ sessions')

    subplot(2,6,5)
    boxplot(bdata_wsession(:,:,n))
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    plot_bp_pnts(bdata_wsession(:,:,n), 5)
    title('within session')

    subplot(2,6,11)
    boxplot(bdata_wsession_PCA(:,:,n))
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    plot_bp_pnts(bdata_wsession_PCA(:,:,n), 5)
    title('within session')

    
    subplot(2,6,6)
    boxplot(bdata_bgroup(:,:,n))
    hold on;
    g = findobj(gca,'Tag','Box');
    for j=1:length(g)
        patch(get(g(j),'XData'),get(g(j),'YData'),colors(j+2,:),'FaceAlpha',.5);
    end
    plot_bp_pnts(bdata_bgroup(:,:,n), 5)
    title('between groups')
    
    subplot(2,6,12)
    boxplot(bdata_bgroup_PCA(:,:,n))
    hold on;
    g = findobj(gca,'Tag','Box');
    for j=1:length(g)
        patch(get(g(j),'XData'),get(g(j),'YData'),colors(j+2,:),'FaceAlpha',.5);
    end
    plot_bp_pnts(bdata_bgroup_PCA(:,:,n), 5)
    title('between groups')
  

end    

%% Iself

% within subject, between sessions

figure; imagesc(bdata_wsubject_bsession(1:14, :, 10))

figure('pos', [100 200 1200 800]);
suptitle({'Iself (Tepper) - all session combinations', ' ', ' '})
for n=1:nr_networks+1
    bdata_wsubject_bsession_all = [bdata_wsubject_bsession(1:10,1,n); bdata_wsubject_bsession(1:10,2,n); bdata_wsubject_bsession(1:10,3,n); bdata_wsubject_bsession(1:10,4, n); bdata_wsubject_bsession(1:10,5, n); ...
        bdata_wsubject_bsession(1:10,6, n)];
    
    subplot(2,5,n)
    boxplot([bdata_wsubject_bsession_all, [bdata_wsubject_bsession(1:14,7,n); NaN*ones(46,1)]]);
    hold on;
    plot_bp_pnts([bdata_wsubject_bsession_all, [bdata_wsubject_bsession(1:14,7,n); NaN*ones(46,1)]], 10);
    title(networks{n})
    set(gca, 'XTickLabel', {'MIG', 'HC'})
    ylim([0 1])
    grid minor
end

figure('pos', [100 200 1200 800]);
suptitle({'Iself (Tepper) - average', ' ', ' '})
for n=1:nr_networks+1
    bdata_wsubject_bsession_all = mean([bdata_wsubject_bsession(1:10,1,n), bdata_wsubject_bsession(1:10,2,n), bdata_wsubject_bsession(1:10,3,n), bdata_wsubject_bsession(1:10,4, n), bdata_wsubject_bsession(1:10,5, n), ...
        bdata_wsubject_bsession(1:10,6, n)],2);
    
    subplot(2,5,n)
    boxplot([[bdata_wsubject_bsession_all; NaN; NaN; NaN; NaN], bdata_wsubject_bsession(1:14,7,n)]);
    hold on;
    plot_bp_pnts([[bdata_wsubject_bsession_all; NaN; NaN; NaN; NaN], bdata_wsubject_bsession(1:14,7,n)], 10);
    title(networks{n})
    set(gca, 'XTickLabel', {'MIG', 'HC'})
    ylim([0 1])
    grid minor
end

%% Deviation from healthy DevHC
colors_aux = [c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle];
colors = colors_aux(end:-1:1,:);

devHC = zeros(nr_subjects, 2, nr_networks+1);
for n=1:nr_networks+1
    for s = 1:nr_subjects
        corr_matrix = corr_networks(:,:,n);
        
        ind_1 = setdiff(41:2:nr_subjects, s);
        ind_2 = setdiff(42:2:nr_subjects, s);
        
        devHC(s,1,n) =1-mean(corr_matrix(ind_1,s));
        devHC(s,2,n) =1-mean(corr_matrix(ind_2,s));
    end
end


figure('pos', [-100 500 1900 900]);
%figure('pos', [50 50 900 300]);
suptitle({'DevfHC (Tepper)', ' ', ' '})
for n=1:nr_networks+1
    devHC_network_premenstrual = devHC(:,1,n);
    devHC_network_midcycle = devHC(:,2,n);
    
    subplot(2,5,n)
    boxplot([[devHC_network_premenstrual(1:4:40); NaN; NaN; NaN; NaN] [devHC_network_premenstrual(2:4:40); NaN; NaN; NaN; NaN]  [devHC_network_premenstrual(3:4:40); NaN; NaN; NaN; NaN] ...
        [devHC_network_premenstrual(4:4:40); NaN; NaN; NaN; NaN] devHC_network_premenstrual(41:2:68) devHC_network_premenstrual(42:2:68) [devHC_network_midcycle(1:4:40); NaN; NaN; NaN; NaN] ...
        [devHC_network_midcycle(2:4:40); NaN; NaN; NaN; NaN] [devHC_network_midcycle(3:4:40); NaN; NaN; NaN; NaN] [devHC_network_midcycle(4:4:40); NaN; NaN; NaN; NaN] devHC_network_midcycle(41:2:68)...
        devHC_network_midcycle(42:2:68)])
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    
    hold on;
    plot_bp_pnts([[devHC_network_premenstrual(1:4:40); NaN; NaN; NaN; NaN] [devHC_network_premenstrual(2:4:40); NaN; NaN; NaN; NaN]  [devHC_network_premenstrual(3:4:40); NaN; NaN; NaN; NaN] ...
        [devHC_network_premenstrual(4:4:40); NaN; NaN; NaN; NaN] devHC_network_premenstrual(41:2:68) devHC_network_premenstrual(42:2:68) [devHC_network_midcycle(1:4:40); NaN; NaN; NaN; NaN] ...
        [devHC_network_midcycle(2:4:40); NaN; NaN; NaN; NaN] [devHC_network_midcycle(3:4:40); NaN; NaN; NaN; NaN] [devHC_network_midcycle(4:4:40); NaN; NaN; NaN; NaN] devHC_network_midcycle(41:2:68)...
        devHC_network_midcycle(42:2:68)], 10)
    line([6.5 6.5], [0.20 0.8], 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1.2)
    ylim([0.20 0.8])
    grid minor
    set(gca, 'XTick', [4 9])
    set(gca, 'XTickLabel', {'premenstrual', 'midcycle'})
    set(gca, 'FontSize', 12)
    title(networks{n})
end


colors_aux = [c_preictal; c_ictal; c_postictal; c_premenstrual; c_interictal; c_midcycle];
colors = colors_aux(end:-1:1,:);

figure('pos', [-100 500 1900 900]);
%figure('pos', [50 50 900 300]);
suptitle({'DevfHC (Tepper)', ' ', ' '})
for n=1:nr_networks+1
    devHC_network_premenstrual = devHC(:,1,n);
    devHC_network_midcycle = devHC(:,2,n);
    
    subplot(2,5,n)
    boxplot([[devHC_network_premenstrual(1:4:40); NaN; NaN; NaN; NaN] [devHC_network_premenstrual(2:4:40); NaN; NaN; NaN; NaN]  [devHC_network_premenstrual(3:4:40); NaN; NaN; NaN; NaN] ...
        devHC_network_premenstrual(41:2:68) [devHC_network_midcycle(4:4:40); NaN; NaN; NaN; NaN] devHC_network_midcycle(42:2:68)])
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    
    hold on;
    plot_bp_pnts([[devHC_network_premenstrual(1:4:40); NaN; NaN; NaN; NaN] [devHC_network_premenstrual(2:4:40); NaN; NaN; NaN; NaN]  [devHC_network_premenstrual(3:4:40); NaN; NaN; NaN; NaN] ...
        devHC_network_premenstrual(41:2:68) [devHC_network_midcycle(4:4:40); NaN; NaN; NaN; NaN] devHC_network_midcycle(42:2:68)], 10)
    line([4.5 4.5], [0.20 0.8], 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1.2)
    ylim([0.20 0.8])
    grid minor
    set(gca, 'XTick', [2.5 5.5])
    set(gca, 'XTickLabel', {'premenstrual', 'midcycle'})
    set(gca, 'FontSize', 12)
    title(networks{n})
end


n = 8;
[p,h,stats] = ranksum(devHC(4:4:40,1,n),devHC(42:2:68,2,n), 'tail', 'right')
%% Iothers

Iothers_MIG = zeros(40, nr_networks+1);
for n=1:nr_networks+1
    for s = 1:40
        corr_matrix = corr_networks(:,:,n);
        
        if ismember(s, [1:4:40])
            e = 1;
        elseif ismember(s, [2:4:40])
            e = 2;
        elseif ismember(s, [3:4:40])
            e = 3;
        elseif ismember(s, [4:4:40])
            e = 4;
        end
        ind = setdiff([e:4:40],s);
        Iothers_MIG(s,n) = mean(corr_matrix(ind, s));
            
    end
end

Iothers_HC = zeros(28, nr_networks+1);
for n=1:nr_networks+1
    for s = 41:nr_subjects
        corr_matrix = corr_networks(:,:,n);
        if ismember(s, [41:2:nr_subjects])
            e = 41;
        elseif ismember(s, [42:2:nr_subjects])
            e = 42;
        end
        ind = setdiff([e:2:nr_subjects],s);
        Iothers_HC(s-40,n) = mean(corr_matrix(ind, s));

    end
end

colors_aux = [c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle];
colors = colors_aux(end:-1:1,:);

figure;
suptitle({'Iothers (Tepper)', ' ', ' '})
for n = 1:nr_networks+1
    subplot(2,5,n)
    boxplot([[Iothers_MIG(1:4:40,n); NaN*ones(4,1)] [Iothers_MIG(2:4:40,n); NaN*ones(4,1)]...
        [Iothers_MIG(3:4:40,n); NaN*ones(4,1)] [Iothers_MIG(4:4:40,n); NaN*ones(4,1)] Iothers_HC(1:2:28,n) Iothers_HC(2:2:28,n)])
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    plot_bp_pnts([[Iothers_MIG(1:4:40,n); NaN*ones(4,1)] [Iothers_MIG(2:4:40,n); NaN*ones(4,1)]...
        [Iothers_MIG(3:4:40,n); NaN*ones(4,1)] [Iothers_MIG(4:4:40,n); NaN*ones(4,1)]...
        Iothers_HC(1:2:28,n) Iothers_HC(2:2:28,n)], 10)
    title(networks{n})
    set(gca, 'XTickLabel', [])
    set(gca, 'FontSize', 12)
    grid minor
    ylim([0.2 0.85])
end

%% Edgewise ICC
edgewise_ICC_session_p = zeros((nr_areas-1)*nr_areas/2, 1);
edgewise_ICC_session_c = zeros((nr_areas-1)*nr_areas/2, 1);
edgewise_ICC_subject_p = zeros((nr_areas-1)*nr_areas/2, 1);
edgewise_ICC_subject_c = zeros((nr_areas-1)*nr_areas/2, 1);

p_ICC_session_p = zeros((nr_areas-1)*nr_areas/2, 1);
p_ICC_session_c = zeros((nr_areas-1)*nr_areas/2, 1);
p_ICC_subject_p = zeros((nr_areas-1)*nr_areas/2, 1);
p_ICC_subject_c = zeros((nr_areas-1)*nr_areas/2, 1);

for edge = 1:(nr_areas-1)*nr_areas/2
  
    data_session_p = [group_data(edge,1:4:40); group_data(edge,2:4:40); group_data(edge,3:4:40); group_data(edge,4:4:40)]';
    data_subject_p = [group_data(edge,1:4:40); group_data(edge,2:4:40); group_data(edge,3:4:40); group_data(edge,4:4:40)];
    
    data_session_c = [group_data(edge,41:2:68); group_data(edge,42:2:68)]';
    data_subject_c = [group_data(edge,41:2:68); group_data(edge,42:2:68)];

    
    [r_session_p, ~, ~, ~, ~, ~, p_session_p] = ICC(data_session_p, '1-1', 0.05, 0);
    edgewise_ICC_session_p(edge,1) = r_session_p;
    p_ICC_session_p(edge,1) = p_session_p;
    
    [r_subject_p, ~, ~, ~, ~, ~, p_subject_p] = ICC(data_subject_p, '1-1', 0.05, 0);
    edgewise_ICC_subject_p(edge,1) = r_subject_p;
    p_ICC_subject_p(edge,1) = p_subject_p;
    
    [r_session_c, ~, ~, ~, ~, ~, p_session_c] = ICC(data_session_c, '1-1', 0.05, 0);
    edgewise_ICC_session_c(edge,1) = r_session_c;
    p_ICC_session_c(edge,1) = p_session_c;
    
    [r_subject_c, ~, ~, ~, ~, ~, p_subject_c] = ICC(data_subject_c, '1-1', 0.05, 0);
    edgewise_ICC_subject_c(edge,1) = r_subject_c;
    p_ICC_subject_c(edge,1) = p_subject_c;
end

% data_in = edgewise_ICC_subject_c;
% figure;
% aux = tril(ones(nr_areas),-1);
% aux(aux > 0) = data_in;
% data_out = (aux + aux')./(eye(nr_areas)+1);
% 
% figure;
% imagesc(tril(data_out))
% 
% figure;
% imagesc(triu(data_out)')
%%
m_edgewise_ICC_session_p =  conversion_vecmat(edgewise_ICC_session_p, nr_areas, 'vec2mat');
m_edgewise_ICC_session_c =  conversion_vecmat(edgewise_ICC_session_c, nr_areas, 'vec2mat');
m_edgewise_ICC_subject_p =  conversion_vecmat(edgewise_ICC_subject_p, nr_areas, 'vec2mat');
m_edgewise_ICC_subject_c =  conversion_vecmat(edgewise_ICC_subject_c, nr_areas, 'vec2mat');

figure;
imagesc(m_edgewise_ICC_subject_c)
m_p_ICC_session_p =  conversion_vecmat(mafdr(p_ICC_session_p), nr_areas, 'vec2mat');
m_p_ICC_session_c =  conversion_vecmat(mafdr(p_ICC_session_c), nr_areas, 'vec2mat');
m_p_ICC_subject_p =  conversion_vecmat(mafdr(p_ICC_subject_p), nr_areas, 'vec2mat');
m_p_ICC_subject_c =  conversion_vecmat(mafdr(p_ICC_subject_c), nr_areas, 'vec2mat');


figure;
subplot(2,2,1)
%m_edgewise_ICC_session_p(m_edgewise_ICC_session_p<prctile(m_edgewise_ICC_session_p(:), 95)) = NaN;
%m_edgewise_ICC_session_p(session_p) = NaN;
imAlpha=ones(size(m_edgewise_ICC_session_p));
imAlpha(isnan(m_edgewise_ICC_session_p))=0;
%imagesc(m_edgewise_ICC_session_p,'AlphaData',imAlpha)
imagesc(m_edgewise_ICC_session_p)
%cmap = gray(256);
%colormap(cmap)
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('session-patient')

subplot(2,2,2)
%m_edgewise_ICC_session_c(m_edgewise_ICC_session_c<prctile(m_edgewise_ICC_session_c(:), 95)) = NaN;
imAlpha=ones(size(m_edgewise_ICC_session_c));
imAlpha(isnan(m_edgewise_ICC_session_c))=0;
%imagesc(m_edgewise_ICC_session_c,'AlphaData',imAlpha)
imagesc(m_edgewise_ICC_session_c)
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('session-control')

subplot(2,2,3)
%m_edgewise_ICC_subject_p(m_edgewise_ICC_subject_p<prctile(m_edgewise_ICC_subject_p(:), 95)) = NaN;
imAlpha=ones(size(m_edgewise_ICC_subject_p));
imAlpha(isnan(m_edgewise_ICC_subject_p))=0;
%imagesc(m_edgewise_ICC_subject_p,'AlphaData',imAlpha)
imagesc(m_edgewise_ICC_subject_p)
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('subject-patient')

subplot(2,2,4)
%m_edgewise_ICC_subject_c(m_edgewise_ICC_subject_c<prctile(m_edgewise_ICC_subject_c(:), 95)) = NaN;
imAlpha=ones(size(m_edgewise_ICC_subject_c));
imAlpha(isnan(m_edgewise_ICC_subject_c))=0;
%imagesc(m_edgewise_ICC_subject_c, 'AlphaData',imAlpha)
imagesc(m_edgewise_ICC_subject_c)
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('subject-control')

%%
figure;
subplot(2,2,1)
imagesc(1-m_p_ICC_session_p, [0.95 1])
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)

subplot(2,2,2)
imagesc(1-m_p_ICC_session_c, [0.95 1])
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)

subplot(2,2,3)
imagesc(1-m_p_ICC_subject_p, [0.95 1])
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)

subplot(2,2,4)
imagesc(1-m_p_ICC_subject_c, [0.95 1])
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)


session_p = 1-m_p_ICC_session_p;
session_p(session_p<0.95) = 1;

%% PCA reconstruction - 19 PCs
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

nComp = 19;
mean_data = mean(data);
PCA_group_data = data-mean_data;
[coeffs, score, latent] = pca(PCA_group_data, 'NumComponents', nComp);
PCA_matrix = score * coeffs';
group_PCA_recon = bsxfun(@plus, PCA_matrix, mean_data);
 
%group_PCA_recon = group_PCA_recon;
min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1700 700])
ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 1:4:40
    subplot(4,10,c)
    imagesc(conversion_vecmat(group_PCA_recon(:,k), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('preictal')
    end
    sub = subjects{c};         
    title(sub(end-2:end))
        
    subplot(4,10,c+10)
    imagesc(conversion_vecmat(group_PCA_recon(:,k+1), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('ictal')
    end
        
    subplot(4,10,c+20)
    imagesc(conversion_vecmat(group_PCA_recon(:,k+2), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('postictal')
    end
    
    subplot(4,10,c+30)
    imagesc(conversion_vecmat(group_PCA_recon(:,k+3), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('interictal')
    end
    
    c= c+1;
end

figure('pos', [50 50 800 700])
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,1:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,6)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,2:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,11)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,3:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,16)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,4:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)


% c = 1;
% figure('pos', [50 50 1700 700])
% ht = suptitle({['40 PCs - Pearson correlation - Controls (complete, N=14) - ', ...
%     atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
%     ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
% ht.Interpreter = 'none';
% for k = 1:4:40
%     subplot(4,10,c)
%     imagesc(conversion_vecmat(group_PCA_recon(:,k), nr_areas, 'vec2mat'), [min_val max_val])
%     hold on;
%     plot_atlas_labels(atlas, nr_areas);
%     set(gca, 'XTickLabel', []);
%     set(gca, 'YTickLabel', []);
%     set(gca, 'FontSize', 14)
%     if k==1
%         ylabel('preictal')
%     end
%     sub = subjects{c};         
%     title(sub(end-2:end))
% end

%% Original data
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1700 700])
ht = suptitle({['Original - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 1:4:40
    subplot(4,10,c)
    imagesc(conversion_vecmat(group_data(:,k), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('preictal')
    end
    sub = subjects{c};         
    title(sub(end-2:end))
        
    subplot(4,10,c+10)
    imagesc(conversion_vecmat(group_data(:,k+1), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('ictal')
    end
        
    subplot(4,10,c+20)
    imagesc(conversion_vecmat(group_data(:,k+2), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('postictal')
    end
    
    subplot(4,10,c+30)
    imagesc(conversion_vecmat(group_data(:,k+3), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('interictal')
    end
    
    c= c+1;
end

figure('pos', [50 50 800 700])
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_data(:,1:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,6)
imagesc(conversion_vecmat(mean(group_data(:,2:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,11)
imagesc(conversion_vecmat(mean(group_data(:,3:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,16)
imagesc(conversion_vecmat(mean(group_data(:,4:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

%% Controls - original data

figure('pos', [50 50 800 700])
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_data(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,6)
imagesc(conversion_vecmat(mean(group_data(:,42:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

%% Controls - PCA

subjects = {'sub-control019','sub-control020','sub-control026','sub-control027','sub-control028','sub-control029',...
    'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 2300 700])
ht = suptitle({['40 PCs - Pearson correlation - Controls (complete, N=14) - ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 41:2:68
    subplot(4,14,c)
    imagesc(conversion_vecmat(group_data(:,k), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('premenstrual')
    end
    sub = subjects{c};         
    title(sub(end-2:end))
        
    subplot(4,14,c+14)
    imagesc(conversion_vecmat(group_data(:,k+1), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('midcycle')
    end
    
    c= c+1;
end

figure('pos', [50 50 800 700])
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

subplot(4,5,6)
imagesc(conversion_vecmat(mean(group_PCA_recon(:,42:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)

%% Comparisons
figure('pos', [50 50 800 700])
ht = suptitle({['AVERAGE - Original - Pearson correlation ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_data(:,1:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('preictal')

subplot(4,5,2)
imagesc(conversion_vecmat(mean(group_data(:,2:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('ictal')

subplot(4,5,3)
imagesc(conversion_vecmat(mean(group_data(:,3:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('postictal')

subplot(4,5,4)
imagesc(conversion_vecmat(mean(group_data(:,4:4:40),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('interictal')

subplot(4,5,7)
imagesc(conversion_vecmat(mean(group_data(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('premenstrual')

subplot(4,5,9)
imagesc(conversion_vecmat(mean(group_data(:,42:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('midcycle')

%%
min_val = -0.5;
max_val = 0.5;

figure('pos', [50 50 800 700])
ht = suptitle({['DIFF - Original - Pearson correlation ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
subplot(4,5,1)
imagesc(conversion_vecmat(mean(group_data(:,1:4:40),2)-mean(group_data(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('preic-prem')

subplot(4,5,2)
imagesc(conversion_vecmat(mean(group_data(:,2:4:40),2)-mean(group_data(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('ic-prem')

subplot(4,5,3)
imagesc(conversion_vecmat(mean(group_data(:,3:4:40),2)-mean(group_data(:,41:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('postic-prem')

subplot(4,5,4)
imagesc(conversion_vecmat(mean(group_data(:,4:4:40),2)-mean(group_data(:,42:2:68),2), nr_areas, 'vec2mat'), [min_val max_val])
hold on;
plot_atlas_labels(atlas, nr_areas);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 14)
title('interic-mid')


%% Distribution
figure;
histogram(group_PCA_recon(:,1:40))
hold on;
histogram(group_data(:,1:40))
legend('40 PCs', 'Original')


%% Distribution - histograms
% PCA reconstruction - 19 PCs
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};

nComp = 19;
mean_data = mean(data);
PCA_group_data = data-mean_data;
[coeffs, score, latent] = pca(PCA_group_data, 'NumComponents', nComp);
PCA_matrix = score * coeffs';
group_PCA_recon = bsxfun(@plus, PCA_matrix, mean_data);
 
%group_PCA_recon = group_PCA_recon;
min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1700 700])
ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 1:4:40
    subplot(4,10,c)
    histogram(group_PCA_recon(:,k), 20)
    %set(gca, 'XTickLabel', []);
    %set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    xlim([-1 1])
    if k==1
        ylabel('preictal')
    end
    sub = subjects{c};         
    title(sub(end-2:end))
        
    subplot(4,10,c+10)
    histogram(group_PCA_recon(:,k+1), 20)
    %set(gca, 'XTickLabel', []);
    %set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    xlim([-1 1])
    if k==1
        ylabel('ictal')
    end
        
    subplot(4,10,c+20)
    histogram(group_PCA_recon(:,k+2), 20)
    %set(gca, 'XTickLabel', []);
    %set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    xlim([-1 1])
    if k==1
        ylabel('postictal')
    end
    
    subplot(4,10,c+30)
    histogram(group_PCA_recon(:,k+3), 20)
    %set(gca, 'XTickLabel', []);
    %set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    xlim([-1 1])
    if k==1
        ylabel('interictal')
    end
    
    c= c+1;
end
