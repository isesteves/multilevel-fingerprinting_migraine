DI_FC_PCA_path = '/home/iesteves/FC/data/results/DI/DI_FC_PCA';
nbs_stats_path = '/home/iesteves/FC/data/results/nbs/stats/';

% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;

colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

%%
nr_patients = 10;
nr_controls = 14; 
nr_sessions = 6;
ncomp = 68;
contrast_type = '-1';
method = 'Extent';
threshold = 4;
    
load([DI_FC_PCA_path, '/FCvec_PCA_recon.mat'])
FCvec = FCvec_PCA_recon(:,:,ncomp);

group_data = FCvec;

%% ICC
nr_areas = 130;

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
   
    data_session_c = [group_data(edge,41:2:68); group_data(edge,42:2:68)]';
  
    [r_session_p, ~, ~, ~, ~, ~, p_session_p] = ICC(data_session_p, '1-1', 0.05, 0);
    edgewise_ICC_session_p(edge,1) = r_session_p;
    p_ICC_session_p(edge,1) = p_session_p;
    
    [r_session_c, ~, ~, ~, ~, ~, p_session_c] = ICC(data_session_c, '1-1', 0.05, 0);
    edgewise_ICC_session_c(edge,1) = r_session_c;
    p_ICC_session_c(edge,1) = p_session_c;
    
end

%%
m_edgewise_ICC_session_p =  conversion_vecmat(edgewise_ICC_session_p, nr_areas, 'vec2mat');
m_edgewise_ICC_session_c =  conversion_vecmat(edgewise_ICC_session_c, nr_areas, 'vec2mat');
m_edgewise_ICC_subject_p =  conversion_vecmat(edgewise_ICC_subject_p, nr_areas, 'vec2mat');
m_edgewise_ICC_subject_c =  conversion_vecmat(edgewise_ICC_subject_c, nr_areas, 'vec2mat');

figure;
ax = subplot(1,2,1);
m_edgewise_ICC_session_p(m_edgewise_ICC_session_p<prctile(m_edgewise_ICC_session_p(:), 95)) = NaN;
imAlpha=ones(size(m_edgewise_ICC_session_p));
imAlpha(isnan(m_edgewise_ICC_session_p))=0;
imagesc(m_edgewise_ICC_session_p,'AlphaData',imAlpha)
%imagesc(m_edgewise_ICC_session_p > 0.6360)
colormap(ax, flipud(gray))
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('session-patient')
caxis([0 1])
%colorbar

ax1 = subplot(1,2,2);
m_edgewise_ICC_session_c(m_edgewise_ICC_session_c<prctile(m_edgewise_ICC_session_c(:), 95)) = NaN;
imAlpha=ones(size(m_edgewise_ICC_session_c));
imAlpha(isnan(m_edgewise_ICC_session_c))=0;
imagesc(m_edgewise_ICC_session_c,'AlphaData',imAlpha)
%imagesc(m_edgewise_ICC_session_c)
hold on
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('session-control')
caxis([0 1])
colormap(ax1, flipud(gray))
%colorbar

%% ICC
nr_areas = 130;

edgewise_ICC_session_p = zeros((nr_areas-1)*nr_areas/2, 6);

p_ICC_session_p = zeros((nr_areas-1)*nr_areas/2, 6);

comparisons = [1 4; 2 4; 3 4; 1 2; 1 3; 2 3];

for edge = 1:(nr_areas-1)*nr_areas/2
    for c = 1:length(comparisons)
        group1 = group_data(edge,comparisons(c,1):4:40);
        group2 = group_data(edge,comparisons(c,2):4:40);
        data_session_p = [group1; group2]';

        [r_session_p, ~, ~, ~, ~, ~, p_session_p] = ICC(data_session_p, '1-1', 0.05, 0);
        edgewise_ICC_session_p(edge,c) = r_session_p;
        p_ICC_session_p(edge,c) = p_session_p;

    end
end

m_edgewise_ICC_session_p =  conversion_vecmat(edgewise_ICC_session_p, nr_areas, 'vec2mat');

title_labels = {'M-pre vs M-inter', 'M-ict vs M-inter', 'M-post vs M-inter', 'M-pre vs M-ict', 'M-pre vs M-post', 'M-ict vs M-post'};

figure('pos', [50 50 1600 900]);
for c = 1:length(comparisons)
    
    m_edgewise_ICC = m_edgewise_ICC_session_p(:,:,c);
    
    ax = subplot(2,4,c);
    %m_edgewise_ICC(m_edgewise_ICC<prctile(m_edgewise_ICC(:), 95)) = NaN;
    %imAlpha=ones(size(m_edgewise_ICC));
    %imAlpha(isnan(m_edgewise_ICC))=0;
    %imagesc(m_edgewise_ICC,'AlphaData',imAlpha)
    imagesc(m_edgewise_ICC)
    %colormap(ax, flipud(gray))
    hold on
    plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
    title(title_labels{c})
    caxis([0.25 0.95])
    %colorbar
end


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