fig_path = '/home/iesteves/dFC/network_metrics';
addpath('/home/iesteves/dFC/2019_03_03_BCT');

% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;

cmap = redbluecmap;
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);


pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB', 'WB'};

nr_patients = 10;
nr_controls = 14;
nr_pcs = 19;

% load Pearson correlation variables (iFC_all - areas x areas x subject_session, ordered by session type)
load(['/home/iesteves/dFC/sFC/MIGcomplete_nonzero-parcels-50/Pearson-results-patient_control-', atlas,'-', pipeline,'.mat'])
nr_areas = size(iFC_all, 1);
nr_subjects = size(iFC_all,3);

% get the vectorized lower triangular matrix with all sessions for all
% subjects
group_data_persession = zeros((N_areas*N_areas-N_areas)/2, n_Subjects);
for k = 1:n_Subjects
    subject_FC =  iFC_all(:,:,k);  
    group_data_persession(:,k) = conversion_vecmat(subject_FC, nr_areas, 'mat2vec');
end

% reorganize correlation matrices by subject (patients: preictal, ictal, postictal, interictal; controls: premenstrual, midcycle)
order = [1:10:31 2:10:32 3:10:33 4:10:34 5:10:35 6:10:36 7:10:37 8:10:38 9:10:39 10:10:40 55  41 56  42 57  43 58  44 59  45 60  46 61  47 62  48 63  49 64  50 65  51 66  52 67  53 68  54];
filenames_persubject = filenames(order);
group_data = group_data_persession(:, order);

%% Networks matrix
matrix_all = zeros(nr_areas, nr_areas, 9);
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 nr_areas 9];
for g = 1:size(areas,1)
    matrix_all(areas(g,1):areas(g,2),areas(g,1):areas(g,2),g) = g;
    matrix_all(:,:,g) = matrix_all(:,:,g) - diag(diag(matrix_all(:,:,g)));
end
matrix_FC_all = sum(matrix_all,3);

% get the vectorized upper triangular matrix of the networks template
v = conversion_vecmat(matrix_FC_all, nr_areas, 'mat2vec');

% FIGURE - Yeo's networks + subcortical + cerebellum matrix
figure; 
imagesc(matrix_FC_all);
hold on;
plot_atlas_labels('SchaeferSubCRB7100', nr_areas)
title('Yeo networks + subcortical + cerebellum')

%% PCA
nr_networks = 9;
components = [8, 12, 19, 25, 35, 40, 68];
mean_data = mean(group_data);
PCA_group_data = group_data-mean_data;
corr_networks = zeros(size(group_data, 2), size(group_data, 2), nr_networks+1, length(components));
group_data_all = zeros(((nr_areas-1)*nr_areas)/2, nr_subjects, length(components));

c = 1;
% FIGURE - corr matrix for different nb of PCs and networks
figure('pos', [100, 200, 1900 1000]);
suptitle({'Identifiability matrices for PCA reconstructed data (networks + whole brain)', ' ', ' '})
for comp = 1:length(components)
    nComp = components(comp);
    [coeffs, score, ~] = pca(PCA_group_data, 'NumComponents', nComp);
    PCA_matrix = score * coeffs';
    group_PCA_recon = bsxfun(@plus, PCA_matrix, mean_data);
    group_data_all(:,:,comp) = group_PCA_recon;
    
    for n = 1:nr_networks+1
        if n <= nr_networks
            network_data = group_PCA_recon(v==n,:); 
        else
            network_data = group_PCA_recon;
        end

        corr_val = corr(network_data, network_data, 'Type', 'Pearson');
        corr_networks(:,:,n,comp) = corr_val;

        subplot(length(components), nr_networks+1, c)
        imagesc(corr_val)
        if comp == 1
            title(networks{n})
        end
        if n == 1
            ylabel([num2str(nComp), ' PCs'])
        end
        c = c+1;
    end   
end

group_data_PCArecon = group_data_all(:,:,3);

%% Thresholding 

figure; 
for  k = 1:nr_subjects
        subject_FC = group_data_PCArecon(:,k);
        subject_FC_perc80 = prctile(subject_FC, 80);
        subject_FC_perc20 = prctile(subject_FC, 20);
        subplot(8,9,k)
        histogram(subject_FC, 20)
        hold on;
        line([subject_FC_perc80 subject_FC_perc80], ylim)
        line([-subject_FC_perc80 -subject_FC_perc80], ylim)
end


% 20% strongest 

thr = 80;
thr_group_data_PCArecon = group_data_PCArecon;
for k = 1:nr_subjects
    subject_FC = group_data_PCArecon(:,k);
    FC_thr = prctile(subject_FC, thr);
    subject_FC(subject_FC<FC_thr) = 0;
    thr_group_data_PCArecon(:,k) = subject_FC;
end

figure('pos', [50 50 1200 900]); 
for  k = 1:nr_subjects
     subplot(8,9,k)
     imagesc(conversion_vecmat(thr_group_data_PCArecon(:,k), nr_areas, 'vec2mat'), [0 1])
     hold on
     plot_atlas_labels(atlas, nr_areas)
     colormap(newCmap)
     %colorbar
     set(gca, 'XTickLabel', [])
     set(gca, 'YTickLabel', [])

end

%% Patient data
data = thr_group_data_PCArecon;
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};
min_val = -1;
max_val = 1;

nComp = 19;
c = 1;
figure('pos', [50 50 1900 750])
ht = suptitle({['Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 1:4:40
    subplot(4,10,c)
    imagesc(conversion_vecmat(data(:,k), nr_areas, 'vec2mat'), [min_val max_val])
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
    imagesc(conversion_vecmat(data(:,k+1), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('ictal')
    end

    subplot(4,10,c+20)
    imagesc(conversion_vecmat(data(:,k+2), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('postictal')
    end

    subplot(4,10,c+30)
    imagesc(conversion_vecmat(data(:,k+3), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==1
        ylabel('interictal')
    end
    colormap(newCmap)

    c= c+1;
end

%% Control data
subjects = {'sub-control019','sub-control020','sub-control026','sub-control027','sub-control028','sub-control029',...
    'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1800 700])
ht = suptitle({['Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation - Controls (complete, N=14) - ', ...
atlas, ', ', nr_areas, ' ROIs - ', pipeline,...
' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for k = 41:2:68
    subplot(5,14,c)
    imagesc(conversion_vecmat(data(:,k), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if k==41
        ylabel('premenstrual')
    end
    sub = subjects{c};         
    title(sub(end-2:end))

    subplot(5,14,c+14)
    imagesc(conversion_vecmat(data(:,k+1), nr_areas, 'vec2mat'), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 13)
    if k==41
        ylabel('midcycle')
    end
    colormap(newCmap)

    c= c+1;
end


%% Binarization

%% Network measures
data = thr_group_data_PCArecon;
 
metric_group_nodal = zeros(nr_areas, nr_subjects);
metric_group_global = zeros(1, nr_subjects);
metric_names = {'nodedegree', 'clustercoef', 'betweencentr', 'charpathlength', 'globaleff', 'transitivity'};
for m = 1:length(metric_names)
    metric_name = metric_names{m};
    for k = 1:nr_subjects
        subject_FC = conversion_vecmat(data(:,k), nr_areas, 'vec2mat');
        
        switch metric_name
            case 'nodedegree'
                %metric = degrees_und(subject_FC)';
                metric = strengths_und(subject_FC)';
            case 'clustercoef'
                metric = clustering_coef_wu(subject_FC);
            case 'betweencentr'   
                metric = betweenness_wei(subject_FC)';
            case 'charpathlength'
                metric = charpath(subject_FC);
            case 'globaleff'
                metric = efficiency_bin(subject_FC);
            case 'transitivity'
                metric = transitivity_wu(subject_FC);
            otherwise
                disp('Unrecognized metric')
        end

        if any(strcmp(metric_name, {'nodedegree', 'clustercoef', 'betweencentr'}))
            metric_group_nodal(:,k) = metric;
            metric_group= metric_group_nodal;
        elseif any(strcmp(metric_name, {'charpathlength', 'globaleff', 'transitivity'}))
            metric_group_global(:,k) = metric;
            metric_group = metric_group_global;
        end
    end
    network.(metric_name) = metric_group;
end

%%
T_atlas = readtable('/home/iesteves/dFC/allnodes_SchaeferSubCRB7100_130.txt');


%% local measures figure 
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

%metric_names = {'charpathlength', 'globaleff', 'transitivity'};
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};
regions = [61, 64, 58, 59, 60, 103, 104, 111, 112];
label_aux = T_atlas{regions, 6};

for r = 1:length(regions)
    
    region = regions(r);
    
    figure('pos', [50, 50, 1300 600]);
    h = suptitle({[label_aux{r,1}, ' - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
    h.Interpreter = 'none';
    for m = 1:length(metric_names)

        metric_name = metric_names{m};

        metric_group = network.(metric_name);

        subplot(1,3,m)

        bmetric_group = [metric_group(region,41:2:68)' metric_group(region,42:2:68)',...
            [metric_group(region,1:4:40)'; NaN*ones(4,1)], [metric_group(region, 2:4:40)'; NaN*ones(4,1)], [metric_group(region,3:4:40)'; NaN*ones(4,1)], [metric_group(region,4:4:40)'; NaN*ones(4,1)]];
        boxplot(bmetric_group);
        hold on;
        plot_bp_pnts(bmetric_group, 10);
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        set(gca, 'XTickLabel', {'premenstrual', 'midcycle', 'preictal', 'ictal', 'postictal', 'interictal'});
        set(gca, 'XTickLabelRotation', 40)
        set(gca, 'FontSize', 13)
        %ylim([0.2 1])
        %ylabel(' ')
        grid minor
        title(metric_name)

    end
    
    print([fig_path, '/local_',label_aux{r,1},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

end

%% local measures figure - per network
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

% stats
vec1 = [1, 3, 3, 3, 4, 4, 5, 1, 1, 1, 2];
vec2 = [2, 4, 5, 6, 5, 6, 6, 3, 4, 5, 6];
ncomp = length(vec1);
hvec = zeros(ncomp, nr_networks, length(metric_names));
pvec = zeros(ncomp, nr_networks, length(metric_names));

min_val = [0, 0, 0.1];
max_val = [45, 480, 0.7];

for m = 1:length(metric_names)

    metric_name = metric_names{m};

    metric_group = network.(metric_name);

    figure('pos', [50, 50, 1000 1000]);
    h = suptitle({[metric_name, ' - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
    h.Interpreter = 'none';
    for n = 1:nr_networks

        network_name = networks{n};
        region = T_atlas{:,4} == n;

        subplot(3,3,n)
        bmetric_group = [mean(metric_group(region,41:2:68))' mean(metric_group(region,42:2:68))',...
            [mean(metric_group(region,1:4:40))'; NaN*ones(4,1)], [mean(metric_group(region, 2:4:40))';...
            NaN*ones(4,1)], [mean(metric_group(region,3:4:40))'; NaN*ones(4,1)], [mean(metric_group(region,4:4:40))'; NaN*ones(4,1)]];
        boxplot(bmetric_group);
        hold on;
        plot_bp_pnts(bmetric_group, 10);       
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        %if ismember(n, [7,8,9])
        %    set(gca, 'XTickLabel', {'premenstrual', 'midcycle', 'preictal', 'ictal', 'postictal', 'interictal'});
        %else 
        set(gca, 'XTickLabel', []);
        %end
        set(gca, 'XTickLabelRotation', 40)
        set(gca, 'FontSize', 13)
        ylim([min_val(m) max_val(m)])
        line([2.5 2.5], ylim, 'Color', 'k', 'LineStyle', '-.')
        grid minor
        title([network_name, ' (N = ',num2str(sum(region)),')']);
        
        % stats
        for i = 1:ncomp
            if ncomp < 8 
                [p,h] = signrank(bmetric_group(:,vec1(i)), bmetric_group(:,vec2(i)));
                hvec(i, n, m) = h;
                pvec(i, n, m) = p;
            elseif ncomp >= 8
                [p,h] = ranksum(bmetric_group(:,vec1(i)), bmetric_group(:,vec2(i)));
                hvec(i, n, m) = h;
                pvec(i, n, m) = p;
            end
        end

    end
    print([fig_path, '/local_',metric_name,'_networks_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

end

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec, 0.05, true);  

%% global measures figure
vec1 = [1, 3, 3, 3, 4, 4, 5, 1, 1, 1, 2];
vec2 = [2, 4, 5, 6, 5, 6, 6, 3, 4, 5, 6];
ncomp = length(vec1);
hvec = zeros(ncomp, length(metric_names));
pvec = zeros(ncomp, length(metric_names));

colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

metric_names = {'charpathlength', 'globaleff', 'transitivity'};
region = 1;
figure('pos', [50, 50, 1200 300]);
h = suptitle({['Global measures - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
h.Interpreter = 'none';
for m = 1:length(metric_names)

    metric_name = metric_names{m};
   
    metric_group = network.(metric_name);

    subplot(1,3,m)

    bmetric_group = [metric_group(region,41:2:68)' metric_group(region,42:2:68)',...
        [metric_group(region,1:4:40)'; NaN*ones(4,1)], [metric_group(region, 2:4:40)'; NaN*ones(4,1)], [metric_group(region,3:4:40)'; NaN*ones(4,1)], [metric_group(region,4:4:40)'; NaN*ones(4,1)]];
    boxplot(bmetric_group);
    hold on;
    plot_bp_pnts(bmetric_group, 10);
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    set(gca, 'XTickLabel', {'premenstrual', 'midcycle', 'preictal', 'ictal', 'postictal', 'interictal'});
    set(gca, 'XTickLabelRotation', 40)
    set(gca, 'FontSize', 13)
    line([2.5 2.5], ylim, 'Color', 'k', 'LineStyle', '-.')
    %ylim([0.2 1])
    %ylabel(' ')
    grid minor
    title(metric_name)

    % stats
    for i = 1:ncomp
        if ncomp < 8 
            [p,h] = signrank(bmetric_group(:,vec1(i)), bmetric_group(:,vec2(i)));
            hvec(i, m) = h;
            pvec(i, m) = p;
        elseif ncomp >= 8
            [p,h] = ranksum(bmetric_group(:,vec1(i)), bmetric_group(:,vec2(i)));
            hvec(i, m) = h;
            pvec(i, m) = p;
        end
    end
    
end
print([fig_path, '/global_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(:), 0.05, true);  
%[FDR] = mafdr(pvec(:),'BHFDR',1);

%%
font = 10;
load('/home/iesteves/mign2treat_sample/MIGN2TREATcontrolsample.mat')
load('/home/iesteves/mign2treat_sample/MIGN2TREATpatientsample.mat')
patients = MIGN2TREATpatientsample;
complete_patients_ind = strcmp(patients.status, 'complete');
complete_patients = patients(complete_patients_ind, :);
preictal = complete_patients(strcmp(complete_patients.session, 'ses-preictal'),:);
postictal = complete_patients(strcmp(complete_patients.session, 'ses-postictal'),:);
ictal = complete_patients(strcmp(complete_patients.session, 'ses-ictal'),:);
interictal = complete_patients(strcmp(complete_patients.session, 'ses-interictal'),:);


controls = MIGN2TREATcontrolsample;
complete_controls_ind = strcmp(controls.status, 'complete');
complete_controls = controls(complete_controls_ind, :);
nr_controls = length(unique(complete_controls.subject))-1; % 025 was excluded
premenstrual = complete_controls(strcmp(complete_controls.session, 'ses-premenstrual'),:);
midcycle = complete_controls(strcmp(complete_controls.session, 'ses-midcycle'),:);

%%
metric_names = {'charpathlength', 'globaleff', 'transitivity'};

figure('pos', [50 50 1000 800])
h = suptitle({['PREICTAL + POSTICTAL - Global measures vs Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
h.Interpreter = 'none';
for m = 1:length(metric_names)

    metric_name = metric_names{m};
   
    metric_group = network.(metric_name);

    m_val_preictal = metric_group(1,1:4:40);
    m_val_postictal = metric_group(1,3:4:40);
    m_val_premenstrual = metric_group(1,41:2:68);

    subplot(3,2,2*m-1)

    scatter(-1*preictal.attack_hrs_before, m_val_preictal, 70, 'o', 'MarkerFaceColor', c_preictal, 'MarkerEdgeColor', 'none')
    lsline
    hold on;
    scatter(zeros(nr_controls,1), m_val_premenstrual, 70, 'o', 'MarkerFaceColor', c_premenstrual, 'MarkerEdgeColor', 'none')
    [r,p] = corrcoef(-1*preictal.attack_hrs_before,  m_val_preictal);
    text(mean(-1*preictal.attack_hrs_before),mean( m_val_preictal), ...
        ['r² = ', num2str(round(r(1,2)^2,3)), '; p = ', num2str(round(p(1,2),4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
    if m == 1
        title('preictal')
    end
    xlabel('hrs to attack onset')
    ylabel(metric_name)

    subplot(3,2,2*m)
    scatter(postictal.attack_hrs_since_attack_end, m_val_postictal, 70, 'o', 'MarkerFaceColor', c_postictal, 'MarkerEdgeColor', 'none')
    lsline
    hold on;
    scatter(zeros(nr_controls,1), m_val_premenstrual, 70, 'o', 'MarkerFaceColor', c_premenstrual, 'MarkerEdgeColor', 'none')
    [r,p] = corrcoef(postictal.attack_hrs_since_attack_end, m_val_postictal);
    text(mean(postictal.attack_hrs_since_attack_end),mean(m_val_postictal), ...
        ['r² = ', num2str(round(r(1,2)^2,3)), '; p = ', num2str(round(p(1,2),4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
    if m == 1      
        title('postictal')
    end
    xlabel('hrs to attack onset')
    
    
end
print([fig_path, '/clinical_ses-preictal-ses-postictal_network-global_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')


%% Global measures vs Ictal
metric_names = {'charpathlength', 'globaleff', 'transitivity'};

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_nausea', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

d = 1;
figure('pos', [50 50 1650 800])
h = suptitle({['ICTAL - Global measures vs Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
h.Interpreter = 'none';
for m = 1:length(metric_names)

    metric_name = metric_names{m};
   
    metric_group = network.(metric_name);

    m_val_ictal = metric_group(1,2:4:40);
    m_val_premenstrual = metric_group(1,41:2:68);

    for c = 1:length(clinical_names)
        
        clinical_name = clinical_names{c};
   
        clinical_data = ictal.(clinical_name);
        subplot(3,6,d)

        scatter(clinical_data, m_val_ictal, 70, 'o', 'MarkerFaceColor', c_ictal, 'MarkerEdgeColor', 'none')
        lsline
        hold on;
        scatter(zeros(nr_controls,1), m_val_premenstrual, 70, 'o', 'MarkerFaceColor', c_premenstrual, 'MarkerEdgeColor', 'none')
        [r,p] = corrcoef(clinical_data, m_val_ictal);
        text(mean(clinical_data),mean(m_val_ictal), ...
            ['r² = ', num2str(round(r(1,2)^2,3)), '; p = ', num2str(round(p(1,2),4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
        if m == 1
              ht = title(clinical_name);
              ht.Interpreter = 'none';
        end
        if c == 1
            ylabel(metric_name)
        end
        d = d +1;
        
        if m ==3
            xlabel(clinical_units{c})
        end
    end
end
print([fig_path, '/clinical_ictal_network-global_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

%% Global measures - Session vs Clinical General

thr = 80;
corr_type = 'Spearman';

metric_names = {'charpathlength', 'globaleff', 'transitivity'};

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'yrs', 'yrs', '#/month', 'hrs', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

vec1 = [1, 3, 3, 3, 4, 4, 5, 1, 1, 1, 2];
vec2 = [2, 4, 5, 6, 5, 6, 6, 3, 4, 5, 6];
ncomp = length(vec1);
hvec = zeros(length(sessions), length(metric_names), length(clinical_names));
pvec = zeros(length(sessions), length(metric_names), length(clinical_names));


for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end

    d = 1;
    figure('pos', [50 50 1450 800])
    h = suptitle({[session_patient, ' - Global measures vs GENERAL Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
    h.Interpreter = 'none';
    for m = 1:length(metric_names)

        metric_name = metric_names{m};

        metric_group = network.(metric_name);

        m_val_patient = metric_group(1,ind_patient:4:40);
        m_val_control = metric_group(1,ind_control:2:68);

        for c = 1:length(clinical_names)

            clinical_name = clinical_names{c};

            clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
            subplot(3,5,d)

            scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
            lsline
            hold on;
            scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
            [r,p] = corr(clinical_data, m_val_patient', 'Type', corr_type);
            %[r,p] = corrcoef(clinical_data, m_val_patient);
            text(mean(clinical_data),mean(m_val_patient), ...
                ['r² = ', num2str(round(r^2,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
            if m == 1
                  ht = title(clinical_name);
                  ht.Interpreter = 'none';
            end
            if c == 1
                ylabel(metric_name)
            end
            d = d +1;

            if m ==3
                xlabel(clinical_units{c})
            end
            
            % stats
            pvec(s, m, c) = p;

            
        end
    end
    print([fig_path, '/general-clinical_ses-', session_patient,'_network-global_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
end

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(4,:,:), 0.05, true);  

%% Local measures - Session vs Clinical General
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

regions = [61, 64, 58, 59, 60, 103, 104, 111, 112];
label_aux = T_atlas{regions, 6};

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    for a = 1:length(regions)
    
        region = regions(a);

        d = 1;
        figure('pos', [50 50 1450 800])
        h = suptitle({[session_patient, ' - ',label_aux{a},' (Local) vs GENERAL Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
        h.Interpreter = 'none';
        for m = 1:length(metric_names)

            metric_name = metric_names{m};

            metric_group = network.(metric_name);

            m_val_patient = metric_group(region,ind_patient:4:40);
            m_val_control = metric_group(region,ind_control:2:68);

            for c = 1:length(clinical_names)

                clinical_name = clinical_names{c};

                clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
                subplot(3,5,d)

                scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
                lsline
                hold on;
                scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
                [r,p] = corrcoef(clinical_data, m_val_patient);
                text(mean(clinical_data),mean(m_val_patient), ...
                    ['r² = ', num2str(round(r(1,2)^2,3)), '; p = ', num2str(round(p(1,2),4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
                if m == 1
                      ht = title(clinical_name);
                      ht.Interpreter = 'none';
                end
                if c == 1
                    ylabel(metric_name)
                end
                d = d +1;

                if m ==3
                    xlabel(clinical_units{c})
                end
            end
        end
        print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
   end
end

%% Local measures vs Ictal
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_nausea', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

regions = [61, 64, 58, 59, 60, 103, 104, 111, 112];
label_aux = T_atlas{regions, 6};

for a = 1:length(regions)
    
    region = regions(a);
    
    d = 1;
    figure('pos', [50 50 1650 800])
    h = suptitle({['ses-ictal - ', label_aux{a},' (Local) vs Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs - Pearson correlation'] , '  ', '  '});
    h.Interpreter = 'none';
    for m = 1:length(metric_names)

        metric_name = metric_names{m};

        metric_group = network.(metric_name);

        m_val_ictal = metric_group(region,2:4:40);
        m_val_premenstrual = metric_group(region,41:2:68);

        for c = 1:length(clinical_names)

            clinical_name = clinical_names{c};

            clinical_data = ictal.(clinical_name);
            subplot(3,6,d)

            scatter(clinical_data, m_val_ictal, 70, 'o', 'MarkerFaceColor', c_ictal, 'MarkerEdgeColor', 'none')
            lsline
            hold on;
            scatter(zeros(nr_controls,1), m_val_premenstrual, 70, 'o', 'MarkerFaceColor', c_premenstrual, 'MarkerEdgeColor', 'none')
            [r,p] = corrcoef(clinical_data, m_val_ictal);
            text(mean(clinical_data),mean(m_val_ictal), ...
                ['r² = ', num2str(round(r(1,2)^2,3)), '; p = ', num2str(round(p(1,2),4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
            if m == 1
                  ht = title(clinical_name);
                  ht.Interpreter = 'none';
            end
            if c == 1
                ylabel(metric_name)
            end
            d = d +1;

            if m ==3
                xlabel(clinical_units{c})
            end
        end
    end
    print([fig_path, '/general-clinical_ictal_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')

end

%% Local metrics vs general clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';
nr_sessions = 4;
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

regions = 1:nr_areas; 
%label_aux = T_atlas{regions, 6};

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

regions_corr = zeros(nr_sessions, nr_areas);
regions_p = zeros(nr_sessions, nr_areas);

   
figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Local metrics (all sessions) vs GENERAL Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 

for m = 1:length(metric_names)

    metric_name = metric_names{m};

    metric_group = network.(metric_name);
  
    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;

        for s = 1:length(sessions)
            session_patient = sessions{s};

            if strcmp(session_patient, 'preictal')
                ind_patient = 1;
                c_patient = c_preictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'ictal')
                ind_patient = 2;
                c_patient = c_ictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'postictal')
                ind_patient = 3;
                c_patient = c_postictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'interictal')
                ind_patient = 4;
                c_patient = c_interictal;
                ind_control = 42;
                c_control = c_midcycle;
            end
          
            for a = 1:length(regions)
    
                region = regions(a);

                m_val_patient = metric_group(region,ind_patient:4:40);            
               
                [r,p] = corr(clinical_data, m_val_patient', 'Type', corr_type);
                regions_corr(s, a) = r;
                regions_p(s, a) = p;
            end % regions
        end % sessions
       
         subplot(6,5,d)
         regions_p(isnan(regions_p)) =Inf;
         imagesc(regions_p, [0 0.05])
         colormap('gray')
         %imagesc(regions_corr, [-1 1])
         %colormap(newCmap)
         d = d+1;
         
         if m == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
         end
         if c == 1
            ylabel(metric_name)
            set(gca, 'YTick', 1:4)
            set(gca, 'YTickLabel', {['\bf \color[rgb]{', num2str(c_preictal),'}', 'preic'], ['\bf \color[rgb]{', num2str(c_ictal),'}', 'ic'],...
                ['\bf \color[rgb]{', num2str(c_postictal),'}', 'postic'], ['\bf \color[rgb]{', num2str(c_interictal),'}', 'interic']})
         else
            set(gca, 'YTickLabel', [])
         end
          set(gca, 'XTickLabel', [])
         % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1=areas(a, 1); x2=areas(a, 2);
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 4.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 4.5])
        %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
     end % clinical variables 
end % metrics 

%% Local metrics vs ictal clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';
xmax = nr_areas;
nr_sessions = 4;
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

regions = 1:nr_areas; 
%label_aux = T_atlas{regions, 6};

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

regions_corr = zeros(nr_sessions, nr_areas);
regions_p = zeros(nr_sessions, nr_areas);

   
figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Local metrics (all sessions) vs Ictal Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 

for m = 1:length(metric_names)

    metric_name = metric_names{m};

    metric_group = network.(metric_name);
  
    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;

        for s = 1:length(sessions)
            session_patient = sessions{s};

            if strcmp(session_patient, 'preictal')
                ind_patient = 1;
                c_patient = c_preictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'ictal')
                ind_patient = 2;
                c_patient = c_ictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'postictal')
                ind_patient = 3;
                c_patient = c_postictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'interictal')
                ind_patient = 4;
                c_patient = c_interictal;
                ind_control = 42;
                c_control = c_midcycle;
            end
          
            for a = 1:length(regions)
    
                region = regions(a);

                m_val_patient = metric_group(region,ind_patient:4:40);            
               
                [r,p] = corr(clinical_data, m_val_patient', 'Type', corr_type);
                regions_corr(s, a) = r;
                regions_p(s, a) = p;
            end % regions
        end % sessions
       
         subplot(6,5,d)
         regions_p(isnan(regions_p)) = Inf;
         %imagesc(regions_p, [0 0.05])
         %colormap('gray')
         imagesc(regions_corr, [-1 1])
         colormap(newCmap)
         d = d+1;
         
         if m == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
         end
         if c == 1
            ylabel(metric_name)
            set(gca, 'YTick', 1:4)
            set(gca, 'YTickLabel', {['\bf \color[rgb]{', num2str(c_preictal),'}', 'preic'], ['\bf \color[rgb]{', num2str(c_ictal),'}', 'ic'],...
                ['\bf \color[rgb]{', num2str(c_postictal),'}', 'postic'], ['\bf \color[rgb]{', num2str(c_interictal),'}', 'interic']})
         else
            set(gca, 'YTickLabel', [])
         end
          set(gca, 'XTickLabel', [])
         % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1=areas(a, 1); x2=areas(a, 2);
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 4.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 4.5])
        
        
        %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
     end % clinical variables 
end % metrics 


%% Avg Local metrics per Network vs ictal clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Pearson';
xmax = nr_areas;
nr_sessions = 4;
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

regions = 1:nr_areas; 
%label_aux = T_atlas{regions, 6};

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

regions_corr = zeros(nr_sessions, nr_networks);
regions_p = zeros(nr_sessions, nr_networks);

   
figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Local metrics (all sessions) vs Ictal Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 

for m = 1:length(metric_names)

    metric_name = metric_names{m};

    metric_group = network.(metric_name);
  
    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;

        for s = 1:length(sessions)
            session_patient = sessions{s};

            if strcmp(session_patient, 'preictal')
                ind_patient = 1;
                c_patient = c_preictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'ictal')
                ind_patient = 2;
                c_patient = c_ictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'postictal')
                ind_patient = 3;
                c_patient = c_postictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'interictal')
                ind_patient = 4;
                c_patient = c_interictal;
                ind_control = 42;
                c_control = c_midcycle;
            end
          
            for k = 1:length(areas)
    
                region = regions(areas(k,1):areas(k,2));

                m_val_patient = mean(metric_group(region,ind_patient:4:40));            
               
                [r,p] = corr(clinical_data, m_val_patient', 'Type', corr_type);
                regions_corr(s, k) = r;
                regions_p(s, k) = p;
            end % regions
        end % sessions
       
         subplot(6,5,d)
         regions_p(isnan(regions_p)) = Inf;
         %imagesc(regions_p(:,:), [0 0.05])
         %colormap('gray')
         imagesc(regions_corr(:,:), [-1 1])
         colormap(newCmap)
         d = d+1;
         
         if m == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
         end
         if c == 1
            ylabel(metric_name)
            set(gca, 'YTick', 1:4)
            set(gca, 'YTickLabel', {['\bf \color[rgb]{', num2str(c_preictal),'}', 'preic'], ['\bf \color[rgb]{', num2str(c_ictal),'}', 'ic'],...
                ['\bf \color[rgb]{', num2str(c_postictal),'}', 'postic'], ['\bf \color[rgb]{', num2str(c_interictal),'}', 'interic']})
         else
            set(gca, 'YTickLabel', [])
         end
          set(gca, 'XTickLabel', [])
         % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1= a-0.5; x2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 4.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 4.5])
        
        
        %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
     end % clinical variables 
end % metrics 

%% Avg Local metrics per Network vs general clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';
xmax = nr_areas;
nr_sessions = 4;
metric_names = {'nodedegree', 'betweencentr', 'clustercoef'};

regions = 1:nr_areas; 
%label_aux = T_atlas{regions, 6};

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

regions_corr = zeros(nr_sessions, nr_networks);
regions_p = zeros(nr_sessions, nr_networks);

   
figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Local metrics (all sessions) vs General Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 

for m = 1:length(metric_names)

    metric_name = metric_names{m};

    metric_group = network.(metric_name);
  
    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;

        for s = 1:length(sessions)
            session_patient = sessions{s};

            if strcmp(session_patient, 'preictal')
                ind_patient = 1;
                c_patient = c_preictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'ictal')
                ind_patient = 2;
                c_patient = c_ictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'postictal')
                ind_patient = 3;
                c_patient = c_postictal;
                ind_control = 41;
                c_control = c_premenstrual;
            elseif strcmp(session_patient, 'interictal')
                ind_patient = 4;
                c_patient = c_interictal;
                ind_control = 42;
                c_control = c_midcycle;
            end
          
            for k = 1:length(areas)
    
                region = regions(areas(k,1):areas(k,2));

                m_val_patient = mean(metric_group(region,ind_patient:4:40));            
               
                [r,p] = corr(clinical_data, m_val_patient', 'Type', corr_type);
                regions_corr(s, k) = r;
                regions_p(s, k) = p;
            end % regions
        end % sessions
       
         subplot(6,5,d)
         regions_p(isnan(regions_p)) = Inf;
         %imagesc(regions_p(:,:), [0 0.05])
         %colormap('gray')
         imagesc(regions_corr(:,:), [-1 1])
         colormap(newCmap)
         d = d+1;
         
         if m == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
         end
         if c == 1
            ylabel(metric_name)
            set(gca, 'YTick', 1:4)
            set(gca, 'YTickLabel', {['\bf \color[rgb]{', num2str(c_preictal),'}', 'preic'], ['\bf \color[rgb]{', num2str(c_ictal),'}', 'ic'],...
                ['\bf \color[rgb]{', num2str(c_postictal),'}', 'postic'], ['\bf \color[rgb]{', num2str(c_interictal),'}', 'interic']})
         else
            set(gca, 'YTickLabel', [])
         end
          set(gca, 'XTickLabel', [])
         % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1= a-0.5; x2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 4.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 4.5])
        
        
        %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
     end % clinical variables 
end % metrics 


%% Avg FC vs general clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';

nr_sessions = 4;

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'yrs', 'yrs', '#/month', 'hrs', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

% stats
hvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));
pvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));

figure('pos', [50 50 1100 800])
h = suptitle({[corr_type, ' correlation - Avg network FC vs General Clinical - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
%h = suptitle({[corr_type, ' correlation - Avg network FC vs General Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 
  
for s = 1:nr_sessions
session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end


    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient_all = zeros(nr_networks, nr_networks, nr_patients);
    for a = 1:nr_networks
        for b = 1:nr_networks                
            area = all_data(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_patient:4:40);              
            for k = 1:10
                area_patient = area(:,:,k);
                if a == b
                    aux = tril(ones(length(area_patient)),-1);
                    area_patient(aux == 0) = NaN;
                end
                m_val_patient_all(a, b, k) = nanmean(area_patient(:));                    
            end
        end
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
        
        g = 1;
        regions_corr = zeros(nr_networks, nr_networks);
        regions_p = zeros(nr_networks, nr_networks);
        for a = 1:nr_networks
            for b = 1:nr_networks   
                [r,p] = corr(clinical_data, squeeze(m_val_patient_all(a,b,:)), 'Type', corr_type);
                regions_corr(a, b) = r;
                regions_p(a, b) = p;             
                
                % stats
                pvec(s, g, c) = p;
                g = g+1;
            end
        end % regions

        subplot(4,5,d)
        %regions_p(isnan(regions_p)) = Inf;
        %imagesc(regions_p(:,:), [0 0.05])
        colormap('gray')
        imagesc(regions_corr(:,:), [-1 1])
        colormap(newCmap)
        d = d+1;

        if s == 1
        ht = title(clinical_name);
        ht.Interpreter = 'none';
        end
        if c == 1
        ylabel(['\bf \color[rgb]{', num2str(c_patient),'}', session_patient]);

        end
        set(gca, 'YTickLabel', [])
        set(gca, 'XTickLabel', [])
        % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1= a-0.5; x2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 9.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 9.5]) 

        % vertical
        for a =1:size(areas, 1)

            x1=-10; x2=0.5;
            y1= a-0.5; y2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([0 9.5], [y2 y2], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 9.5])
        xlim([0 9.5])
        
        if s==length(sessions)
            xlabel(clinical_units{c})
        end

    end % clinical variables
     
    %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
end % sessions

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(1,:,:), 0.05, true);  

[dim1, dim2, dim3] = size(pvec(1,:,:));

[i, j, k] = ind2sub([dim1, dim2, dim3], 152);

[i1, j1] = ind2sub([9 9], 71);
%% Avg FC vs ictal clinical features
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';

nr_sessions = 4;

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

% stats
hvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));
pvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));

figure('pos', [50 50 1250 950])
h = suptitle({[corr_type, ' correlation - Avg network FC vs Ictal Clinical - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
%h = suptitle({[corr_type, ' correlation - Avg network FC vs General Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 
  
for s = 1:nr_sessions
session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end


    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient_all = zeros(nr_networks, nr_networks, nr_patients);
    for a = 1:nr_networks
        for b = 1:nr_networks                
            area = all_data(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_patient:4:40);              
            for k = 1:10
                area_patient = area(:,:,k);
                if a == b
                    aux = tril(ones(length(area_patient)),-1);
                    area_patient(aux == 0) = NaN;
                end
                m_val_patient_all(a, b, k) = nanmean(area_patient(:));                    
            end
        end
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;

        g = 1;
        regions_corr = zeros(nr_networks, nr_networks);
        regions_p = zeros(nr_networks, nr_networks);
        for a = 1:nr_networks
            for b = 1:nr_networks   
                [r,p] = corr(clinical_data, squeeze(m_val_patient_all(a,b,:)), 'Type', corr_type);
                regions_corr(a, b) = r;
                regions_p(a, b) = p;
                                
                % stats
                pvec(s, g, c) = p;
                g = g+1;
            end
        end % regions

        subplot(4,5,d)
        regions_p(isnan(regions_p)) = Inf;
        imagesc(regions_p(:,:), [0 0.05])
        colormap('gray')
        %imagesc(regions_corr(:,:), [-1 1])
        %colormap(newCmap)
        d = d+1;

        if s == 1
        ht = title(clinical_name);
        ht.Interpreter = 'none';
        end
        if c == 1
        ylabel(['\bf \color[rgb]{', num2str(c_patient),'}', session_patient]);

        end
        set(gca, 'YTickLabel', [])
        set(gca, 'XTickLabel', [])
        % horizontal
        for a =1:size(areas, 1)

            y1=-10; y2=0.5;
            x1= a-0.5; x2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([x2 x2], [0 9.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 9.5]) 

        % vertical
        for a =1:size(areas, 1)

            x1=-10; x2=0.5;
            y1= a-0.5; y2=a+0.5;
            v=[x1 y1; x1 y2; x2 y2; x2 y1];
            f=[1 2 3 4];
            hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
            h(a) = hi(1);
            hold on;
            line([0 9.5], [y2 y2], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', '-.')
        end
        ylim([0 9.5])
        xlim([0 9.5])
        
        if s==length(sessions)
            xlabel(clinical_units{c})
        end

    end % clinical variables
     
    %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
end % sessions

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec, 0.05, true); 

%% Avg sig edges vs General Clinical
corr_type = 'Spearman';
results_path = '/home/iesteves/dFC/NBS1.2/NBS1.2/mydata/results_pca19';

contrast_type = 1;
method = 'Extent';
threshold = 4;

clinical_names = {'mig_age_onset', 'mig_disease_duration', 'mig_freq_monthly', 'mig_duration', 'mig_pain_intensity'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'ictal', 'postictal'};

% stats
hvec = zeros(length(sessions), length(clinical_names));
pvec = zeros(length(sessions), length(clinical_names));

figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Avg sig FC edges vs General Clinical - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;
for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'ic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'postic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    load([results_path, '/nbs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
    adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
    
%     figure('pos', [50, 50, 300 450]);
%     subplot(3,1,[1 2])
%     imAlpha=ones(size(adj_matrix));
%     imAlpha(adj_matrix==0)=0;
%     imagesc(tril(adj_matrix,0),'AlphaData', imAlpha)
%     hold on;
%     plot_atlas_labels(atlas, nr_areas);
%     cmap = gray(256);
%     colormap(flipud(cmap))
%     set(gca, 'XTickLabel', '')
%     set(gca, 'YTickLabel', '')
%     title(comparison)

    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient = zeros(nr_patients, 1);
    area = all_data(:,:,ind_patient:4:40);
    for k = 1:nr_patients
        area_patient = area(:,:,k);
        area_sig = area_patient(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_patient(k) = val;      
    end
    
    m_val_control = zeros(nr_controls, 1);
    area = all_data(:,:,ind_control:2:68);
    for k = 1:nr_controls
        area_control = area(:,:,k);
        area_sig = area_control(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_control(k) = val;      
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};
       
        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
        subplot(2,5,d)

        scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
        lsline
        hold on;
        scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
        [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
         % stats
        pvec(s, c) = p;
        text(mean(clinical_data),mean(m_val_patient), ...
            ['r² = ', num2str(round(r^2,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
        ylim([0.15 0.75])
        if s == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
        end
        if c == 1
            ylabel(['\bf \color[rgb]{', num2str(c_patient),'}', session_patient]);
        end

        if s == 2
            xlabel(clinical_units{c})
        end
        
        d = d +1;
    end
end

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec(2,:), 0.05, true); 

[i,j] = ind2sub([2,5],6)
%% Avg sig edges vs Ictal Clinical
corr_type = 'Spearman';
results_path = '/home/iesteves/dFC/NBS1.2/NBS1.2/mydata/results_pca19';

contrast_type = 1;
method = 'Extent';
threshold = 4;

clinical_names = {'attack_hrs_since_onset', 'attack_pain', 'attack_photophobia', 'attack_phonophobia', 'attack_motion_intolerance'};
clinical_units = {'hrs', 'score', 'score', 'score', 'score', 'score'};

sessions = {'ictal', 'postictal'};

figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Avg sig FC edges vs General Clinical - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;
for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'ic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'postic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    load([results_path, '/nbs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
    adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
    
%     figure('pos', [50, 50, 300 450]);
%     subplot(3,1,[1 2])
%     imAlpha=ones(size(adj_matrix));
%     imAlpha(adj_matrix==0)=0;
%     imagesc(tril(adj_matrix,0),'AlphaData', imAlpha)
%     hold on;
%     plot_atlas_labels(atlas, nr_areas);
%     cmap = gray(256);
%     colormap(flipud(cmap))
%     set(gca, 'XTickLabel', '')
%     set(gca, 'YTickLabel', '')
%     title(comparison)

    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient = zeros(nr_patients, 1);
    area = all_data(:,:,ind_patient:4:40);
    for k = 1:nr_patients
        area_patient = area(:,:,k);
        area_sig = area_patient(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_patient(k) = val;      
    end
    
    m_val_control = zeros(nr_controls, 1);
    area = all_data(:,:,ind_control:2:68);
    for k = 1:nr_controls
        area_control = area(:,:,k);
        area_sig = area_control(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_control(k) = val;      
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};
       
        clinical_data = ictal.(clinical_name); % values are the same for all session types of each subject;
        subplot(2,5,d)

        scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
        lsline
        hold on;
        scatter(zeros(nr_controls,1), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
        [r,p] = corr(clinical_data, m_val_patient, 'Type', corr_type);
        text(mean(clinical_data),mean(m_val_patient), ...
            ['r² = ', num2str(round(r^2,3)), '; p = ', num2str(round(p,4))], 'vert', 'middle', 'horiz', 'center', 'FontSize', font)
        ylim([0.15 0.75])
        if s == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
        end
        if c == 1
            ylabel(['\bf \color[rgb]{', num2str(c_patient),'}', session_patient]);
        end

        if s == 2
            xlabel(clinical_units{c})
        end
        
        d = d +1;
    end
end

%% Avg sig FC vs cycle day
corr_type = 'Spearman';
results_path = '/home/iesteves/dFC/NBS1.2/NBS1.2/mydata/results_pca19';

contrast_type = 1;
method = 'Extent';
threshold = 4;

clinical_names = {'cycle_day'};
clinical_units = {'days'};

sessions = {'ictal', 'postictal', 'preictal'};

figure('pos', [50 50 1800 800])
h = suptitle({[corr_type, ' correlation - Avg sig FC edges vs Cycle Day - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;
for s = 1:length(sessions)
    session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        ind_control = 41;
        c_control = c_premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'ic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        ind_control = 41;
        c_control = c_premenstrual;
        comparison = 'postic-prem'; % 'postic-prem';
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        ind_control = 42;
        c_control = c_midcycle;
    end
    
    if s~=3
        load([results_path, '/nbs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
        adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
    else
        continue
    end
    
%     figure('pos', [50, 50, 300 450]);
%     subplot(3,1,[1 2])
%     imAlpha=ones(size(adj_matrix));
%     imAlpha(adj_matrix==0)=0;
%     imagesc(tril(adj_matrix,0),'AlphaData', imAlpha)
%     hold on;
%     plot_atlas_labels(atlas, nr_areas);
%     cmap = gray(256);
%     colormap(flipud(cmap))
%     set(gca, 'XTickLabel', '')
%     set(gca, 'YTickLabel', '')
%     title(comparison)

    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    m_val_patient = zeros(nr_patients, 1);
    area = all_data(:,:,ind_patient:4:40);
    for k = 1:nr_patients
        area_patient = area(:,:,k);
        area_sig = area_patient(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_patient(k) = val;      
    end
    
    m_val_control = zeros(nr_controls, 1);
    area = all_data(:,:,ind_control:2:68);
    for k = 1:nr_controls
        area_control = area(:,:,k);
        area_sig = area_control(logical(tril(adj_matrix, 0)));
        val = mean(area_sig);
        m_val_control(k) = val;      
    end

    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};
        if strcmp(session_patient, 'preictal')
            clinical_data = preictal.(clinical_name); 
        elseif strcmp(session_patient, 'ictal')
            clinical_data = ictal.(clinical_name);
        elseif strcmp(session_patient, 'postictal')
            clinical_data = postictal.(clinical_name);
        end
        clinical_data_controls = premenstrual.(clinical_name);
        %subplot(1,2,d)

        scatter(clinical_data, m_val_patient, 70, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
        %lsline
        hold on;
        if s == 1
            scatter(clinical_data_controls([1,2,4:end]), m_val_control, 70, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
        end
        ylim([0.15 0.75])
        if s == 1
            ht = title(clinical_name);
            ht.Interpreter = 'none';
        end
        if c == 1
            ylabel('Avg FC');
        end

        if s == 2
            xlabel(clinical_units{c})
        end
        
        d = d +1;
    end
end


%% Avg FC vs cycle day
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

corr_type = 'Spearman';

nr_sessions = 4;

clinical_names = {'cycle_day'};
clinical_units = {'days'};

sessions = {'preictal', 'ictal', 'postictal', 'interictal'};

% stats
hvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));
pvec = zeros(nr_sessions, nr_networks*nr_networks, length(clinical_names));

figure('pos', [50 50 1800 1000])
h = suptitle({[corr_type, ' correlation - Avg network FC vs Cycle Day - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
%h = suptitle({[corr_type, ' correlation - Avg network FC vs General Clinical - Thresh ', num2str(thr), '% - ', num2str(nComp), ' PCs'] , '  ', '  '});
h.Interpreter = 'none';
d = 1;                 
for s = 1:nr_sessions
session_patient = sessions{s};

    if strcmp(session_patient, 'preictal')
        ind_patient = 1;
        c_patient = c_preictal;
        data_patient = preictal;
        ind_control = 41;
        c_control = c_premenstrual;
        data_control = premenstrual;
    elseif strcmp(session_patient, 'ictal')
        ind_patient = 2;
        c_patient = c_ictal;
        data_patient = ictal;
        ind_control = 41;
        c_control = c_premenstrual;
        data_control = premenstrual;
    elseif strcmp(session_patient, 'postictal')
        ind_patient = 3;
        c_patient = c_postictal;
        data_patient = postictal;
        ind_control = 41;
        c_control = c_premenstrual;
        data_control = premenstrual;
    elseif strcmp(session_patient, 'interictal')
        ind_patient = 4;
        c_patient = c_interictal;
        data_patient = interictal;
        ind_control = 42;
        c_control = c_midcycle;
        data_control = midcycle;
    end


    all_data = conversion_vecmat(group_data_PCArecon, nr_areas, 'vec2mat');
    % patients
    m_val_patient_all = zeros(nr_networks, nr_networks, nr_patients);
    for a = 1:nr_networks
        for b = 1:nr_networks                
            area = all_data(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_patient:4:40);              
            for k = 1:nr_patients
                area_patient = area(:,:,k);
                if a == b
                    aux = tril(ones(length(area_patient)),-1);
                    area_patient(aux == 0) = NaN;
                end
                m_val_patient_all(a, b, k) = nanmean(area_patient(:));                    
            end
        end
    end
   
    % controls
    m_val_control_all = zeros(nr_networks, nr_networks, nr_patients);
    for a = 1:nr_networks
        for b = 1:nr_networks                
            area = all_data(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_control:2:68);              
            for k = 1:nr_controls
                area_patient = area(:,:,k);
                if a == b
                    aux = tril(ones(length(area_patient)),-1);
                    area_patient(aux == 0) = NaN;
                end
                m_val_control_all(a, b, k) = nanmean(area_patient(:));                    
            end
        end
    end


    for c = 1:length(clinical_names)

        clinical_name = clinical_names{c};

        clinical_data = data_patient.(clinical_name); % values are the same for all session types of each subject;
        clinical_data_controls = data_control.(clinical_name); % values are the same for all session types of each subject;
        d = 1;
        for a = 1:nr_networks
            for b = 1:nr_networks
                subplot(9,9,d)
                scatter(clinical_data, squeeze(m_val_patient_all(a, b,:)), 5, 'o', 'MarkerFaceColor', c_patient, 'MarkerEdgeColor', 'none')
                %lsline
                hold on;
                scatter(clinical_data_controls([1,2,4:end]), squeeze(m_val_control_all(a, b, :)), 5, 'o', 'MarkerFaceColor', c_control, 'MarkerEdgeColor', 'none')
                xlim([-10 37])
                ylim([-0.1 0.9])
                d = d + 1;
                
                if a~= nr_networks
                    set(gca, 'XTickLabel', [])
                end
                if b >1 
                    set(gca, 'YTickLabel', [])
                end
                if a == 1 
                    title(['\bf \color[rgb]{', num2str(colors(b,:)),'}', networks{b}]);
                end
                if b == 1
                    ylabel(['\bf \color[rgb]{', num2str(colors(a,:)),'}', networks{a}]);
                end
            end
        end

    end % clinical variables
     
    %print([fig_path, '/general-clinical_ses-', session_patient,'_network-local_',label_aux{a},'_thresh-', num2str(thr), '_PCArecon-', num2str(nComp), '_', atlas,'-', pipeline], '-dpng')
end % sessions
