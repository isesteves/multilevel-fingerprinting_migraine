addpath('/home/iesteves/stats/MultipleTestingToolbox/MultipleTestingToolbox');

DI_path = '/home/iesteves/FC/data/results/DI';
fig_path = '/home/iesteves/FC/figures/DI_FC_PCA';

managefolders(fig_path, 'create');

load([DI_path, '/DI_FC_PCA/FCvec_PCA_recon.mat']);
load([DI_path, '/DI_FC_PCA/FC_PCA_recon.mat']);
load([DI_path, '/DI_FC_PCA/corr_PCA_recon_networks.mat']);

pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB', 'WB'};

% number of networks
nr_networks = length(networks);

nr_areas = size(FC_PCA_recon, 1);
nr_subjects = size(FC_PCA_recon, 3);

% colormap
cmap = redbluecmap;
newCmap = imresize(cmap, [64, 3]);  % original color map contain just 11 colors, this increase it to 64
newCmap = min(max(newCmap, 0), 1);

% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;

colors_networks = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];

components = [2, 8, 12, 19, 25, 35, 40, 68];


%% Patient data - individual
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};
min_val = -1;
max_val = 1;

for comp = 1:length(components)
    nComp = components(comp);
 
    c = 1;
    figure('pos', [50 50 1900 700])
    ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Patients (complete, N=10) - ', ...
        atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
    ht.Interpreter = 'none';
    for k = 1:4:40
        subplot(4,10,c)
        imagesc(FC_PCA_recon(:,:,k,nComp), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_preictal),'} preic'])
        end
        sub = subjects{c};         
        title(sub(end-2:end))
        colormap(newCmap)


        subplot(4,10,c+10)
        imagesc(FC_PCA_recon(:,:,k+1,nComp), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_ictal),'} ic'])
        end
        colormap(newCmap)

        subplot(4,10,c+20)
        imagesc(FC_PCA_recon(:,:,k+2,nComp), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_postictal),'} postic'])
        end
        colormap(newCmap)


        ax = subplot(4,10,c+30);
        imagesc(FC_PCA_recon(:,:,k+3,nComp), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_interictal),'} interic'])
        end
        colormap(newCmap)
        caxis([min_val max_val])
        cb = colorbar(ax); 
        cb.Position = [.92 .4 .01 .3];


        c= c+1;
    end
      print([fig_path, '/DI_FC_PCA_FCmatrix_individual-patient_', num2str(nComp), 'PCs'], '-dpng')
end

%% Patient data - session average

min_val = -1;
max_val = 1;

components = [2, 8, 12, 19, 25, 35, 40, 68];

c = 1;
figure('pos', [50 50 1300 700])
ht = suptitle({['AVERAGE - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
   
    subplot(4,length(components),c)
    imagesc(mean(FC_PCA_recon(:,:,1:4:40,nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_preictal),'} preic'])
    end
    title([num2str(nComp), ' PCs'])

    subplot(4,length(components),c+length(components))
    imagesc(mean(FC_PCA_recon(:,:,2:4:40,nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_ictal),'} ic'])
    end

    subplot(4,length(components),c+length(components)*2)
    imagesc(mean(FC_PCA_recon(:,:,3:4:40,nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_postictal),'} postic'])
    end

    ax1 = subplot(4,length(components),c+length(components)*3);
    imagesc(mean(FC_PCA_recon(:,:,4:4:40,nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_interictal),'} interic'])
    end
    colormap(newCmap)
    caxis([min_val max_val])
    cb = colorbar(ax1); 
    cb.Position = [.92 .4 .01 .3];
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_FCmatrix_average-patient'], '-dpng')
    
%% Patient data - session std

min_val = -0.3;
max_val = 0.3;

c = 1;
figure('pos', [50 50 1300 700])
ht = suptitle({['STD - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
   
    subplot(4,length(components),c)
    imagesc(std(FC_PCA_recon(:,:,1:4:40,nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_preictal),'} preic'])
    end
    title([num2str(nComp), ' PCs'])

    subplot(4,length(components),c+length(components))
    imagesc(std(FC_PCA_recon(:,:,2:4:40,nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_ictal),'} ic'])
    end

    subplot(4,length(components),c+length(components)*2)
    imagesc(std(FC_PCA_recon(:,:,3:4:40,nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_postictal),'} postic'])
    end

    ax1 = subplot(4,length(components),c+length(components)*3);
    imagesc(std(FC_PCA_recon(:,:,4:4:40,nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_interictal),'} interic'])
    end
    colormap(newCmap)
    caxis([min_val max_val])
    cb = colorbar(ax1); 
    cb.Position = [.92 .4 .01 .3];
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_FCmatrix_std-patient'], '-dpng')
 
%% Control data
subjects = {'sub-control019','sub-control020','sub-control026','sub-control027','sub-control028','sub-control029',...
    'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};

min_val = -1;
max_val = 1;


for comp = 1:length(components)
    nComp = components(comp);
   
    c = 1;
    figure('pos', [50 50 1600 600])
    ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Controls (complete, N=14) - ', ...
        atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
    ht.Interpreter = 'none';
    for k = 41:2:68
        subplot(5,14,c)
        imagesc(FC_PCA_recon(:,:,k), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==41
            ylabel(['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens'])
        end
        sub = subjects{c};         
        title(sub(end-2:end))

        ax2 = subplot(5,14,c+14);
        imagesc(FC_PCA_recon(:,:,k+1), [min_val max_val])
        hold on;
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 13)
        if k==41
            ylabel(['\bf \color[rgb]{', num2str(c_midcycle),'} postov'])
        end
        colormap(newCmap)
        caxis([min_val max_val])
        cb = colorbar(ax2); 
        cb.Position = [.92 .63 .01 .3];
        c= c+1;
    end
    print([fig_path, '/DI_FC_PCA_FCmatrix_individual-control_', num2str(nComp), 'PCs_FCmatrix'], '-dpng')
end

%% Controls - session average

min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1700 700])
ht = suptitle({['AVERAGE - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    
    subplot(3,length(components),c)
    imagesc(mean(FC_PCA_recon(:,:,41:2:68, nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens'])
    end
    title([num2str(nComp), ' PCs'])

    ax3 = subplot(3,length(components),c+length(components));
    imagesc(mean(FC_PCA_recon(:,:,42:2:68, nComp),3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_midcycle),'} postov'])
    end
    colormap(newCmap)
    caxis([min_val max_val])
    cb = colorbar(ax3); 
    cb.Position = [.92 .4 .01 .3]; 
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_FCmatrix_avg-control'], '-dpng')

%% Controls - session std

min_val = -0.3;
max_val = 0.3;

c = 1;
figure('pos', [50 50 1700 700])
ht = suptitle({['STD - Pearson correlation - Patients (complete, N=10) - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    
    subplot(3,length(components),c)
    imagesc(std(FC_PCA_recon(:,:,41:2:68, nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens'])
    end
    title([num2str(nComp), ' PCs'])

    ax3 = subplot(3,length(components),c+length(components));
    imagesc(std(FC_PCA_recon(:,:,42:2:68, nComp),0,3), [min_val max_val])
    hold on;
    plot_atlas_labels(atlas, nr_areas);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'FontSize', 14)
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_midcycle),'} postov'])
    end
    colormap(newCmap)
    caxis([min_val max_val])
    cb = colorbar(ax3); 
    cb.Position = [.92 .4 .01 .3]; 
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_FCmatrix_std-control'], '-dpng')

%% Histograms
min_val = -1;
max_val = 1;

c = 1;
figure('pos', [50 50 1300 700])
ht = suptitle({['HISTOGRAM - Pearson correlation - All subjects and edges - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
   
    subplot(6,length(components),c)
    histogram(FCvec_PCA_recon(:,41:2:68,nComp), 50, 'facecolor', c_premenstrual, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens'])
    end
    title([num2str(nComp), ' PCs'])
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    
    subplot(6,length(components),c+length(components))
    histogram(FCvec_PCA_recon(:,42:2:68,nComp), 50, 'facecolor', c_midcycle, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_midcycle),'} postov'])
    end
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    
    subplot(6,length(components),c+length(components)*2)
    histogram(FCvec_PCA_recon(:,1:4:40,nComp), 50, 'facecolor', c_preictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_preictal),'} preic'])
    end
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)

    subplot(6,length(components),c+length(components)*3)
    histogram(FCvec_PCA_recon(:,2:4:40,nComp), 50, 'facecolor', c_ictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_ictal),'} ic'])
    end
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)

    subplot(6,length(components),c+length(components)*4)
    histogram(FCvec_PCA_recon(:,3:4:40,nComp), 50, 'facecolor', c_postictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_postictal),'} postic'])
    end
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)

    subplot(6,length(components),c+length(components)*5);
    histogram(FCvec_PCA_recon(:,4:4:40,nComp), 50, 'facecolor', c_interictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'YTickLabel', []);
    if c == 1
       ylabel(['\bf \color[rgb]{', num2str(c_interictal),'} interic'])
    end
    xlim([min_val max_val])
    ylim([0 10000])  
    hold on;
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_histograms_FCmatrix_edges-all_subjects-all'], '-dpng')


%% average across subjects for each connection; sessions overlapping per group
min_val = -0.2;
max_val = 1;
nr_bins = 100;

c = 1;
figure('pos', [50 50 900 900])
ht = suptitle({['HISTOGRAM - Pearson correlation - Edges (subjects avg) - Per group'], ...
    [atlas, ', ', num2str(nr_areas),' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
   
    subplot(length(components),2,c)
    histogram(mean(FCvec_PCA_recon(:,41:2:68,nComp),2), nr_bins, 'facecolor', c_premenstrual, 'facealpha',.4, 'edgealpha', 0)
    hold on;
    histogram(mean(FCvec_PCA_recon(:,42:2:68,nComp),2), nr_bins, 'facecolor', c_midcycle, 'facealpha',.4, 'edgealpha', 0)
    if comp < length(components)
       set(gca, 'XTickLabel', []); 
    end
    set(gca, 'YTickLabel', []);
    if c == 1
       title('HC')
    end
    ylabel(['\bf \color{black}', num2str(nComp), ' PCs'])
    xlim([min_val max_val])
    ylim([0 300])  
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,41:2:68,nComp),2)) mean(mean(FCvec_PCA_recon(:,41:2:68,nComp),2))], ylim, 'Color', c_premenstrual, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,42:2:68,nComp),2)) mean(mean(FCvec_PCA_recon(:,42:2:68,nComp),2))], ylim, 'Color', c_midcycle, 'LineStyle', '-.', 'LineWidth', 1.1)
    
    c = c+1;
    subplot(length(components),2, c)
    histogram(mean(FCvec_PCA_recon(:,1:4:40,nComp),2), nr_bins, 'facecolor', c_preictal, 'facealpha',.4, 'edgealpha', 0)
    hold on;
    histogram(mean(FCvec_PCA_recon(:,2:4:40,nComp),2), nr_bins, 'facecolor', c_ictal, 'facealpha',.4, 'edgealpha', 0)
    histogram(mean(FCvec_PCA_recon(:,3:4:40,nComp),2), nr_bins, 'facecolor', c_postictal, 'facealpha',.4, 'edgealpha', 0)
    histogram(mean(FCvec_PCA_recon(:,4:4:40,nComp),2), nr_bins, 'facecolor', c_interictal, 'facealpha',.4, 'edgealpha', 0)
    if comp < length(components)
       set(gca, 'XTickLabel', []); 
    end
    set(gca, 'YTickLabel', []);
    if c == 2
       title('MIG')
    end
    xlim([min_val max_val])
    ylim([0 300])  
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,1:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,1:4:40,nComp),2))], ylim, 'Color', c_preictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,2:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,2:4:40,nComp),2))], ylim, 'Color', c_ictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,3:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,3:4:40,nComp),2))], ylim, 'Color', c_postictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,4:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,4:4:40,nComp),2))], ylim, 'Color', c_interictal, 'LineStyle', '-.', 'LineWidth', 1.1)
 
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_histograms_FCmatrix_edges-all_subjects-avg_pergroup'], '-dpng')


%% average across subjects for each connection; sessions overlapping per comparison
min_val = -0.2;
max_val = 1;
nr_bins = 100;

c = 1;
figure('pos', [50 50 900 900])
ht = suptitle({['HISTOGRAM - Pearson correlation - Edges (subjects avg) - per cycle'], ...
    [atlas, ', ', num2str(nr_areas),' ROIs - ', pipeline,...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
   
    subplot(length(components),2,c)
    histogram(mean(FCvec_PCA_recon(:,41:2:68,nComp),2), nr_bins, 'facecolor', c_premenstrual, 'facealpha',.4, 'edgealpha', 0)
    hold on;
    histogram(mean(FCvec_PCA_recon(:,1:4:40,nComp),2), nr_bins, 'facecolor', c_preictal, 'facealpha',.4, 'edgealpha', 0)
    hold on;
    histogram(mean(FCvec_PCA_recon(:,2:4:40,nComp),2), nr_bins, 'facecolor', c_ictal, 'facealpha',.4, 'edgealpha', 0)
    histogram(mean(FCvec_PCA_recon(:,3:4:40,nComp),2), nr_bins, 'facecolor', c_postictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'YTickLabel', []);
    if c == 1
       title('perimenstrual')
    end
    if comp < length(components)
       set(gca, 'XTickLabel', []); 
    end
    ylabel(['\bf \color{black}', num2str(nComp), ' PCs'])
    xlim([min_val max_val])
    ylim([0 300])  
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,41:2:68,nComp),2)) mean(mean(FCvec_PCA_recon(:,41:2:68,nComp),2))], ylim, 'Color', c_premenstrual, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,1:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,1:4:40,nComp),2))], ylim, 'Color', c_preictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,2:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,2:4:40,nComp),2))], ylim, 'Color', c_ictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,3:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,3:4:40,nComp),2))], ylim, 'Color', c_postictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    
    c = c+1;
    subplot(length(components),2, c)
    histogram(mean(FCvec_PCA_recon(:,42:2:68,nComp),2), nr_bins, 'facecolor', c_midcycle, 'facealpha',.4, 'edgealpha', 0)
    hold on;
    histogram(mean(FCvec_PCA_recon(:,4:4:40,nComp),2), nr_bins, 'facecolor', c_interictal, 'facealpha',.4, 'edgealpha', 0)
    set(gca, 'YTickLabel', []);
    if c == 2
       title('postovulation')
    end
    xlim([min_val max_val])
    ylim([0 300])  
    line([0 0], ylim, 'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,42:2:68,nComp),2)) mean(mean(FCvec_PCA_recon(:,42:2:68,nComp),2))], ylim, 'Color', c_midcycle, 'LineStyle', '-.', 'LineWidth', 1.1)
    line([mean(mean(FCvec_PCA_recon(:,4:4:40,nComp),2)) mean(mean(FCvec_PCA_recon(:,4:4:40,nComp),2))], ylim, 'Color', c_interictal, 'LineStyle', '-.', 'LineWidth', 1.1)
    if comp < length(components)
       set(gca, 'XTickLabel', []); 
    end
    
    c = c+1;
end

print([fig_path, '/DI_FC_PCA_histograms_FCmatrix_edges-all_subjects-avg_percycle'], '-dpng')

%% Statistics with univariate tests - Patients and Controls
session1 = [1, 1, 1, 3, 2, 3, 41];
session2 = [2, 3, 4, 2, 4, 4, 42];
comparisons = {'preic-ic', 'preic-postic', 'preic-interic', 'postic-ic', 'ic-interic', 'postic-interic', 'prem-mid'};
p_vec = zeros((nr_areas-1)*nr_areas/2, 6, length(components));
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        for edge = 1:(nr_areas-1)*nr_areas/2          
            if c==7
               p = signrank(FCvec_PCA_recon(edge,session1(c):2:68, comp),FCvec_PCA_recon(edge,session2(c):2:68, comp));
            else
               p = signrank(FCvec_PCA_recon(edge,session1(c):4:40, comp),FCvec_PCA_recon(edge,session2(c):4:40, comp));
            end
            p_vec(edge,c, comp) = p;
        end 
    end
end

min_val = 0.95;
max_val = 1;
d = 1;
figure('pos', [50 50 1400 1300])
ht = suptitle({['Edgewise Wilcoxon Signed Rank test (1-p-value, uncorrected) - Patients and Controls - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs',...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        ax4 = subplot(length(components),length(comparisons),d);
        imagesc(conversion_vecmat(1-p_vec(:,c, comp), nr_areas, 'vec2mat'), [min_val max_val])
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        title(comparisons{c})
        caxis([min_val max_val])
        cb = colorbar(ax4); 
        cb.Position = [.92 .4 .01 .3]; 
    
        if any(d==[1,8,15,22,29,36,43])
            ylabel([num2str(nComp), ' PCs'])
        end
        d = d+1;
    end
end

print([fig_path, '/DI_FC_PCA_stats-edgewise_uncorrected_FCmatrix_pergroup'], '-dpng')


min_val = 0.95;
max_val = 1;
d = 1;
figure('pos', [50 50 1400 1300])
ht = suptitle({['Edgewise Wilcoxon Signed Rank test (1-p-value, FDR corrected) - Patients and Controls - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs',...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        ax5 = subplot(length(components),length(comparisons),d);
        [c_pvalues, c_alpha, h, extra] = fdr_BH(p_vec(:,c, comp), 0.05, false);
        imagesc(conversion_vecmat(1-c_pvalues', nr_areas, 'vec2mat'), [min_val max_val])
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        title(comparisons{c})
        caxis([min_val max_val])
        cb = colorbar(ax5); 
        cb.Position = [.92 .4 .01 .3]; 
       
        if any(d==[1,8,15,22,29,36,43])
            ylabel([num2str(nComp), ' PCs'])
        end
        d = d+1;
    end
end

print([fig_path, '/DI_FC_PCA_stats-edgewise_FDRcorrected-permatrix_FCmatrix_pergroup'], '-dpng')

%% Statistics with univariate tests - Patients vs Controls
session1 = [1, 2, 3, 4];
session2 = [41, 41, 41, 42];
comparisons = {'preic-prem', 'ic-prem', 'postic-prem', 'interic-mid'};
p_vec = zeros((nr_areas-1)*nr_areas/2, 4, length(components));
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        for edge = 1:(nr_areas-1)*nr_areas/2

            p = ranksum(FCvec_PCA_recon(edge,session1(c):4:40, comp),FCvec_PCA_recon(edge,session2(c):2:68, comp));
            p_vec(edge,c, comp) = p;
        end 
    end
end

min_val = 0.95;
max_val = 1;
d = 1;
figure('pos', [50 50 1000 1300])
ht = suptitle({['Edgewise Wilcoxon Rank Sum test (1-p-value, uncorrected) - Patients vs Controls - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs',...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        ax6 = subplot(length(components),length(comparisons),d);
        imagesc(conversion_vecmat(1-p_vec(:,c, comp), nr_areas, 'vec2mat'), [0.95 1])
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        title(comparisons{c})
        caxis([min_val max_val])
        cb = colorbar(ax6); 
        cb.Position = [.92 .4 .01 .3]; 

        
        if any(d==[1,5,9,13,17,21])
            ylabel([num2str(nComp), ' PCs'])            
        end

        d = d+1;
    end
end

print([fig_path, '/DI_FC_PCA_stats-edgewise_uncorrected_FCmatrix_percycle'], '-dpng')


min_val = 0.95;
max_val = 1;
d = 1;
figure('pos', [50 50 1000 1300])
ht = suptitle({['Edgewise Wilcoxon Rank Sum test (1-p-value, FDR correction) - Patients vs Controls - ', ...
    atlas, ', ', num2str(nr_areas), ' ROIs',...
    ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    for c = 1:length(comparisons)
        ax6 = subplot(length(components),length(comparisons),d);
        [c_pvalues, c_alpha, h, extra] = fdr_BH(p_vec(:,c, comp), 0.05, false);
        imagesc(conversion_vecmat(1-c_pvalues', nr_areas, 'vec2mat'), [0.95 1])
        plot_atlas_labels(atlas, nr_areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        title(comparisons{c})
        caxis([min_val max_val])
        cb = colorbar(ax6); 
        cb.Position = [.92 .4 .01 .3]; 

        
        if any(d==[1,5,9,13,17,21])
            ylabel([num2str(nComp), ' PCs'])            
        end

        d = d+1;
    end
end

print([fig_path, '/DI_FC_PCA_stats-edgewise_FDRcorrected-permatrix_FCmatrix_percycle'], '-dpng')

%% Average FC per network
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 130 9];

average_FC = zeros(nr_networks-1, nr_networks-1, nr_subjects, length(components));
average_FC_WB = zeros(nr_subjects, length(components));
for comp = 1:length(components)
    nComp = components(comp);
    for k = 1:nr_subjects
        subject_FC = FC_PCA_recon(:,:, k, nComp);
        for a = 1:nr_networks-1
            for b = 1:nr_networks-1                
                area = subject_FC(areas(a,1):areas(a,2), areas(b,1):areas(b,2));
                average_FC(a, b, k, comp) = mean(area(:));
            end
        end
        average_FC_WB(k, comp) = mean(subject_FC(:));
    end
end

%% Patient data
subjects = {'sub-patient002', 'sub-patient003', 'sub-patient005', 'sub-patient006', 'sub-patient008',...
    'sub-patient009', 'sub-patient034', 'sub-patient038', 'sub-patient041', 'sub-patient045'};
min_val = -1;
max_val = 1;

for comp = 1:length(components)
    nComp = components(comp);
    data = average_FC(:,:,:,comp);

    c = 1;
    figure('pos', [50 50 1700 700])
    ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Patients (complete, N=10) - ', ...
        atlas, ', ', ' 9 networks (Avg)- ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
    ht.Interpreter = 'none';
    for k = 1:4:40
        subplot(4,10,c)
        imagesc(data(:,:,k), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_preictal),'} preic'])
        end
        sub = subjects{c};         
        title(sub(end-2:end))

        subplot(4,10,c+10)
        imagesc(data(:,:,k+1), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_ictal),'} ic'])
        end

        subplot(4,10,c+20)
        imagesc(data(:,:,k+2), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_postictal),'} postic'])
        end

        ax7 = subplot(4,10,c+30);
        imagesc(data(:,:,k+3), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 14)
        if k==1
            ylabel(['\bf \color[rgb]{', num2str(c_interictal),'} interic'])
        end
        colormap(newCmap)
        caxis([min_val max_val])
        cb = colorbar(ax7); 
        cb.Position = [.92 .4 .01 .3]; 

        c= c+1;
    end
    print([fig_path, '/DI_FC_PCA_FCnetworks_individual-patient_', num2str(nComp), 'PCs'], '-dpng')
end

%% Control data
subjects = {'sub-control019','sub-control020','sub-control026','sub-control027','sub-control028','sub-control029',...
    'sub-control030', 'sub-control031', 'sub-control033', 'sub-control044', 'sub-control046', 'sub-control048', 'sub-control049', 'sub-control051'};
min_val = -1;
max_val = 1;

for comp = 1:length(components)
    nComp = components(comp);
    data = average_FC(:,:,:,comp);

    c = 1;
    figure('pos', [50 50 1700 600])
    ht = suptitle({[num2str(nComp), ' PCs - Pearson correlation - Controls (complete, N=14) - ', ...
        atlas, ', ', ' 9 networks (Avg)- ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
    ht.Interpreter = 'none';
    for k = 41:2:68
        subplot(5,14,c)
        imagesc(data(:,:,k), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 12)
        if k==41
            ylabel(['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens'])
        end
        sub = subjects{c};         
        title(sub(end-2:end))

        ax8 = subplot(5,14,c+14);
        imagesc(data(:,:,k+1), [min_val max_val])
        hold on;
        plot_networks_labels(areas);
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'FontSize', 12)
        if k==41
            ylabel(['\bf \color[rgb]{', num2str(c_midcycle),'} postov'])
        end
        colormap(newCmap)
        caxis([min_val max_val])
        cb = colorbar(ax8); 
        cb.Position = [.92 .6 .01 .3]; 
        
        c= c+1;
    end
    print([fig_path, '/DI_FC_PCA_FCnetworks_individual-control_', num2str(nComp), 'PCs'], '-dpng')
end

%% Average WB
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

d = 1;
figure('pos', [50, 50, 1200 800]);
ht = suptitle({['Whole-brain Avg FC - Pearson correlation - ', ...
        atlas, ' - ', pipeline], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    subplot(2,4,d)
    b= boxplot([average_FC_WB(41:2:68,comp) average_FC_WB(42:2:68,comp)...
        [average_FC_WB(1:4:40,comp); NaN*ones(4,1)]...
        [average_FC_WB(2:4:40,comp); NaN*ones(4,1)]...
        [average_FC_WB(3:4:40,comp); NaN*ones(4,1)]...
        [average_FC_WB(4:4:40,comp); NaN*ones(4,1)]...
        ], 'Color', 'k');
    set(b, 'linew', 1.1)
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
    end
    plot_bp_pnts([average_FC_WB(41:2:68,comp) average_FC_WB(42:2:68,comp) ...
        [average_FC_WB(1:4:40,comp); NaN*ones(4,1)]...
        [average_FC_WB(2:4:40,comp); NaN*ones(4,1)]...
        [average_FC_WB(3:4:40,comp); NaN*ones(4,1)] ...
        [average_FC_WB(4:4:40,comp); NaN*ones(4,1)]], 10)
    
    ylim([0.1 0.85])
    line([2.5 2.5], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1)
    title([num2str(nComp), ' PCs'])
    set(gca, 'XTickLabel', [])
    set(gca, 'FontSize', 12)
    grid minor
    set(gca, 'XTickLabel', {'perimens', 'postov', 'preic', 'ic', 'postic', 'interic'})
    set(gca, 'XTickLabelRotation', 30)
    hAx=gca;
    hAx.XAxis.TickLabelInterpreter='tex';       
    hAx.XTickLabel(1)={['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens']};
    hAx.XTickLabel(2)={['\bf \color[rgb]{', num2str(c_midcycle),'} postov']};
    hAx.XTickLabel(3)={['\bf \color[rgb]{', num2str(c_preictal),'} preic']};
    hAx.XTickLabel(4)={['\bf \color[rgb]{', num2str(c_ictal),'} ic']};
    hAx.XTickLabel(5)={['\bf \color[rgb]{', num2str(c_postictal),'} postic']};
    hAx.XTickLabel(6)={['\bf \color[rgb]{', num2str(c_interictal),'} interic']};
    if any(d ==[1,5])
        ylabel('Average FC')
    end
    
    d= d+1;
    print([fig_path, '/DI_FC_PCA_boxplot_FCavg-WB'], '-dpng')
end

%% Average networks
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

d = 1;
figure('pos', [50 50 2000 1200]);
ht = suptitle({['Network Avg FC - Pearson correlation - ', ...
        atlas, ' - ', pipeline], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
    
    for n = 1:nr_networks-1
        subplot(length(components),nr_networks-1,d)
        boxplot([squeeze(average_FC(n,n,41:2:68,comp)) squeeze(average_FC(n,n,42:2:68,comp)), ...
            [[squeeze(average_FC(n,n,1:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,2:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,3:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,4:4:40,comp)); NaN*ones(4,1)]]...
            ]);
        hold on;
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        plot_bp_pnts([squeeze(average_FC(n,n,41:2:68,comp)) squeeze(average_FC(n,n,42:2:68,comp)), ...
            [[squeeze(average_FC(n,n,1:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,2:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,3:4:40,comp)); NaN*ones(4,1)]...
            [squeeze(average_FC(n,n,4:4:40,comp)); NaN*ones(4,1)]]...
            ], 5)
        if comp == 1
            title(networks{n})
        end
        if n == 1
            ylabel(['\bf \color{black}', num2str(nComp), ' PCs'])
        end
        set(gca, 'XTickLabel', [])
        if n > 1
            set(gca, 'YTickLabel', [])
        end
        set(gca, 'FontSize', 12)
        grid minor
        ylim([0.05 0.85])
        line([2.5 2.5], ylim, 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1)
%         if comp == length(components)
%             set(gca, 'XTickLabel', {'perimens', 'postov', 'preic', 'ic', 'postic', 'interic'})
%             set(gca, 'XTickLabelRotation', 30)
%             hAx=gca;
%             hAx.XAxis.TickLabelInterpreter='tex';       
%             hAx.XTickLabel(1)={['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens']};
%             hAx.XTickLabel(2)={['\bf \color[rgb]{', num2str(c_midcycle),'} postov']};
%             hAx.XTickLabel(3)={['\bf \color[rgb]{', num2str(c_preictal),'} preic']};
%             hAx.XTickLabel(4)={['\bf \color[rgb]{', num2str(c_ictal),'} ic']};
%             hAx.XTickLabel(5)={['\bf \color[rgb]{', num2str(c_postictal),'} postic']};
%             hAx.XTickLabel(6)={['\bf \color[rgb]{', num2str(c_interictal),'} interic']};
%             hAx.XAxis.FontSize = 8;
%         end
        d= d+1;
    end
end

print([fig_path, '/DI_FC_PCA_boxplot_FCavg-networks'], '-dpng')