DI_path = '/home/iesteves/FC/data/results/DI';
fig_path = '/home/iesteves/FC/figures/DI_FC_Imatrix_PCA';

managefolders(fig_path, 'create');

%load([DI_path, '/DI_FC_PCA/FCvec_PCA_recon.mat']);
%load([DI_path, '/DI_FC_PCA/FC_PCA_recon.mat']);
load([DI_path, '/DI_FC_PCA/corr_PCA_recon_networks.mat']);

pipeline = 'preprocessed_rp_mo_csf_wm_nui';
atlas = 'SchaeferSubCRB7100';

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB', 'WB'};

% number of networks
nr_networks = length(networks);
nr_areas = 130;

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
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0; 0 0 0];

components = [2, 8, 12, 19, 25, 35, 40, 68];

%%
min_val = 0;
max_val = 1;
c = 1;
% FIGURE - corr matrix for different nb of PCs and networks
figure('pos', [100, 200, 1600 1200]);
ht = suptitle({['Identifiability matrices for PCA reconstructed data (networks + whole brain)'], [ atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
  
    for n = 1:nr_networks
        ax = subplot(length(components), nr_networks, c);
        imagesc(corr_PCA_recon_networks(:,:,n,nComp))
        if comp == 1
            title(['\color[rgb]{', num2str(colors_networks(n,:)),'}', networks{n}])
        end
        if n == 1
            ylabel(['\bf', num2str(nComp), ' PCs'])
            
        end
        c = c+1;
        if comp == length(components)
            h = xlabel('[Subjects x Sessions]');           
            set(h, 'FontSize', 8)
        end
        set(gca, 'XTickLabel', [])
        set(gca, 'YTickLabel', [])
        caxis([min_val max_val])
        cb = colorbar(ax); 
        cb.Position = [.92 .4 .02 .3];
       
    end   
end

print([fig_path, '/DI_FC-Imatrix_PCA_networks-WB'], '-dpng')

%% WB
min_val = 0;
max_val = 1;
c = 1;

n = 10; 

% FIGURE - corr matrix for different nb of PCs and networks
figure('pos', [100, 200, 1900 800]);
ht = suptitle({['Identifiability matrices for PCA reconstructed data (whole brain)'], [ atlas, ', ', num2str(nr_areas), ' ROIs - ', pipeline,...
        ' - [', num2str(min_val), ' ', num2str(max_val), ']'], ' ', ' '});
ht.Interpreter = 'none';
for comp = 1:length(components)
    nComp = components(comp);
  
    ax = subplot(2,4,comp);
    imagesc(corr_PCA_recon_networks(:,:,n,nComp))

    title(['\bf', num2str(nComp), ' PCs'])
    h = xlabel('[Subjects x Sessions]');           
    g = ylabel('[Subjects x Sessions]');       
    set(gca, 'XTickLabel', [])
    set(gca, 'YTickLabel', [])
    caxis([min_val max_val])
    cb = colorbar(ax); 
    cb.Position = [.92 .4 .02 .3];
    set(gca, 'FontSize', 10)
    set(gca, 'FontSize', 14)

end

print([fig_path, '/DI_FC-Imatrix_PCA_WB'], '-dpng')