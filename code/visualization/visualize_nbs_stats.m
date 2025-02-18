nbs_stats_path = '/home/iesteves/FC/data/results/nbs/stats';
fig_path = '/home/iesteves/FC/figures/nbs_stats';

managefolders(fig_path, 'create');

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB'};

atlas = 'SchaeferSubCRB7100';
nr_areas = 130;

%% Figure - adjacency matrices for each comparison as a function of method and threshold
% only saves the figure if there is any significant edge
comparisons = {'preic-ic', 'preic-postic', 'preic-interic', 'postic-ic', 'ic-interic', 'postic-interic', 'prem-mid', 'preic-prem', 'ic-prem', 'postic-prem', 'interic-mid'};
methods = {'Extent', 'Intensity'};
contrasts = {'1', '-1'};
thresholds = [2, 3.1, 4, 5, 6];
components = [19, 68];

for comp = 1:length(components)
    nComp = components(comp);
 
    for c = 1:length(comparisons)
        comparison = comparisons{c};

        for a = 1:length(contrasts)
            contrast_type = contrasts{a};   

            d = 0;
            figure('pos', [50 50 350 1100]);
            suptitle({[num2str(nComp), 'PCs - ',comparison], ' ', ' '})
            for t = 1:length(thresholds)
                threshold = thresholds(t);
                for m = 1:length(methods)
                    method= methods{m};

                    d = d +1;
                    try

                    load([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(thresholds(t)))]);
                    adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
                    subplot(5,2,d)
                    imagesc(adj_matrix)
                    hold on;
                    plot_atlas_labels(atlas, nr_areas);
                    set(gca, 'XTickLabel', '')
                    set(gca, 'YTickLabel', '')
                    title(method)
                    ylabel(['thr = ',num2str(threshold)])

                    print([fig_path, '/nbs_stats_thresholds_', num2str(nComp),'PCs_', comparison], '-dpng')
                    catch
                        continue
                    end
                end
            end
        end
    end
end



