nbs_stats_path = '/home/iesteves/FC/data/results/nbs/stats';
nbs_analysis_path = '/home/iesteves/FC/data/results/nbs/analysis';
fig_path = '/home/iesteves/FC/figures/nbs_analysis/';

managefolders(fig_path, 'create');

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB'};
atlas = 'SchaeferSubCRB7100';
nr_areas = 130;
nr_networks = 9;

colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0.25 0.25 0.25];


%% Figure - Sig edges matrices + % significant edges

components = [19, 68];
thresholds = [3.1, 4];
edge_plots = {'count_sig_FC', 'perc_sig_FC', 'percarea_sig_FC'}; 
method = 'Extent';

for  t = 1:length(thresholds)
    threshold = thresholds(t);
   
    for comp = 1:length(components)
        nComp = components(comp);

        load([nbs_analysis_path, '/analysis-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_sig-metrics_networks_method-',method, '_cluster-', num2str(round(threshold))])
        comparisons = sigmetrics.comparisons;
        contrast_types = sigmetrics.contrast_types;
        for e = 1:length(edge_plots)
             data_plot = edge_plots{e};
             edge_data = sigmetrics.(data_plot);

            if strcmp(data_plot, 'count_sig_FC')
                plot_name = 'sig-count';  
                edge_data_ylabel = {'# Significant', 'Edges'};    
            elseif strcmp(data_plot, 'perc_sig_FC')
                plot_name = 'sig-perc';
                edge_data_ylabel = {'% Significant', 'Edges'};
            elseif strcmp(data_plot, 'percarea_sig_FC')
                plot_name = 'sig-percarea';
                edge_data_ylabel = {'% Area with ', 'Significant Edges'};
            end
            
            for c = 1:length(comparisons)
                comparison = comparisons{c};
                contrast_type = contrast_types{c};

                try
                    load([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
                catch ME
                    continue
                end
                adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');
                matrix = tril(round(edge_data(:,:,c),1));

                % Sort the elements of the matrix in ascending order
                [sortedValues, linearIndices] = sort(matrix(:), 'descend');

                % Get the row and column indexes corresponding to the sorted values
                [rowIndices, colIndices] = ind2sub(size(matrix), linearIndices);
                indexes = [rowIndices colIndices];

                sig_networks = indexes(1:5, :);
                connection_labels = cell(1,length(sig_networks));
                for s = 1:length(sig_networks)
                    if indexes(s,1)==indexes(s,2)        
                        connection = ['\bf \color[rgb]{', num2str(colors(indexes(s,1),:)),'} ' networks{indexes(s,1)}];
                    else
                        net1 = ['\bf \color[rgb]{', num2str(colors(indexes(s,1),:)),'} ' networks{indexes(s,1)}];
                        net2 = ['\color[rgb]{', num2str(colors(indexes(s,2),:)),'} ' networks{indexes(s,2)}];
                        connection = [net1, '\color{black} -', net2];

                    end
                    connection_labels{s} = connection;
                end

                figure('pos', [50, 50, 300 450]);
                subplot(3,1,[1 2])
                imAlpha=ones(size(adj_matrix));
                imAlpha(adj_matrix==0)=0;
                imagesc(tril(adj_matrix,0),'AlphaData', imAlpha)
                hold on;
                plot_atlas_labels(atlas, nr_areas);
                cmap = gray(256);
                colormap(flipud(cmap))
                set(gca, 'XTickLabel', '')
                set(gca, 'YTickLabel', '')
                title(comparison)

                subplot(3,1,3)
                bar(sortedValues(1:5), 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
                set(gca, 'XTickLabel', connection_labels)
                set(gca, 'FontSize', 11)
                set(gca, 'XTickLabelRotation', 30)
                ylim([0 1.1*sortedValues(1)])
                xlim([0.5 5.5])
                ylabel(edge_data_ylabel)

                print([fig_path, '/nbs_', num2str(nComp),'PCs_FCmatrix-sig_', plot_name, '_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))], '-dpng')
            
            end % comparisons
        end % significant edges metrics
    end % components
end % thresholds

%% Figure - network matrix - #, % and %network of significant edges
components = [19, 68];
thresholds = [3.1, 4];
edge_plots = {'count_sig_FC', 'perc_sig_FC', 'percarea_sig_FC'}; 
method = 'Extent';

for t = 1:length(thresholds)
    threshold = thresholds(t);
   
    for comp = 1:length(components)
        nComp = components(comp);

        for e = 1:length(edge_plots)
             data_plot = edge_plots{e};
             edge_data = sigmetrics.(data_plot);
    
            if strcmp(data_plot, 'count_sig_FC')
                edge_data_title = {'Number of significant edges per network/network-network', ' ', ' '};
                plot_name = 'sig-count'; 
            elseif strcmp(data_plot, 'perc_sig_FC')
                plot_name = 'sig-perc';
                edge_data_title = {'Percentage of the total significant edges belonging to each network/network-network', ' ', ' '};
            elseif strcmp(data_plot, 'percarea_sig_FC')
                plot_name = 'sig-percarea';
                edge_data_title = {'Percentage of edges of each network/network-network considered significantly different', ' ', ' '};
            end

            figure('pos', [100 100 1900 800]);
            suptitle(edge_data_title)
            for c = 1:length(comparisons)
                comparison = comparisons{c};

                matrix = tril(round(edge_data(:,:,c),1));

                subplot(2, 4, c)
                imagesc(matrix)
                hold on;    
                title(comparison)
                set(gca, 'XTick', 1:nr_networks);
                set(gca, 'YTick', 1:nr_networks);
                set(gca, 'XTickLabel', networks);
                set(gca, 'YTickLabel', networks);
                set(gca,'XTickLabelRotation', 30);
                colorbar

                print([fig_path, '/nbs_', num2str(nComp),'PCs_networks_',  plot_name, '_method-',method, '_cluster-', num2str(round(threshold))], '-dpng')
            
            end % comparisons
        end % significant edges metrics
    end % components
end % thresholds


%% NETWORKS LEGEND
networks_extended = {'Fronto-parietal Network (FPN)', 'Default Mode Network (DMN)', 'Dorsal Attention Network (DAN)', 'Limbic Network (LN)', 'Ventral Attention Network (VAN)',...
    'Somatomotor Network (SMN)', 'Visual Network (VN)', 'Subcortical (SUB)', 'Cerebellum (CRB)'};

line_w = 10;
figure;
L1 = plot(1:50, 'color', colors(1,:),'LineWidth', line_w);
hold on;
L2 = plot(1:50+1, 'color', colors(2,:),'LineWidth', line_w);
L3 = plot(1:50+2, 'color', colors(3,:),'LineWidth', line_w);
L4 = plot(1:50+3, 'color', colors(4,:),'LineWidth', line_w);
L5 = plot(1:50+4, 'color', colors(5,:),'LineWidth', line_w);
L6 = plot(1:50+5, 'color', colors(6,:),'LineWidth', line_w);
L7 = plot(1:50+6, 'color', colors(7,:),'LineWidth', line_w);
L8 = plot(1:50+7, 'color', colors(8,:),'LineWidth', line_w);
L9 = plot(1:50+8, 'color', colors(9,:),'LineWidth', line_w);
hl = legend([L1, L2, L3, L4, L5, L6, L7, L8, L9], networks_extended, 'Orientation', 'Vertical');
set(gca,'FontSize',13)
%hl.Position = [0.5 0.01 0.01 0.01];
print([fig_path, '/nbs_networks_legend'], '-dpng')