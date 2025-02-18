addpath('/home/iesteves/toolboxes/BrainNetViewer_20191031');
network_metrics_path = '/home/iesteves/FC/data/results/network_metrics';
nbs_stats_path = '/home/iesteves/FC/data/results/nbs/stats';
fig_path = '/home/iesteves/FC/figures/nbs_brainnetviewer';
brainviewerfiles_path = '/home/iesteves/FC/files/brainnetviewer';
atlas_path = '/home/iesteves/FC/files/atlases';
nr_areas = 130;

surface = '/home/iesteves/toolboxes/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2withCerebellum.nv';

comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};
thresholds = [3.1, 4];
method = 'Extent';
metric_names = {'nodedegree', 'clustercoef', 'betweencentr'};
components = 68;
  
%% Networks matrix

matrix_all = zeros(nr_areas, nr_areas, 9);
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 130 9];
for g = 1:size(areas,1)
    matrix_all(areas(g,1):areas(g,2),areas(g,1):areas(g,2),g) = g;
    matrix_all(:,:,g) = matrix_all(:,:,g) - diag(diag(matrix_all(:,:,g)));
end

matrix_FC_all = sum(matrix_all,3);

T1 = table(matrix_FC_all);
writetable(T1,[brainviewerfiles_path, '/edges_SchaeferSubCRB7100_130_networks.txt'],'Delimiter',' ', 'WriteVariableNames', 0)  
file = fullfile([brainviewerfiles_path, '/edges_SchaeferSubCRB7100_130_networks.txt']);
[tempDir, tempFile] = fileparts(file); 
status = copyfile(file, fullfile(tempDir, [tempFile, '.edge']));
delete(file)

%% create files 
length_Schaefer = 100;
t1 = readtable([atlas_path, '/coord_SchaeferSubCRB7100_130.txt']);
full_labels = t1{1:length_Schaefer,6};
labels_aux = cellfun(@(x) strsplit(x, '_'), full_labels, 'UniformOutput', 0);
labels = cellfun(@(x) [x{:,end-1},'-',x{:,end}], labels_aux, 'UniformOutput', 0);
t1{1:length_Schaefer,6} = labels;
writetable(t1, [brainviewerfiles_path, '/nodes_SchaeferSubCRB7100_130.txt'], 'Delimiter', ' ', 'WriteVariableNames', 0)

for t = 1:length(thresholds)
    threshold = thresholds(t);
    for comp = 1:length(components)
        nComp = components(comp);

        load([network_metrics_path, '/metrics_nbs-', num2str(nComp), 'PCs/network-metrics_nbs_', num2str(nComp), 'PCs_method-',method, '_cluster-', num2str(round(threshold))])

        for c = 1:length(comparisons)
            comparison = comparisons{c};
            contrast_type = contrast_types{c};
            
            try
                load([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
            catch ME
                continue
            end
            adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');

            managefolders([brainviewerfiles_path, '/edges-', num2str(nComp), 'PCs/'], 'create')
            dlmwrite([brainviewerfiles_path, '/edges-', num2str(nComp), 'PCs', ...
                '/edges_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '.txt'], adj_matrix, 'delimiter', ' ');

            file = fullfile([brainviewerfiles_path, '/edges-', num2str(nComp), 'PCs/'], ['edges_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '.txt']);
            [tempDir, tempFile] = fileparts(file); 
            status = copyfile(file, fullfile(tempDir, [tempFile, '.edge']));
            %delete(file)

            for m = 1:length(metric_names)
                metric_name = metric_names{m};

                managefolders([brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/', metric_name], 'create');
                metric = round(normalize_zero_one(network.(metric_name))*100, 2)+0.5;

                t_nodes = readtable([brainviewerfiles_path, '/nodes_SchaeferSubCRB7100_130.txt']);
                t_nodes{:,5} = metric(:,c);
                writetable(t_nodes, [brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/', metric_name, '/nodes_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name, '.txt'], 'Delimiter', ' ', 'WriteVariableNames', 0)

                file = fullfile([brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/', metric_name], ['nodes_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name, '.txt']);
                [tempDir, tempFile] = fileparts(file); 
                status = copyfile(file, fullfile(tempDir, [tempFile, '.node']));
                delete(file)

            end % metrics
        end % comparisons
    end % components
end % thresholds


%% figures only nodes

metric_names = {'nodedegree'};

for t = 1:length(thresholds)
    threshold = thresholds(t);
    for comp = 1:length(components)
        nComp = components(comp);
        for m = 1:length(metric_names)
            metric_name = metric_names{m};
            for c = 1:length(comparisons)
                comparison = comparisons{c};
                contrast_type = contrast_types{c};
                
               
                    nodes = [brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/',...
                        metric_name, '/nodes_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name, '.node'];

                    config = [brainviewerfiles_path, '/config-new_alpha-035.mat'];
                    figure_name = [fig_path, '/nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold)), '.png'];
                    tic
                    try
                        BrainNet_MapCfg(surface, nodes, config, figure_name);
                    catch ME
                        continue
                    end
                    toc

                    figure;
                    im1 = imread(figure_name);
                    imshow(im1)
                    h = title(['SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name]);
                    h.Interpreter = 'none';
                    set(gca, 'FontSize', 14);
                    print([fig_path, '/nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold))], '-dpng')
                    close all
               
            end % metrics
        end % comparisons
    end % components
end % thresholds

%% figures only nodes - coronal view

metric_names = {'nodedegree'};

for t = 1:length(thresholds)
    threshold = thresholds(t);
    for comp = 1:length(components)
        nComp = components(comp);
        for m = 1:length(metric_names)
            metric_name = metric_names{m};
            for c = 1:length(comparisons)
                comparison = comparisons{c};
                contrast_type = contrast_types{c};
                
               
                    nodes = [brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/',...
                        metric_name, '/nodes_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name, '.node'];

                    config = [brainviewerfiles_path, '/config-new_alpha-035_coronal.mat'];
                    figure_name = [fig_path, '/nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold)), '_coronal.png'];
                    tic
                    try
                        BrainNet_MapCfg(surface, nodes, config, figure_name);
                    catch ME
                        continue
                    end
                    toc

                    figure;
                    im1 = imread(figure_name);
                    imshow(im1)
                    h = title(['SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name]);
                    h.Interpreter = 'none';
                    set(gca, 'FontSize', 14);
                    print([fig_path, '/nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold)), '_coronal'], '-dpng')
                    close all
               
            end % metrics
        end % comparisons
    end % components
end % thresholds

%% figures nodes and edges

metric_names = {'nodedegree'};

for t = 1:length(thresholds)
    threshold = thresholds(t);        
    for comp = 1:length(components)
        nComp = components(comp);
        for m = 1:length(metric_names)
            metric_name = metric_names{m};
            for c = 1:length(comparisons)
                comparison = comparisons{c};
                contrast_type = contrast_types{c};

                nodes = [brainviewerfiles_path, '/nodes-', num2str(nComp), 'PCs/',...
                    metric_name, '/nodes_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name, '.node'];

                edges = [brainviewerfiles_path, '/edges-', num2str(nComp), 'PCs',...
                    '/edges_SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '.edge'];
                
                config = [brainviewerfiles_path, '/config-new_alpha-035.mat'];
                figure_name = [fig_path, '/edges-nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold)), '.png'];
                tic
                try
                    BrainNet_MapCfg(surface, nodes, edges, config, figure_name);
                catch ME
                     continue
                end
                toc

                figure;
                im1 = imread(figure_name);
                imshow(im1)
                h = title(['SchaeferSubCRB7100_130_nbs_',num2str(nComp), 'PCs_',comparison, '_', num2str(round(threshold)), '_', metric_name]);
                h.Interpreter = 'none';
                set(gca, 'FontSize', 14);
                print([fig_path, '/edges-nodes_SchaeferSubCRB7100_130_',metric_name, '_nbs_', num2str(nComp), 'PCs_',comparison, '_cluster-', num2str(round(threshold))], '-dpng')
                close all
                
            end % metrics
        end % comparisons
    end % components
end % thresholds


%% tests brainviewer
% surface = '/home/iesteves/toolboxes/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_Ch2withCerebellum.nv';
% nodes = '/home/iesteves/FC/files/brainnetviewer/nodes_SchaeferSubCRB7100_130.node';
% edges = '/home/iesteves/FC/files/brainnetviewer/edges_SchaeferSubCRB7100_130_networks.edge';
% config = '/home/iesteves/FC/files/brainnetviewer/config-new.mat';     
% BrainNet_MapCfg(surface, nodes, edges, config);