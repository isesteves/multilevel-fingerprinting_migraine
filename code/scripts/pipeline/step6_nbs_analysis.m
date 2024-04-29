%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 6 - NETWORK-BASED STATISTIC (NBS) ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

main_path = '/home/iesteves/FC';

nbs_stats_path = [main_path, '/data/results/nbs/stats'];
nbs_analysis_path = [main_path, '/data/results/nbs/analysis'];
atlas_path = [main_path, '/files/atlases'];
pipelinelog_path = [main_path, '/logs/pipeline_logs'];

% log file name
NBSanalysis_log_filename = [pipelinelog_path,'/step6_NBSanalysis-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%% #, % and %network of significant edges

components = [19, 68];
thresholds = [3.1, 4];
comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};
method = 'Extent';
nr_networks = 9;
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 130 9];

% Mapping of region hemispheres
t_labels_H = readtable([atlas_path, '/network_labels_hemisphere.txt'], 'HeaderLines', 0, 'ReadVariableNames', false);
mapping_H =t_labels_H{:,1}; % network and hemisphere for each region
labels_H = unique(mapping_H); % unique labels: L and R for each network + CRB vermis (19)

% store relevant information as struct fields
sigmetrics.comparisons = comparisons;
sigmetrics.contrast_types = contrast_types;
sigmetrics.method = method;
sigmetrics.nr_networks = nr_networks;
sigmetrics.areas = areas;
sigmetrics.dim = 'networks-networks-comparisons';
sigmetrics.labels_H = labels_H;

%%

% for each threshold
for t = 1:length(thresholds)
    threshold = thresholds(t);
    sigmetrics.threshold = threshold; % store threshold

    % for each number of PCs
    for comp = 1:length(components)
        nComp = components(comp);

        managefolders([nbs_analysis_path, '/analysis-', num2str(nComp), 'PCs'], 'create');
        
        count_sig_FC = zeros(nr_networks, nr_networks, length(comparisons)); % number of significant edges
        perc_sig_FC = zeros(nr_networks, nr_networks, length(comparisons)); % percentage of significant edges considering the total number of significant edges
        percarea_sig_FC = zeros(nr_networks, nr_networks, length(comparisons)); % percentage of significant edges considering the total number of possible edges for that area
        countH_sig_FC = zeros(length(labels_H), length(labels_H), length(comparisons)); % number of significant edges, split by hemisphere
        percH_sig_FC = zeros(length(labels_H), length(labels_H), length(comparisons)); % percentage of significant edges from the total of significant ones, split per hemisphere

        % for each comparison
        for c = 1:length(comparisons)
            comparison = comparisons{c};
            contrast_type = contrast_types{c};
            
            % Log file
            NBSanalysis_log_msg = ['------ ', 'Computing significant edge metrics for: ', num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(thresholds(t)))];
            disp(NBSanalysis_log_msg);
            logCustom(NBSanalysis_log_msg, NBSanalysis_log_filename)

            try
                load([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);
                adj_matrix=full(nbs.NBS.con_mat{1});
                
                % since the adjacency matrix is an upper triangular binary
                % matrix, summing gives the total # of sig edges
                nr_sig_edges = sum(adj_matrix(:));
                
                % for each network combination
                for a = 1:nr_networks
                    for b = 1:nr_networks                
                        area = adj_matrix(areas(a,1):areas(a,2), areas(b,1):areas(b,2)); % define network pair area
                        % get the number of possible edges (within-network, there is redundancy)
                        if a==b
                            size_area = size(area,1)*(size(area,1)-1)/2; % exclude repeated pairs and the diagonal (between the same regions)
                        else
                            size_area = length(area(:)); 
                        end
                        count_sig_FC(a, b, c) = sum(area(:)); 
                        perc_sig_FC(a, b, c) = (sum(area(:))/nr_sig_edges)*100; % divide by the total number of sig edges
                        percarea_sig_FC(a, b, c) = (sum(area(:))/size_area)*100; % divide by the number of possible edges within or between networks
                    end
                end
                
                % results per hemisphere
                [row,col] = find(adj_matrix);
                
                % for each significant edge, get the network/hemisphere
                % label
                for k = 1:nr_sig_edges
                    % get the networks/hemispheres corresponding to the
                    % atlas regions
                    label1 = mapping_H(row(k));
                    label2 = mapping_H(col(k));
                    
                    % get matrix indices corresponding to the involved
                    % networks/hemispheres
                    ind1 = find(strcmp(labels_H, label1));
                    ind2 = find(strcmp(labels_H, label2));
                    
                    % sort to get always the same order, for a triangular
                    % matrix
                    sorted = sort([ind1,ind2], 'descend');

                    % adds 1 to the corresponding position, to get total
                    % number of edges
                    countH_sig_FC(sorted(1), sorted(2), c) = countH_sig_FC(sorted(1), sorted(2), c)+1;
                end
                
                percH_sig_FC(:,:,c) = (countH_sig_FC(:,:,c)/nr_sig_edges)*100; %divide by the number of significant edges
            catch ME
                disp(ME)
                continue             
            end 
        end % comparisons
        
        % assign values to each field
        sigmetrics.count_sig_FC = count_sig_FC;
        sigmetrics.perc_sig_FC = perc_sig_FC;
        sigmetrics.percarea_sig_FC = percarea_sig_FC; 
        sigmetrics.countH_sig_FC = countH_sig_FC;
        sigmetrics.percH_sig_FC = percH_sig_FC;

        % store sigmetrics
        save([nbs_analysis_path, '/analysis-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_sig-metrics_networks_method-',method, '_cluster-', num2str(round(threshold))], 'sigmetrics')
        
        % Log file
        NBSanalysis_log_msg = ['------ ', 'Saved significant edge metrics for all comparisons: ', num2str(nComp), 'PCs_cluster-', num2str(round(thresholds(t)))];
        disp(NBSanalysis_log_msg);
        logCustom(NBSanalysis_log_msg, NBSanalysis_log_filename)
        
    end % components
end % thresholds

