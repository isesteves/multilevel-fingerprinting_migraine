%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% STEP 7 - NETWORK METRICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

main_path = '/home/iesteves/FC';

addpath('/home/iesteves/toolboxes/BrainNetViewer_20191031') 
addpath('/home/iesteves/toolboxes/2019_03_03_BCT')

nbs_stats_path = [main_path, '/data/results/nbs/stats'];
fig_path = [main_path, '/figures/nbs_brainnetviewer'];
brainviewerfiles_path = [main_path, '/files/brainnetviewer'];
atlas_path = [main_path, '/files/atlases'];
DI_path = [main_path, '/data/results/DI']; 
network_metrics_path = '/home/iesteves/FC/data/results/network_metrics';
pipelinelog_path = [main_path, '/logs/pipeline_logs'];

managefolders(fig_path, 'create');
managefolders(network_metrics_path, 'create');

% log file name
netmetrics_log_filename = [pipelinelog_path,'/step7_netmetrics-log_', datestr(now, 'yyyy-mm-dd'),'.txt'];

%% NBS sig edges - compute local measures and store as .mat - node degree, betweenness centrality and clustering coefficient

nr_areas = 130;
metric_names = {'nodedegree', 'clustercoef', 'betweencentr'};
comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};
method = 'Extent';
dim = 'areas-comparisons';

components = [19, 68];
thresholds = [3.1, 4];

network.comparisons = comparisons;
network.contrast_types = contrast_types;
network.method = method;
network.dim_local = dim;

for  t = 1:length(thresholds)
    threshold = thresholds(t);
    network.threshold = threshold;

    for comp = 1:length(components)
        nComp = components(comp);
        
        managefolders([network_metrics_path, '/metrics_nbs-', num2str(nComp), 'PCs/'], 'create')

        for m = 1:length(metric_names)
            metric_name = metric_names{m};
            metric = zeros(nr_areas, length(comparisons));
            
            for c = 1:length(comparisons)
                comparison = comparisons{c};
                contrast_type = contrast_types{c};
                try
                    load([nbs_stats_path, '/stats-', num2str(nComp), 'PCs/nbs_',num2str(nComp), 'PCs_', comparison, '_contrast', num2str(contrast_type), '_method-',method, '_cluster-', num2str(round(threshold))]);

                    adj_matrix=full(nbs.NBS.con_mat{1}+nbs.NBS.con_mat{1}');

                    if strcmp(metric_name, 'nodedegree')
                        metric(:,c) = degrees_und(adj_matrix)';
                    elseif strcmp(metric_name, 'clustercoef')
                        metric(:,c) = clustering_coef_bu(adj_matrix);
                    elseif strcmp(metric_name, 'betweencentr')    
                        metric(:,c) = betweenness_bin(adj_matrix)';
                    end

                catch ME
                    continue
                end
            end
            network.(metric_name) = metric;
        end
        save([network_metrics_path, '/metrics_nbs-', num2str(nComp), 'PCs/network-metrics_nbs_', num2str(nComp), 'PCs_method-',method, '_cluster-', num2str(round(threshold))], 'network')
        
        % Log file
        netmetrics_log_msg = ['------ ', 'Saved network metrics for significant edges for all comparisons: ', num2str(nComp), 'PCs_cluster-', num2str(round(thresholds(t)))];
        disp(netmetrics_log_msg);
        logCustom(netmetrics_log_msg, netmetrics_log_filename)
    end
end


%% FC PCA recon - compute local and global measures
% LOCAL: 'nodedegree', 'clustercoef', 'betweencentr'
% GLOBAL: 'charpathlength', 'globaleff', 'transitivity'
load([DI_path,'/DI_FC_PCA/FC_PCA_recon.mat'])

components = [19, 68];
nr_subjects = 68;

% Thresholding - 80th percentile (only keeping 20% strongest)
thr = 80;
for comp = 1:length(components)
    nComp = components(comp);
    thr_FC_PCA_recon = FC_PCA_recon(:,:,:,nComp);
    for k = 1:nr_subjects       
        subject_FC = FC_PCA_recon(:,:,k,nComp);
        FC_thr = prctile(subject_FC, thr);
        subject_FC(subject_FC<FC_thr) = 0;
        thr_FC_PCA_recon(:,:,k) = subject_FC;
    end
    save([network_metrics_path, '/thr_FC_PCA_recon_', num2str(nComp),'PCs.mat'], 'thr_FC_PCA_recon')
    
    % Log file
    netmetrics_log_msg = ['------ ', 'Saved thresholded FC matrices keeping the ', num2str(100-thr), '% strongest connections'];
    disp(netmetrics_log_msg);
    logCustom(netmetrics_log_msg, netmetrics_log_filename)
end

% Network measures
metric_names = {'nodedegree', 'clustercoef', 'betweencentr', 'charpathlength', 'globaleff', 'transitivity'};
network.dim_local = 'areas-subjects';
networ.dim_global = 'subjects';
metric_group_nodal = zeros(nr_areas, nr_subjects);
metric_group_global = zeros(1, nr_subjects);
for comp = 1:length(components)
    nComp = components(comp);
    load([network_metrics_path, '/thr_FC_PCA_recon_', num2str(nComp),'PCs.mat']);
    
    managefolders([network_metrics_path, '/metrics_FC-', num2str(nComp), 'PCs/'], 'create')
  
    for m = 1:length(metric_names)
        metric_name = metric_names{m};
        for k = 1:nr_subjects
            subject_FC = thr_FC_PCA_recon(:,:,k);

            switch metric_name
                case 'nodedegree'
                    metric = degrees_und(subject_FC)';
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
    save([network_metrics_path, '/metrics_FC-', num2str(nComp), 'PCs/network-metrics_FC_', num2str(nComp), 'PCs'], 'network')
    
    % Log file
    netmetrics_log_msg = ['------ ', 'Saved network metrics for thresholded FC matrices reconstructed with ', num2str(nComp), 'PCs'];
    disp(netmetrics_log_msg);
    logCustom(netmetrics_log_msg, netmetrics_log_filename)
end

