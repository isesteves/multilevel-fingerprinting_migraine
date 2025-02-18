load('/home/iesteves/FC/data/results/network_metrics/metrics_nbs-19PCs/network-metrics_nbs_19PCs_method-Extent_cluster-4.mat')

figure; bar([network.nodedegree(:,5),network.nodedegree(:,6), network.nodedegree(:,7)], 'stacked')
legend('ict-pm', 'postic-pm', 'pre-ict')
ylabel('node degree')
xlabel('regions')
set(gca, 'FontSize', 14)

atlas = readtable('/home/iesteves/FC/files/atlases/coord_SchaeferSubCRB7100_130.txt');
labels = atlas{:,6};

%
labels([50,76,120,79,49]);

%
labels([42, 67, 96, 87, 15]);

%
labels([15, 38, 42, 65, 87, 96, 97, 99]);

% pre-ict
labels([113, 90, 114, 119, 122])
