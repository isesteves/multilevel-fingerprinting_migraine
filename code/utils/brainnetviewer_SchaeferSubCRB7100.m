nr_areas = 138;

% load Schaefer (reorganized by network)
atlas = readtable('/home/iesteves/dFC/parcellation/Schaefer2018_100Parcels_7Networks_order_new.txt');

% load original Schaefer centroid MNI coordinates and reorganize them by
% network
centroids = readtable('/home/iesteves/dFC/parcellation/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');
order = atlas.Var1;
new_centroids = centroids(order, :);

% load AAL116 coordinates
atlas_aal = readtable('/home/iesteves/dFC/parcellation/AAL116_coordinates_brainGraph.xlsx'); 
SubCRB_areas = [37 38 41 42 71:78 91:116];

atlas_aal.name  = cellfun(@(x) replace(x, '.', '_'), atlas_aal.name, 'UniformOutput', false);

labels = atlas.Var5;
labels_all = cell(nr_areas,1);
labels_all(1:100) = labels;
labels_all(101:end,:) = atlas_aal{SubCRB_areas, 2};

coordinates_all = zeros(nr_areas,3);
coordinates_all(1:100,:) = centroids{:,3:5};
coordinates_all(101:end,:) = atlas_aal{SubCRB_areas, 3:5};

areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 138 9];

colors_all = zeros(nr_areas, 1);
for k = 1:size(areas, 1)
    colors_all(areas(k,1):areas(k,2)) = areas(k,3);
end

size_all = 3*ones(nr_areas, 1);

T = table(coordinates_all, colors_all, size_all, labels_all);
writetable(T,'/home/iesteves/dFC/allnodes.txt','Delimiter',' ', 'WriteVariableNames', 0);

ind_50 =[115, 116, 123, 124, 125, 126, 127, 128];%[119, 120, 121, 122, 123, 124];
exclude_ind = ind_50;
include_ind = setdiff(1:nr_areas, exclude_ind);

T_included = table(coordinates_all(include_ind,:), colors_all(include_ind,:), size_all(include_ind,:), labels_all(include_ind,:));
writetable(T_included,'/home/iesteves/dFC/allnodes_SchaeferSubCRB7100_130.txt','Delimiter',' ', 'WriteVariableNames', 0)  

%% Networks matrix

matrix_all = zeros(nr_areas, nr_areas, 9);
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 138 9];
for g = 1:size(areas,1)
    matrix_all(areas(g,1):areas(g,2),areas(g,1):areas(g,2),g) = g;
    matrix_all(:,:,g) = matrix_all(:,:,g) - diag(diag(matrix_all(:,:,g)));
end

matrix_FC_all = sum(matrix_all,3);

T1 = table(matrix_FC_all);
writetable(T1,'alledges.txt','Delimiter',' ', 'WriteVariableNames', 0)  

