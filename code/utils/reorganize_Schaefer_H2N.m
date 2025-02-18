%% Paths and variables 
atlas_path = '/home/iesteves/FC/files/atlases';

disp('=== Reorganizing Schaefer atlas according to networks (instead of hemispheres)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reorganize Schaefer labels by network (instead of hemisphere)
% read .txt created from Yeo's lookup table corresponding to per hemisphere
% organization
t = readtable([atlas_path, '/Schaefer/Schaefer2018_100Parcels_7Networks_order.txt']);
t_region = t(:,5);

% split the each label 
C = table2cell(t_region);
C_aux = cellfun(@(x) strsplit(x, '_'), C, 'UniformOutput', 0);
C_All = cellfun(@(x) [x{3}, '_', x{4:end}], C_aux, 'UniformOutput', 0);

% reorganize according to network (alphabetically)
[C_ordered, ind] = sort({C_All{:}});

% create .txt file organized by network instead of hemisphere
new_t_region = t_region(ind,1);
new_t = t(ind,:);
writetable(new_t, [atlas_path, '/Schaefer/Schaefer2018_100Parcels_7Networks_order_new.txt'], 'WriteVariableNames', 0);

save([atlas_path, '/Schaefer/Schaefer_H2N_ind.mat'], 'ind')

