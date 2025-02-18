atlas_path = '/home/iesteves/FC/files/atlases';

length_Schaefer = 100;
t1 = readtable([atlas_path, '/coord_SchaeferSubCRB7100_130.txt']);
full_labels = t1{1:length_Schaefer,6};
labels_aux = cellfun(@(x) strsplit(x, '_'), full_labels, 'UniformOutput', 0);
labels_H = cellfun(@(x) [x{:,3},'-',x{:,2}], labels_aux, 'UniformOutput', 0);

full_labels_AAL = t1{101:112,6};
labels_aux = cellfun(@(x) strsplit(x, '_'), full_labels_AAL, 'UniformOutput', 0);
labels_H_SUB = cellfun(@(x) ['SUB-',x{:,2}], labels_aux, 'UniformOutput', 0);

full_labels_AAL = t1{113:122,6};
labels_aux = cellfun(@(x) strsplit(x, '_'), full_labels_AAL, 'UniformOutput', 0);
labels_H_CRB = cellfun(@(x) ['CRB-',x{:,2}], labels_aux, 'UniformOutput', 0);

full_labels_AAL = t1{123:130,6};
labels_aux = cellfun(@(x) strsplit(x, '_'), full_labels_AAL, 'UniformOutput', 0);
labels_H_CRB_vermis = cellfun(@(x) ['CRB-VERMIS'], labels_aux, 'UniformOutput', 0);


t2 = {labels_H{:}, labels_H_SUB{:}, labels_H_CRB{:}, labels_H_CRB_vermis{:}}';
t2{117,1}='CRB-L';
t2{118,1}='CRB-R';

t3 = cell2table(t2);
writetable(t3, [atlas_path, '/network_labels_nichord.txt'], 'Delimiter', ' ', 'WriteVariableNames', 0)
