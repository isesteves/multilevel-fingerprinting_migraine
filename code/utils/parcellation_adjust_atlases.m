atlas_path = '/home/iesteves/FC/files/atlases/';

%% Desikan
gunzip([atlas_path, 'Desikan/Desikan_MNI152-2mm.nii.gz']);
atlas = load_untouch_nii([atlas_path, 'Desikan/Desikan_MNI152-2mm.nii']);

Desikan_original = atlas.img;

% Set white matter as zero (1 - white matter L; 36 - white matter R)
new_Desikan = Desikan_original;
new_Desikan(Desikan_original==1) = 0;
new_Desikan(Desikan_original==36) = 0;

inds = unique(new_Desikan);
for k = 2:inds(end)
    if any(new_Desikan(Desikan_original==k) == k) & (k<37)
        new_Desikan(Desikan_original==k) = k-1;
    elseif any(new_Desikan(Desikan_original==k) == k) & (k>=37)
        new_Desikan(Desikan_original==k) = k-2;
    end
end

figure;
imagesc(new_Desikan(:, :, 40))

adjustedDesikan = new_Desikan;
save([atlas_path, 'Desikan/adjustedDesikan.mat'], 'adjustedDesikan')

%% AAL90

gunzip([atlas_path, 'AAL/AAL116_2mm_woSchaefer.nii.gz']);
atlas = load_untouch_nii([atlas_path, 'AAL/AAL116_2mm_woSchaefer.nii']);

AAL_original = atlas.img;

AAL90 = AAL_original;
AAL90(AAL_original>90) = 0;

figure;
imagesc(AAL90(:, :, 40))

save([atlas_path, 'AAL/AAL90_2mm_woSchaefer.mat'], 'AAL90')

