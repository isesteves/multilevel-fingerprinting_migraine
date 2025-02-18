
config_path = '/home/iesteves/FC/code/config';

config_details = dir('/home/iesteves/FC/code/config/step1a_parcellation_config.mat');
disp(['--- Step 1 - Config file last modified on ', datestr(config_details.date)]);

% Load configuration settings
config = load([config_path, '/step1a_parcellation_config.mat']);
config = config.config;  % Extract the structure from the loaded file

% Access script-specific parameters
parcellation_inpath = config.in_path;
atlas_path = config.atlas_path;
pipelinelog_path = config.pipelinelog_path;
fig_path = config.fig_path;

subjects = config.subjects;
sessions = config.sessions;
cleanup_list = config.pipelines;
atlas_list = config.atlases;
atlas_list_N_areas = config.atlases_Nareas;

plt = config.plt;

%%

for a = 1:length(atlas_list)
    atlas = atlas_list{a};
    N_areas = atlas_list_N_areas(a);
    
    for p = 1:length(cleanup_list)
        pipeline = cleanup_list{p};

        auxfilenames = cell(length(sessions),1);
        for k = 1:length(sessions)
            kfiles = subjects{k};
            files = cellfun(@(x) ['/home/iesteves/FC/data/parcellated_data/',x,'/',sessions{k},'/percent-', pipeline, '-parcellated', atlas, '.mat'],kfiles,'UniformOutput',false);
            auxfilenames{k,1} = files; 
        end
        filenames = [auxfilenames{:}];
        
        percentage_included_all = zeros(N_areas, length(filenames));
        for s = 1:length(filenames) 
            filename = filenames{s};
            load(filename);
            percentage_included_all(:,s) = percentage_included;
        end
        
    end
end

%% Figures

figure;
plot(percentage_included_all, 'o')
ylabel('Percentage included')
xlabel('Parcel number')
line(xlim, [50 50], 'Color', 'k')
grid minor
print([fig_path, '/parcellation_nullROI_percentage_included_scatter'], '-dpng')

figure;
imagesc(percentage_included_all)
ylabel('Parcels')
xlabel('Data Samples')
print([fig_path, '/parcellation_nullROI_percentage_included_imagesc'], '-dpng')

%%
nr_subjects = length([subjects{:}]);

percentage_100 = sum(percentage_included_all==100,2);
ind_100 = find(percentage_100 < nr_subjects);
percentage_90 = sum(percentage_included_all>90,2);
ind_90 = find(percentage_90 < nr_subjects);
percentage_50 = sum(percentage_included_all>50,2);
ind_50 = find(percentage_50 < nr_subjects);

ind_100_list = cellfun(@(x) [x,', '], cellstr(num2str(ind_100)),'UniformOutput',false);
ind_90_list = cellfun(@(x) [x,', '], cellstr(num2str(ind_90)),'UniformOutput',false);
ind_50_list = cellfun(@(x) [x,', '], cellstr(num2str(ind_50)),'UniformOutput',false);


figure('pos', [50 50 1200 900]);
subplot(3,1,1)
bar(percentage_100)
xlim([1 138])
ylabel('Data Samples')
%ind_100_list_all = [ind_100_list{1:25}];
title({'Percentage included equal to 100% for all subjects NOT verified for parcels:',[ind_100_list{1:25}], [ind_100_list{25:49}]})

subplot(3,1,2)
bar(percentage_90)
xlim([1 138])
ylabel('Data Samples')
title({'Percentage included larger than 90% for all subjects NOT verified for parcels:', [ind_90_list{:}]})

subplot(3,1,3)
bar(percentage_50)
xlim([1 138])
xlabel('Parcels')
ylabel('Data Samples')
title({'Percentage included larger than 50% for all subjects NOT verified for parcels:', [ind_50_list{:}]})
print([fig_path, '/parcellation_nullROI_percentage_included_criteria'], '-dpng')

parcel2exclude = ind_50;

disp(['=== Storing parcel indices to exclude from ',atlas,' in ../file/atlases: ',  [ind_50_list{:}]])
save([atlas_path, '/parcellation_nullROI_parcel2exclude_',atlas,'_idx'], 'parcel2exclude');
