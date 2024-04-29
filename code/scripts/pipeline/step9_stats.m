addpath('/home/iesteves/stats/MultipleTestingToolbox/MultipleTestingToolbox')

fig_path = '/home/iesteves/FC/figures/figures_paper';
template_path = '/home/iesteves/FC/files/template_matrices';
DI_path = '/home/iesteves/FC/data/results/DI';

managefolders(fig_path, 'create');

nr_subjects = 68;
nr_patients = 10;
nr_controls = 14;

n = 10;
comp = 19;

%% Load correlation matrices
load([DI_path, '/DI_corr/corr_val_wsubject_bsession.mat']);
load([DI_path, '/DI_corr/corr_val_wsession.mat']);
load([DI_path, '/DI_corr/corr_val_bgroup.mat']);

%% within subject 

corr_wsubject = corr_val_wsubject_bsession{n, comp};
corr_wsubject_patients = mean(corr_wsubject(1:nr_patients, 1:6), 2);
corr_wsubject_controls = corr_wsubject(1:nr_controls,7);

[p,h] = ranksum(corr_wsubject_patients, corr_wsubject_controls);

% %% within session comparisons - inside HC and inside MIG
% 
% corr_wsession = corr_val_wsession{n, comp};
% 
% vec1 = [5, 1, 1, 1, 2, 2, 3];
% vec2 = [6, 2, 3, 4, 3, 4, 4];
% ncomp = length(vec1);
% hvec = zeros(ncomp, 1);
% pvec = zeros(ncomp, 1);
% sig = zeros(ncomp, 2);
% for i = 1:ncomp
%     [h,p] = ttest(corr_wsession(:,vec1(i)), corr_wsession(:,vec2(i)));
%     hvec(i) = h;
%     pvec(i) = p;
% end
% sig(:,1) = pvec < 0.05/ncomp;
% sig(:,2) = pvec < 0.01/ncomp;
% disp(sig)
% 
% col_names = {'sig', 'verysig'};
% row_names = {'pm-postov', 'pre-ict', 'pre-post', 'pre-inter', 'ict-post', 'ict-inter', 'post-inter'};
% array2table(sig, 'VariableNames', col_names, 'RowNames', row_names)

%% between session comparisons HC-MIG vs HC-MIG - not averaged

corr_bgroup = corr_val_bgroup{n, comp};

vec1 = [1, 1, 1, 2, 2, 3]; 
vec2 = [2, 3, 4, 3, 4, 4];
ncomp = length(vec1);
hvec = zeros(ncomp, 1);
pvec = zeros(ncomp, 1);
sig = zeros(ncomp, 2);
for i = 1:ncomp
    [h,p] = ttest(corr_bgroup(:,vec1(i)), corr_bgroup(:,vec2(i)));
    hvec(i) = h;
    pvec(i) = p;
end
sig(:,1) = pvec < 0.05/ncomp;
sig(:,2) = pvec < 0.01/ncomp;
disp(sig)
col_names = {'sig', 'verysig'};
row_names = {'pre-ict', 'pre-post', 'pre-inter', 'ict-post', 'ict-inter', 'post-inter'};
array2table(sig, 'VariableNames', col_names, 'RowNames', row_names)

%% between session - average

corr_bgroup_aux = corr_val_bgroup{n, comp};

corr_bgroup = zeros(10,4);

for s =1:4

    corr_session_vec = corr_bgroup_aux(:,s);
    corr_session = reshape(corr_session_vec, 14, 10);

    m_val_patient = mean(corr_session)';
    corr_bgroup(:,s) = m_val_patient;
end

vec1 = [1, 1, 1, 2, 2, 3]; 
vec2 = [2, 3, 4, 3, 4, 4];
ncomp = length(vec1);
hvec = zeros(ncomp, 1);
pvec = zeros(ncomp, 1);
sig = zeros(ncomp, 2);
for i = 1:ncomp
    [p,h] = signrank(corr_bgroup(:,vec1(i)), corr_bgroup(:,vec2(i)));
    hvec(i) = h;
    pvec(i) = p;
end
sig(:,1) = pvec < 0.05/ncomp;
sig(:,2) = pvec < 0.01/ncomp;
disp(sig)
col_names = {'sig', 'verysig'};
row_names = {'pre-ict', 'pre-post', 'pre-inter', 'ict-post', 'ict-inter', 'post-inter'};
array2table(sig, 'VariableNames', col_names, 'RowNames', row_names)

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec, 0.05, true);  

array2table([pvec, c_pvalues', (c_pvalues<0.05)'], 'VariableNames',{'pval', 'corrpval', 'FDR'}, 'RowNames', row_names)

%% stats within session - not averaged

corr_wsession = corr_val_wsession{n, comp};

vec1 = [5, 1, 1, 1, 2, 2, 3, 5, 5, 5, 6];
vec2 = [6, 2, 3, 4, 3, 4, 4, 1, 2, 3, 4];
ncomp = length(vec1);
hvec = zeros(ncomp,1);
pvec = zeros(ncomp,1);
sig = zeros(ncomp,2);
% stats
for i = 1:ncomp
    if ncomp < 8 
        [h,p] = ttest(corr_wsession(:,vec1(i)), corr_wsession(:,vec2(i)));
        hvec(i) = h;
        pvec(i) = p;
    elseif ncomp >= 8
        [h,p] = ttest2(corr_wsession(:,vec1(i)), corr_wsession(:,vec2(i)));
        hvec(i) = h;
        pvec(i) = p;
    end
end

sig(:,1) = pvec < 0.05/ncomp;
sig(:,2) = pvec < 0.01/ncomp;
disp(sig)

col_names = {'sig', 'verysig'};
row_names = {'pm-postov', 'pre-ict', 'pre-post', 'pre-inter', 'ict-post', 'ict-inter', 'post-inter', 'pre-pm', 'ict-pm', 'post-pm', 'inter-postov'};
array2table(sig, 'VariableNames', col_names, 'RowNames', row_names)


%% stats within session - average

corr_wsession = corr_val_wsession{n, comp};

nr_val_patients = (nr_patients*nr_patients-nr_patients)/2;
nr_val_controls = (nr_controls*nr_controls-nr_controls)/2;

avg_corr_patients = NaN*ones(nr_controls, 4);
for s = 1:4
    corr_wsession_m = conversion_vecmat(corr_wsession(1:nr_val_patients,s), nr_patients, 'vec2mat');
    corr_wsession_m(1:size(corr_wsession_m, 1)+1:end) = NaN;
    avg_corr_patients(1:nr_patients,s) = nanmean(corr_wsession_m);
end

avg_corr_controls = zeros(nr_controls, 2);

for s = 1:2
    corr_wsession_m = conversion_vecmat(corr_wsession(:,s+4), nr_controls, 'vec2mat');
    corr_wsession_m(1:size(corr_wsession_m, 1)+1:end) = NaN;
    avg_corr_controls(:,s) = nanmean(corr_wsession_m);
end

corr_wsession = [avg_corr_patients, avg_corr_controls];

vec1 = [5, 1, 1, 1, 2, 2, 3, 5, 5, 5, 6];
vec2 = [6, 2, 3, 4, 3, 4, 4, 1, 2, 3, 4];
ncomp = length(vec1);
hvec = zeros(ncomp,1);
pvec = zeros(ncomp,1);
sig = zeros(ncomp,3);
% stats
for i = 1:ncomp
    if ncomp < 8 
        [p,h] = signrank(corr_wsession(:,vec1(i)), corr_wsession(:,vec2(i)));
        hvec(i) = h;
        pvec(i) = p;
    elseif ncomp >= 8
        [p,h] = ranksum(corr_wsession(:,vec1(i)), corr_wsession(:,vec2(i)));
        hvec(i) = h;
        pvec(i) = p;
    end
end

sig(:,1) = pvec < 0.05;
sig(:,2) = pvec < 0.05/ncomp;
sig(:,3) = pvec < 0.01/ncomp;
disp(sig)

col_names = {'uncorr', 'sig', 'verysig'};
row_names = {'pm-postov', 'pre-ict', 'pre-post', 'pre-inter', 'ict-post', 'ict-inter', 'post-inter', 'pre-pm', 'ict-pm', 'post-pm', 'inter-postov'};
array2table(sig, 'VariableNames', col_names, 'RowNames', row_names)

[c_pvalues, c_alpha, h_BH] = fdr_BH(pvec, 0.05, true);  

array2table([pvec, c_pvalues', (c_pvalues<0.05)'], 'VariableNames',{'pval', 'corrpval', 'FDR'}, 'RowNames', row_names)
