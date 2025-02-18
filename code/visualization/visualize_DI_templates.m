fig_path = '/home/iesteves/FC/figures/DI_templates';
template_path = '/home/iesteves/FC/files/template_matrices';

managefolders(fig_path, 'create');

nr_subjects = 68;

%% Load template matrices
% general (for Idiff)
load([template_path, '/within_group'])
load([template_path, '/within_subject'])
load([template_path, '/within_session'])
load([template_path, '/within_menstrual_session'])
load([template_path, '/between_group_within_session'])

% subdivisions (to explore more specific aspects, for example, for specific sessions)
load([template_path, '/wgroup_corr'])
load([template_path, '/wgroup_bsession_corr'])
load([template_path, '/wsubject_bsession_corr'])
load([template_path, '/wsession_corr'])
load([template_path, '/bgroup_corr'])

%% Figure - Idiff templates
figure('pos', [200 1000 1800 300]);
subplot(1,5,1)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_group==0.5)=0;
imagesc(within_group,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Within group')

subplot(1,5,2)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_subject==0.5)=0;
imagesc(within_subject,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Within subject')

subplot(1,5,3)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_session==0.5)=0;
imagesc(within_session,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Within session')

subplot(1,5,4)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_menstrual_session==0.5)=0;
imagesc(within_menstrual_session,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title({'Within menstrual', 'session'})

subplot(1,5,5)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(between_group_within_session==0.5)=0;
imagesc(between_group_within_session,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title({'Between group','within menstrual session'})

print([fig_path, '/DI_templates_matrices_Idiff'], '-dpng')

%% Figure - correlation templates

figure('pos', [200 1000 1800 300])
subplot(1,5,1)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(wgroup_corr==0.5)=0;
imagesc(wgroup_corr,'AlphaData',imAlpha);
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Within group')

subplot(1,5,2)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(wgroup_bsession_corr==0.5)=0;
imagesc(wgroup_bsession_corr,'AlphaData',imAlpha)
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title({'Within group', 'between sessions'})

subplot(1,5,3)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(wsubject_bsession_corr==0.5)=0;
imagesc(wsubject_bsession_corr,'AlphaData',imAlpha)
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title({'Within subject', 'between sessions'})

subplot(1,5,4)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(wsession_corr==0.5)=0;
imagesc(wsession_corr,'AlphaData',imAlpha)
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Within session')

subplot(1,5,5)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(bgroup_corr==0.5)=0;
imagesc(bgroup_corr,'AlphaData',imAlpha)
xlabel('[Subjects x Sessions]')
ylabel('[Subjects x Sessions]')
title('Between groups')

print([fig_path, '/DI_templates_matrices_correlation'], '-dpng')

%% new version of Idiff templates
%figure('pos', [200 1000 1100 300]);
figure('pos', [200 900 1900 450])
subplot(1,3,1)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_subject==0.5)=0;
imagesc(within_subject,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
ylabel('Subjects x Sessions')
title('Within subject')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
colormap([0.9020 0.9020 0.9020; 0 0 0; 0.4980 0.4980 0.4980])
set(gca, 'FontSize', 14)

subplot(1,3,2)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_session==0.5)=0;
imagesc(within_session,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
%ylabel('Subjects x Sessions')
title('Within session')
colormap([0.9020 0.9020 0.9020; 0 0 0; 0.4980 0.4980 0.4980])
set(gca, 'FontSize', 14)

subplot(1,3,3)
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(within_group==0.5)=0;
imagesc(within_group,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
%ylabel('Subjects x Sessions')
title('Within group')
colormap([0.9020 0.9020 0.9020; 0 0 0; 0.4980 0.4980 0.4980])
set(gca, 'FontSize', 14)
print([fig_path, '/DI_templates_matrices_Idiff_updated'], '-dpng')

%% new version of correlation templates

%test_cmap_within_subject = [0.5059 0 0; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 
%    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

%test_cmap_between_groups = [1.0000    0.8000    0.2000; 1.0000 0.2000 0.8000; 1.0000 0.4000 0.6000; 0 1.0000 0.5000];

cmap_within_subject = [113 121 126; 132 136 132; 137 148 153; 169 169 169; 211 211 211; 229 228 226; 54 69 79]/255;

c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;
cmap_within_session = [c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle];

cmap_between_groups = [112 128 144; 115 147 179; 178 180 191; 74 80 100]/255;

figure('pos', [200 900 2100 450])
ax1 = subplot(1,3,1);
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(ismember(wsubject_bsession_corr, [0, 0.5]))=0;
imagesc(wsubject_bsession_corr,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
ylabel('Subjects x Sessions')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
title({'Within subject'})
set(gca, 'FontSize', 14);
caxis([1, max(wsubject_bsession_corr(:))]); 
colormap(ax1, cmap_within_subject);
c = colorbar;
c.Ticks = 1.4:0.85:7.2;
c.TickLabels = {['\bf \color{black}M-\color[rgb]{',num2str(c_preictal),'}pre\color{black}; \color[rgb]{',num2str(c_ictal),'}ict'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_preictal),'}pre\color{black}; \color[rgb]{',num2str(c_postictal),'}post'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_preictal),'}pre\color{black}; \color[rgb]{',num2str(c_interictal),'}inter'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_ictal),'}ict\color{black}; \color[rgb]{',num2str(c_postictal),'}post'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_ictal),'}ict\color{black}; \color[rgb]{',num2str(c_interictal),'}inter'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_postictal),'}post\color{black}; \color[rgb]{',num2str(c_interictal),'}inter'],...
    ['\bf \color{black}HC-\color[rgb]{',num2str(c_premenstrual),'}peri\color{black}; \color[rgb]{',num2str(c_midcycle),'}postov']};

ax2 = subplot(1,3,2);
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(ismember(wsession_corr, [0, 0.5]))=0;
imagesc(wsession_corr,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
%ylabel('[Subjects x Sessions]')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
title('Within session')
colormap(cmap_within_session)
set(gca, 'FontSize', 14);
caxis([1, max(wsession_corr(:))]); 
colormap(ax2, cmap_within_session); 
c = colorbar;
c.Ticks = 1.4:0.85:6.2;
c.TickLabels = {['\bf \color{black}M-\color[rgb]{',num2str(c_preictal),'}pre'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_ictal),'}ict'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_postictal),'}post'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_interictal),'}inter'],...
    ['\bf \color{black}HC-\color[rgb]{',num2str(c_premenstrual),'}peri'],...
    ['\bf \color{black}HC-\color[rgb]{',num2str(c_midcycle),'}postov']};


ax3 = subplot(1,3,3);
imAlpha=ones(nr_subjects, nr_subjects);
imAlpha(ismember(bgroup_corr, [0,0.5]))=0;
imagesc(bgroup_corr,'AlphaData',imAlpha);
xlabel('Subjects x Sessions')
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
%ylabel('[Subjects x Sessions]')
title('Between groups')
set(gca, 'FontSize', 14);
colormap(ax3, cmap_between_groups);
caxis([1, max(bgroup_corr(:))]); % Set color axis limits for the first subplot
c = colorbar;
c.Ticks = 1.4:0.75:4.2;
c.TickLabels = {['\bf \color{black}M-\color[rgb]{',num2str(c_preictal),'}pre\color{black}; HC-\color[rgb]{',num2str(c_premenstrual),'}peri'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_ictal),'}ict\color{black}; HC-\color[rgb]{',num2str(c_premenstrual),'}peri'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_postictal),'}post\color{black}; HC-\color[rgb]{',num2str(c_premenstrual),'}peri'],...
    ['\bf \color{black}M-\color[rgb]{',num2str(c_interictal),'}inter\color{black}; HC-\color[rgb]{',num2str(c_midcycle),'}postov']};

print([fig_path, '/DI_templates_matrices_correlation_updated'], '-dpng')

