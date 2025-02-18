DI_path = '/home/iesteves/FC/data/results/DI';
fig_path = '/home/iesteves/FC/figures/DI_correlation';

managefolders(fig_path, 'create');

nr_subjects = 68;

% number of networks (including the whole brain - WB)
networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'Subcortical', 'Cerebellum', 'WB'};
nr_networks = length(networks);

% number of PCs, from 1 to the number of data samples
components = 1:nr_subjects;
nr_components = length(components);

load([DI_path, '/DI_FC_PCA/corr_PCA_recon_networks.mat']);

load([DI_path, '/DI_values/within_session_DI.mat']);
load([DI_path, '/DI_values/within_subject_DI.mat']);
load([DI_path, '/DI_values/within_group_DI.mat']);
load([DI_path, '/DI_values/within_menstrual_session_DI.mat']);
load([DI_path, '/DI_values/between_group_within_session_DI.mat']);

load([DI_path, '/DI_corr/corr_val_wgroup.mat']);
load([DI_path, '/DI_corr/corr_val_wgroup_bsession.mat']);
load([DI_path, '/DI_corr/corr_val_wsubject_bsession.mat']);
load([DI_path, '/DI_corr/corr_val_wsession.mat']);
load([DI_path, '/DI_corr/corr_val_bgroup.mat']);

% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;

colors_aux = [c_preictal; c_ictal; c_postictal; c_interictal; c_premenstrual; c_midcycle];
colors = colors_aux(end:-1:1,:);

%%
for n = 1:nr_networks
    fig = figure('pos', [200, 1000, 2000, 600]);
    suptitle({[networks{n}, ' - Differential Identifiability: Correlation values'], ' ', ' '})
    for m_ind = 2:nr_components

        corr_matrix_original = corr_PCA_recon_networks(:,:,n, nr_subjects);
        corr_matrix = corr_PCA_recon_networks(:,:,n, m_ind);

        within_session_Idiff = squeeze(within_session_DI(1,n,:));
        within_subject_Idiff = squeeze(within_subject_DI(1,n,:));
        within_group_Idiff = squeeze(within_group_DI(1,n,:));
        within_menstrual_session_Idiff = squeeze(within_menstrual_session_DI(1,n,:));
        between_group_within_session_Idiff = squeeze(between_group_within_session_DI(1,n,:));

        [y_ma_within_group, ma_within_group] = max(within_group_Idiff);
        [y_ma_within_subject, ma_within_subject] = max(within_subject_Idiff);
        [y_ma_within_session, ma_within_session] = max(within_session_Idiff);
        [~, mi_within_menstrual_session] = min(within_menstrual_session_Idiff);
        [~, mi_between_group_within_session] = min(within_menstrual_session_Idiff);

        subplot(2,6,1)   
        imagesc(corr_matrix);
        xlabel('[Subjects x Sessions]')
        ylabel('[Subjects x Sessions]')
        title('PCA Reconstructed')

        subplot(2,6,7)
        imagesc(corr_matrix_original);
        xlabel('[Subjects x Sessions]')
        ylabel('[Subjects x Sessions]')
        title('Original')

        subplot(2,6,2)
        plot(within_group_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        hold on;
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([1 nr_components]);
        xlabel('# PCs')
        ylabel('Idiff')
        title('group')
        hold off;

        subplot(2,6,8)
        boxplot(corr_val_wgroup{n, m_ind})
        hold on;
        plot_bp_pnts(corr_val_wgroup{n, m_ind}, 5);
        hold off;
        title('within group')
        ylabel('Correlation')
        set(gca, 'XTickLabel', {'MIG', 'HC'})

        subplot(2,6,3)
        plot(within_subject_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        title('subject')
        hold on;
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([1 nr_components])
        xlabel('# PCs')
        hold off;

        subplot(2,6,9)
        y1 = 0.2;
        y2 = 0.6;
        patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        hold on;
        patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
        boxplot(corr_val_wgroup_bsession{n,m_ind})
        hold on;
        plot_bp_pnts(corr_val_wgroup_bsession{n,m_ind}, 5)
        ylim([0.2 1]);
        set(gca, 'XTickLabel', {'ic-preic', 'postic-preic', 'interic-preic', 'ic-postic', 'interic-ic','interic-postic', 'postov-perimens'})
        set(gca, 'XTickLabelRotation', 30)
        title('w/ group, b/ sessions')
        hold off;
        
        subplot(2,6,4)
        plot(within_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([1 nr_components])
        xlabel('# PCs')
        title('session')

        subplot(2,6,10)  
        y1 = 0.2;
        y2 = 0.6;
        patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        hold on;
        patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_preictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [4.5 y1; 5.5 y1; 5.5 y2; 4.5 y2], 'FaceColor', c_ictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [4.5 y2; 5.5 y2; 5.5 1; 4.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [5.5 y1; 6.5 y1; 6.5 y2; 5.5 y2], 'FaceColor', c_postictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [5.5 y2; 6.5 y2; 6.5 1; 5.5 1], 'FaceColor', c_interictal, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [6.5 y1; 7.5 y1; 7.5 y2; 6.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', .5, 'EdgeColor', 'none')
        patch('Faces', [1 2 3 4], 'Vertices', [6.5 y2; 7.5 y2; 7.5 1; 6.5 1], 'FaceColor', c_midcycle, 'FaceAlpha', .5, 'EdgeColor', 'none')
        boxplot(corr_val_wsubject_bsession{n,m_ind})
        hold on;
        plot_bp_pnts(corr_val_wsubject_bsession{n,m_ind}, 5)
        ylim([0.2 1])
        set(gca, 'XTickLabel', {'ic-preic', 'postic-preic', 'interic-preic', 'ic-postic', 'interic-ic','interic-postic', 'postov-perimens'})
        set(gca, 'XTickLabelRotation', 30)
        title('w/ subject, b/ sessions')
        hold off;
        
        subplot(2,6,5)
        plot(within_menstrual_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        hold on;
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([1 nr_components])
        title('cycle')
        xlabel('# PCs')
        hold off;

        subplot(2,6,11)
        boxplot(corr_val_wsession{n,m_ind})
        hold on;
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
        end
        plot_bp_pnts(corr_val_wsession{n,m_ind}, 5)
        set(gca, 'XTickLabel', {'preic', 'ic', 'postic', 'interic', 'perimens','postov'})
        set(gca, 'XTickLabelRotation', 30)
        title('within session')
        hold off

        subplot(2,6,6)
        plot(between_group_within_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        xlim([1 nr_components])
        xlabel('# PCs')
        hold on;
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        title('group-cycle')
        hold off;

        subplot(2,6,12)
        boxplot(corr_val_bgroup{n,m_ind})
        hold on;
        g = findobj(gca,'Tag','Box');
        for j=1:length(g)
            patch(get(g(j),'XData'),get(g(j),'YData'),colors(j+2,:),'FaceAlpha',.5);
        end
        plot_bp_pnts(corr_val_bgroup{n,m_ind}, 5)
        set(gca, 'XTickLabel', {'preic-perimens', 'ic-perimens', 'postic-perimens', 'interic-postov'})
        set(gca, 'XTickLabelRotation', 30)
        title('between groups')
        hold off

        %pause(0.1)
        drawnow %limitrate
        frame = getframe(fig);
        im{m_ind-1} = frame2im(frame);

    end
    filename = [fig_path, '/DI_correlation_boxplots-gif_', networks{n},'.gif']; % Specify the output file name
    for idx = 1:nr_components-1
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end


%% NEW Figure - panel B - changed order to start with controls

%%%%%%%%%%%%%%%%%%%%%%%%%% Note: Correlation instead of average correlation

n = 10; % whole brain
m_ind = 19;
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

% within subject
corr_val_wsubject_bsession_aux = corr_val_wsubject_bsession{n,m_ind};
corr_val_wsubject_bsession_patients = corr_val_wsubject_bsession_aux(:, 1:6);
corr_val_wsubject_bsession_controls_aux = corr_val_wsubject_bsession_aux(:, 7);
corr_val_wsubject_bsession_patients_bp = corr_val_wsubject_bsession_patients(:);
corr_val_wsubject_bsession_controls_bp = [corr_val_wsubject_bsession_controls_aux; NaN*ones(length(corr_val_wsubject_bsession_patients_bp)-length(corr_val_wsubject_bsession_controls_aux),1)];
corr_val_wsubject_bsession_bp = [corr_val_wsubject_bsession_controls_bp, corr_val_wsubject_bsession_patients_bp];

% within session
corr_val_wsession_aux = corr_val_wsession{n,m_ind};
corr_val_wsession_bp = corr_val_wsession_aux(:,[5:6 1:4]);

% between group
corr_val_bgroup_bp = corr_val_bgroup{n,m_ind};

figure('pos', [50, 50, 250 1000]);
subplot(4,1,1)
boxplot(corr_val_wsubject_bsession_bp)
hold on;
plot_bp_pnts(corr_val_wsubject_bsession_bp, 10);
set(gca, 'XTickLabel', {'HC', 'MIG'})
set(gca, 'FontSize', 13)
ylim([0 1])
ylabel(' ')
grid minor
title('within subject')

subplot(4,1,2)
boxplot(corr_val_wsession_bp)
hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end
plot_bp_pnts(corr_val_wsession_bp, 5)
set(gca, 'XTickLabel', [])
set(gca, 'FontSize', 13)
ylim([0 1])
grid minor
title('within session')
ylabel('Average correlation')

colors_aux_2 = [c_preictal; c_ictal; c_postictal; c_interictal];
colors_2 = colors_aux_2(end:-1:1,:);

subplot(4,1,3)
boxplot(corr_val_bgroup_bp)
hold on;
g = findobj(gca,'Tag','Box');
for j=1:length(g)
    patch(get(g(j),'XData'),get(g(j),'YData'),colors_2(j,:),'FaceAlpha',.5);
end
plot_bp_pnts(corr_val_bgroup_bp, 5)
set(gca, 'XTickLabel', [])
set(gca, 'FontSize', 13)
ylim([0 1])
ylabel(' ')
grid minor
title('between groups')

print([fig_path, '/DI_correlation_BrainFunctionWorkshop-panelB_Correlation'],'-dpng');

%% 
n = 10; 
m_ind = 19;
colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

y1 = 0;
y2 = 0.5;
face_alpha = 0.5;

% within subject
corr_val_wsubject_bsession_aux = corr_val_wsubject_bsession{n,m_ind};
corr_val_wsubject_bsession_patients = corr_val_wsubject_bsession_aux(:, 1:6);
corr_val_wsubject_bsession_controls_aux = corr_val_wsubject_bsession_aux(:, 7);
corr_val_wsubject_bsession_patients_bp = corr_val_wsubject_bsession_patients(:);
corr_val_wsubject_bsession_controls_bp = [corr_val_wsubject_bsession_controls_aux; NaN*ones(length(corr_val_wsubject_bsession_patients_bp)-length(corr_val_wsubject_bsession_controls_aux),1)];
corr_val_wsubject_bsession_bp = [corr_val_wsubject_bsession_controls_bp, corr_val_wsubject_bsession_patients_bp];

% within session
corr_val_wsession_aux = corr_val_wsession{n,m_ind};
corr_val_wsession_bp = corr_val_wsession_aux(:,[5:6 1:4]);

% between group
corr_val_bgroup_bp = corr_val_bgroup{n,m_ind};

figure('pos', [50, 50, 1200 500]);
subplot(1,3,1)
b = boxplot(corr_val_wsubject_bsession_bp, 'Color', 'k');
set(b,{'linew'},{1.5})
hold on;
plot_bp_pnts(corr_val_wsubject_bsession_bp, 10);
set(gca, 'XTickLabel', {'HC', 'MIG'})
set(gca, 'FontSize', 13)
ylim([0 1])
ylabel('Correlation')
grid minor
title('within subject')
hAx=gca;
hAx.XAxis.TickLabelInterpreter='tex';            % make%hAx.XTickLabel(1)={['\color{magenta} perimens']};
hAx.XTickLabel(1)={['\bf ', hAx.XTickLabel{1}]};
hAx.XTickLabel(2)={['\bf ', hAx.XTickLabel{2}]};

subplot(1,3,2)
b = boxplot(corr_val_wsession_bp, 'Color', 'k');
set(b,{'linew'},{1.5})
hold on;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',face_alpha);
end
plot_bp_pnts(corr_val_wsession_bp, 5)
set(gca, 'XTickLabel', [])
set(gca, 'FontSize', 13)
ylim([0 1])
grid minor
title('within session')
set(gca, 'XTickLabel', {'perimens', 'postov', 'preic', 'ic', 'postic', 'interic'})
set(gca, 'XTickLabelRotation', 30)
hAx=gca;
hAx.XAxis.TickLabelInterpreter='tex';      
hAx.XTickLabel(1)={['\bf \color[rgb]{', num2str(c_premenstrual),'} perimens']};
hAx.XTickLabel(2)={['\bf \color[rgb]{', num2str(c_midcycle),'} postov']};
hAx.XTickLabel(3)={['\bf \color[rgb]{', num2str(c_preictal),'} preic']};
hAx.XTickLabel(4)={['\bf \color[rgb]{', num2str(c_ictal),'} ic']};
hAx.XTickLabel(5)={['\bf \color[rgb]{', num2str(c_postictal),'} postic']};
hAx.XTickLabel(6)={['\bf \color[rgb]{', num2str(c_interictal),'} interic']};

colors_aux_2 = [c_preictal; c_ictal; c_postictal; c_interictal];
colors_2 = colors_aux_2(end:-1:1,:);

subplot(1,3,3)
patch('Faces', [1 2 3 4], 'Vertices', [0 y1; 1.5 y1; 1.5 y2; 0 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
hold on;
patch('Faces', [1 2 3 4], 'Vertices', [0 y2; 1.5 y2; 1.5 1; 0 1], 'FaceColor', c_preictal, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [1.5 y1; 2.5 y1; 2.5 y2; 1.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [1.5 y2; 2.5 y2; 2.5 1; 1.5 1], 'FaceColor', c_ictal, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [2.5 y1; 3.5 y1; 3.5 y2; 2.5 y2], 'FaceColor', c_premenstrual, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [2.5 y2; 3.5 y2; 3.5 1; 2.5 1], 'FaceColor', c_postictal, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [3.5 y1; 4.5 y1; 4.5 y2; 3.5 y2], 'FaceColor', c_midcycle, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
patch('Faces', [1 2 3 4], 'Vertices', [3.5 y2; 4.5 y2; 4.5 1; 3.5 1], 'FaceColor', c_interictal, 'FaceAlpha', face_alpha, 'EdgeColor', 'none')
b = boxplot(corr_val_bgroup_bp, 'Color', 'k');
set(b,{'linew'},{1.5})
plot_bp_pnts(corr_val_bgroup_bp, 5)
set(gca, 'XTickLabel', {'preic-perimens', 'ic-perimens', 'postic-perimens', 'interic-postov'})
hAx=gca;
hAx.XAxis.TickLabelInterpreter='tex';         
hAx.XTickLabel(1)={['\bf \color[rgb]{', num2str(c_preictal),'}preic' , ' \color{black}- ',  '\bf \color[rgb]{', num2str(c_premenstrual),'}perimens']};
hAx.XTickLabel(2)={['\bf \color[rgb]{', num2str(c_ictal),'}ic' , ' \color{black}- ',  '\bf \color[rgb]{', num2str(c_premenstrual),'}perimens']};
hAx.XTickLabel(3)={['\bf \color[rgb]{', num2str(c_postictal),'}postic' , ' \color{black}- ',  '\bf \color[rgb]{', num2str(c_premenstrual),'}perimens']};
hAx.XTickLabel(4)={['\bf \color[rgb]{', num2str(c_interictal),'}interic' , ' \color{black}- ',  '\bf \color[rgb]{', num2str(c_midcycle),'}postov']};
set(gca, 'XTickLabelRotation', 30)
set(gca, 'FontSize', 13)
ylim([0 1])
ylabel(' ')
grid minor
title('between groups')
hold off;

print([fig_path, '/DI_correlation_results'],'-dpng');