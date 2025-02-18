DI_path = '/home/iesteves/FC/data/results/DI';
fig_path = '/home/iesteves/FC/figures/DI_values';
template_path = '/home/iesteves/FC/files/template_matrices';

managefolders(fig_path, 'create');

nr_subjects = 68;

% number of networks (including the whole brain - WB)
networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'Subcortical', 'Cerebellum', 'WB'};
nr_networks = length(networks);

% number of PCs, from 1 to the number of data samples
components = 1:nr_subjects;
nr_components = length(components);

%%
% general (for Idiff)
load([template_path, '/within_group'])
load([template_path, '/within_subject'])
load([template_path, '/within_session'])
load([template_path, '/within_menstrual_session'])
load([template_path, '/between_group_within_session'])

load([DI_path, '/DI_FC_PCA/corr_PCA_recon_networks.mat']);

load([DI_path, '/DI_values/within_session_DI.mat']);
load([DI_path, '/DI_values/within_subject_DI.mat']);
load([DI_path, '/DI_values/within_group_DI.mat']);
load([DI_path, '/DI_values/within_menstrual_session_DI.mat']);
load([DI_path, '/DI_values/between_group_within_session_DI.mat']);

%% % FIGURE Idiff, IA and IB for each template
% get Idiff values for each template and corresponding max/min

line_color =[0 0.4470 0.7410];

for n = 1:nr_networks

    within_session_Idiff = squeeze(within_session_DI(1,n,:));
    within_session_Iself = squeeze(within_session_DI(2,n,:));
    within_session_Iothers  = squeeze(within_session_DI(3,n,:));

    within_subject_Idiff = squeeze(within_subject_DI(1,n,:));
    within_subject_Iself = squeeze(within_subject_DI(2,n,:));
    within_subject_Iothers = squeeze(within_subject_DI(3,n,:));

    within_group_Idiff = squeeze(within_group_DI(1,n,:));
    within_group_Iself = squeeze(within_group_DI(2,n,:));
    within_group_Iothers = squeeze(within_group_DI(3,n,:));

    within_menstrual_session_Idiff = squeeze(within_menstrual_session_DI(1,n,:));
    within_menstrual_session_Iself = squeeze(within_menstrual_session_DI(2,n,:));
    within_menstrual_session_Iothers = squeeze(within_menstrual_session_DI(3,n,:));

    between_group_within_session_Idiff= squeeze(between_group_within_session_DI(1,n,:));
    between_group_within_session_Iself = squeeze(between_group_within_session_DI(2,n,:));
    between_group_within_session_Iothers = squeeze(between_group_within_session_DI(3,n,:));
    
    [y_ma_within_group, ma_within_group] = max(within_group_Idiff);

    [y_ma_within_subject, ma_within_subject] = max(within_subject_Idiff);

    [y_ma_within_session, ma_within_session] = max(within_session_Idiff);

    [y_mi_within_menstrual_session, mi_within_menstrual_session] = min(within_menstrual_session_Idiff);
    [y_ma_within_menstrual_session, ma_within_menstrual_session] = max(within_menstrual_session_Idiff);

    [y_mi_between_group_within_session, mi_between_group_within_session] = min(within_menstrual_session_Idiff);
    [y_ma_between_group_within_session, ma_between_group_within_session] = max(within_menstrual_session_Idiff);

    figure('pos', [200, 1000, 2000, 900]);
    suptitle({[networks{n}, ' - Differential Identifiability: Idiff, Iself and Iothers'], ' ', ' '})
    subplot(2,5,1)
    plot(within_group_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_subjects, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
    xlim([2 nr_components])
    xlabel('# PCs')
    ylabel('Differential Identifiability I_{diff} (%)')
    grid minor
    title('group')
    hold on;
    ylim([0 y_ma_within_group+y_ma_within_group*0.1])
    line([ma_within_group ma_within_group], [0 y_ma_within_group+y_ma_within_group*0.1], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'k')
    text(double(ma_within_group), y_ma_within_group*0.5, {strcat('\leftarrow #PCs =  ',num2str(ma_within_group))}, 'FontSize', 13)
    hold off;
    set(gca, 'FontSize', 13)

    subplot(2,5,6)
    plot(within_group_Iself, 'g', 'LineWidth', 1.2)
    hold on;
    plot(within_group_Iothers, 'r', 'LineWidth', 1.2)
    legend('Iself', 'Iothers')
    xlim([2 nr_components])
    xlabel('# PCs')
    ylabel('Average Pearson Correlation')
    grid minor
    set(gca, 'FontSize', 13)

    subplot(2,5,2)
    plot(within_subject_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_subjects, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    title('subject')
    hold on;
    ylim([0 y_ma_within_subject+y_ma_within_subject*0.1])
    line([ma_within_subject ma_within_subject], [0 y_ma_within_subject+y_ma_within_subject*0.1], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'k')
    text(double(ma_within_subject), y_ma_within_subject*0.5, {strcat('\leftarrow #PCs =  ',num2str(ma_within_subject))}, 'FontSize', 13)
    hold off;
    set(gca, 'FontSize', 13)

    subplot(2,5,7)
    plot(within_subject_Iself, 'g', 'LineWidth', 1.2)
    hold on;
    plot(within_subject_Iothers, 'r', 'LineWidth', 1.2)
    legend('Iself', 'Iothers')
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    set(gca, 'FontSize', 13)

    subplot(2,5,3)
    plot(within_session_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_subjects, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    title('session')
    ylim([0 2])
    hold on;
    line([ma_within_session ma_within_session], [0 2], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'k')
    text(double(ma_within_session), y_ma_within_session*0.5, {strcat('\leftarrow #PCs =  ',num2str(ma_within_session))}, 'FontSize', 13)
    hold off;
    set(gca, 'FontSize', 13)

    subplot(2,5,8)
    plot(within_session_Iself, 'g', 'LineWidth', 1.2)
    hold on;
    plot(within_session_Iothers, 'r', 'LineWidth', 1.2)
    legend('Iself', 'Iothers')
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    set(gca, 'FontSize', 13)

    subplot(2,5,4)
    plot(within_menstrual_session_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    title('cycle')
    ylim([y_mi_within_menstrual_session-y_mi_within_menstrual_session*0.1 y_ma_within_menstrual_session+y_ma_within_menstrual_session*0.1])
    hold on;
    line([mi_within_menstrual_session mi_within_menstrual_session], ...
        [y_mi_within_menstrual_session-y_mi_within_menstrual_session*0.1 y_ma_within_menstrual_session+y_ma_within_menstrual_session*0.1], ...
        'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
    hold off;
    set(gca, 'FontSize', 13)

    subplot(2,5,9)
    plot(within_menstrual_session_Iself, 'g', 'LineWidth', 1.2)
    hold on;
    plot(within_menstrual_session_Iothers, 'r', 'LineWidth', 1.2)
    legend('Iself', 'Iothers')
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    set(gca, 'FontSize', 13)

    subplot(2,5,5)
    plot(between_group_within_session_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    title('group-cycle')
    hold on;
    line([mi_between_group_within_session mi_between_group_within_session], ylim, 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
    hold off;
    set(gca, 'FontSize', 13)

    subplot(2,5,10)
    plot(between_group_within_session_Iself, 'g', 'LineWidth', 1.2)
    hold on;
    plot(between_group_within_session_Iothers, 'r', 'LineWidth', 1.2)
    legend('Iself', 'Iothers')
    xlim([2 nr_components])
    xlabel('# PCs')
    grid minor
    set(gca, 'FontSize', 13)
    
    print([fig_path, '/DI_values_Iself-Iothers_', networks{n}], '-dpng');
end


%% Histograms
for n = 1:nr_networks

    fig = figure('pos', [200, 1000, 2000, 600]);
    suptitle({[networks{n}, ' - Differential Identifiability: Idiff and Correlation histograms: ', '\color{green}Self',' \color{black}and ', '\color{red}Others'], ' ', ' '})
    for m_ind = 2:nr_components

        corr_matrix_original = corr_PCA_recon_networks(:,:,n,nr_components);
        corr_matrix = corr_PCA_recon_networks(:,:,n, m_ind);

        within_session_Idiff = squeeze(within_session_DI(1,n,:));
        within_subject_Idiff = squeeze(within_subject_DI(1,n,:));
        within_group_Idiff = squeeze(within_group_DI(1,n,:));
        within_menstrual_session_Idiff = squeeze(within_menstrual_session_DI(1,n,:));
        between_group_within_session_Idiff = squeeze(between_group_within_session_DI(1,n,:));
        
        [y_ma_within_group,~] = max(within_group_Idiff);
        [y_ma_within_subject, ~] = max(within_subject_Idiff);
        [y_ma_within_session, ~] = max(within_session_Idiff);
     
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
        line([m_ind m_ind], [0 y_ma_within_group], 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([2 nr_components])
        ylim([0 y_ma_within_group])
        xlabel('# PCs')
        ylabel('Idiff')
        title('group')
        hold off;    

        subplot(2,6,8)
        histogram(corr_matrix(within_group==1), 20, 'facecolor', [0 1 0], 'facealpha',.3, 'edgealpha', 0.5)
        hold on;
        histogram(corr_matrix(within_group==0), 'facecolor', [1 0 0], 'facealpha',.3, 'edgealpha', 0.5);
        xlim([0 1])
        hold off;
        ylabel('Counts')
        xlabel('Correlation')

        subplot(2,6,3)
        plot(within_subject_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        title('subject')
        hold on;
        line([m_ind m_ind], [0 y_ma_within_subject], 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([2 nr_components])
        ylim([0 y_ma_within_subject])
        xlabel('# PCs')
        hold off;

        subplot(2,6,9)
        histogram(corr_matrix(within_subject==1), 'facecolor', [0 1 0], 'facealpha',.3, 'edgealpha', 0.5)
        hold on;
        histogram(corr_matrix(within_subject==0), 'facecolor', [1 0 0], 'facealpha',.3, 'edgealpha', 0.5);
        xlim([0 1])
        xlabel('Correlation')
        hold off;

        subplot(2,6,4)
        plot(within_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        line([m_ind m_ind], [0 y_ma_within_session], 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([2 nr_components])
        ylim([0 y_ma_within_session])
        xlabel('# PCs')
        title('session')

        subplot(2,6,10)
        histogram(corr_matrix(within_session==1), 20, 'facecolor', [0 1 0], 'facealpha',.3, 'edgealpha', 0.5)
        hold on;
        histogram(corr_matrix(within_session==0), 'facecolor', [1 0 0], 'facealpha',.3, 'edgealpha', 0.5);
        xlim([0 1])
        xlabel('Correlation')
        hold off;

        subplot(2,6,5)
        plot(within_menstrual_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        hold on;
        line([m_ind m_ind], [-1.5 0.5], 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        xlim([2 nr_components])
        ylim([-1.5 0.5])
        title('cycle')
        xlabel('# PCs')
        hold off;

        subplot(2,6,11)
        histogram(corr_matrix(within_menstrual_session==1), 20, 'facecolor', [0 1 0], 'facealpha',.3, 'edgealpha', 0.5)
        hold on;
        histogram(corr_matrix(within_menstrual_session==0), 'facecolor', [1 0 0], 'facealpha',.3, 'edgealpha', 0.5);
        xlim([0 1])
        xlabel('Correlation')
        hold off;

        subplot(2,6,6)
        plot(between_group_within_session_Idiff, '-o','MarkerSize',3, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[0 0 0], 'MarkerIndices', m_ind, 'LineWidth', 1.2)
        xlim([2 nr_components])
        xlabel('# PCs')
        hold on;
        line([m_ind m_ind], ylim, 'Color', 'k', 'LineWidth', 1.2);
        grid minor
        title('group-cycle')
        hold off;

        subplot(2,6,12)
        histogram(corr_matrix(between_group_within_session==1), 20, 'facecolor', [0 1 0], 'facealpha',.3, 'edgealpha', 0.5)
        hold on;
        histogram(corr_matrix(between_group_within_session==0), 'facecolor', [1 0 0], 'facealpha',.3, 'edgealpha', 0.5);
        xlim([0 1])
        xlabel('Correlation')
        hold off;

        %pause(0.1)
        drawnow limitrate
        frame = getframe(fig);
        im{m_ind-1} = frame2im(frame);
    end
    
    filename = [fig_path, '/DI_values_histograms-gif_', networks{n},'.gif']; % Specify the output file name
    for idx = 1:nr_components-1
        [A,map] = rgb2ind(im{idx},256);
        if idx == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
        end
    end
end


%% FIGURE Brain Function Workshop - panel A
line_color =[0 0.4470 0.7410];
comp = 68; % reconstruction with 68 PCs to show original matrix
nr_pcs = 19; % optimal number of PCs
n = 10; % whole-brain

original_corr = corr_PCA_recon_networks(:,:,n,comp);

within_session_Idiff = squeeze(within_session_DI(1,n,:));
within_subject_Idiff = squeeze(within_subject_DI(1,n,:));
within_group_Idiff = squeeze(within_group_DI(1,n,:));

[y_ma_within_group, ma_within_group] = max(within_group_Idiff);
[y_ma_within_subject, ma_within_subject] = max(within_subject_Idiff);
[y_ma_within_session, ma_within_session] = max(within_session_Idiff);

% FIGURE Idiff, IA and IB for each template
figure('pos', [200, 1000, 300, 1000]);
subplot(4,1,1)
imagesc(original_corr);
colorbar
title({'Identifiability','(Original)'})
set(gca, 'FontSize', 13)
ylabel('Data Samples')
xlabel('Data Samples')

subplot(4,1,2)
plot(within_subject_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
%xlabel('# PCs')
grid minor
title('subject')
hold on;
ylim([0 25])
line([ma_within_subject ma_within_subject], [0 25], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 25], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_subject), y_ma_within_subject*0.5, {strcat('\leftarrow #PCs =  ',num2str(ma_within_subject))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_subject*0.7, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

subplot(4,1,3)
plot(within_session_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
%xlabel('# PCs')
ylabel('Differential Identifiability I_{diff} (%)')
grid minor
title('session')
ylim([0 3])
hold on;
line([ma_within_session ma_within_session], [0 3], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 3], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_session), y_ma_within_session*1.2, {strcat('\leftarrow #PCs =  ',num2str(ma_within_session))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_session*1.6, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

subplot(4,1,4)
plot(within_group_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
xlabel('# PCs')
grid minor
title('group')
hold on;
ylim([0 3])
line([ma_within_group ma_within_group], [0 3], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 3], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_group), y_ma_within_group*0.4, {strcat('\leftarrow #PCs =  ',num2str(ma_within_group))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_group*0.6, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

print([fig_path, '/DI_values_BrainFunctionWorkshop-panelA_Idiff'],'-dpng');

%% FIGURE panel A - modified
line_color =[0 0.4470 0.7410];
comp = 68; % reconstruction with 68 PCs to show original matrix
nr_pcs = 19; % optimal number of PCs
n = 10; % whole-brain

original_corr = corr_PCA_recon_networks(:,:,n,comp);


within_session_Idiff = squeeze(within_session_DI(1,n,:));
within_subject_Idiff = squeeze(within_subject_DI(1,n,:));
within_group_Idiff = squeeze(within_group_DI(1,n,:));

[y_ma_within_group, ma_within_group] = max(within_group_Idiff);
[y_ma_within_subject, ma_within_subject] = max(within_subject_Idiff);
[y_ma_within_session, ma_within_session] = max(within_session_Idiff);

% FIGURE Idiff, IA and IB for each template

figure('pos', [200, 1000, 1600, 500]);
subplot(2,4,8)
imagesc(corr_PCA_recon_networks(:,:,n,nr_pcs));
colorbar
title({['Identifiability (19 PCs)']})
set(gca, 'FontSize', 13)
ylabel('[Subjects x Sessions]')
xlabel('[Subjects x Sessions]')

subplot(2,4,[1 5])
plot(within_subject_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
%xlabel('# PCs')
grid minor
title('subject')
hold on;
ylim([0 25])
line([ma_within_subject ma_within_subject], [0 25], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 25], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_subject), y_ma_within_subject*0.5, {strcat('\leftarrow #PCs =  ',num2str(ma_within_subject))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_subject*0.7, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

subplot(2,4,[2 6])
plot(within_session_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
%xlabel('# PCs')
ylabel('Differential Identifiability I_{diff} (%)')
grid minor
title('session')
ylim([0 3])
hold on;
line([ma_within_session ma_within_session], [0 3], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 3], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_session), y_ma_within_session*1.2, {strcat('\leftarrow #PCs =  ',num2str(ma_within_session))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_session*1.6, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

subplot(2,4,[3 7])
plot(within_group_Idiff, '-o','MarkerSize',2, 'MarkerEdgeColor', line_color ,'MarkerFaceColor',line_color, 'MarkerIndices', 1:2:nr_components, 'Color', line_color, 'LineWidth', 1.2, 'LineWidth', 1.2)
xlim([2 nr_components])
xlabel('# PCs')
grid minor
title('group')
hold on;
ylim([0 3])
line([ma_within_group ma_within_group], [0 3], 'LineStyle', ':', 'LineWidth', 1.1, 'Color', 'k')
line([nr_pcs nr_pcs], [0 3], 'LineStyle', '-.', 'LineWidth', 1.1, 'Color', 'r')
text(double(ma_within_group), y_ma_within_group*0.4, {strcat('\leftarrow #PCs =  ',num2str(ma_within_group))}, 'FontSize', 13)
text(nr_pcs, y_ma_within_group*0.6, {strcat('\leftarrow #PCs =  ',num2str(nr_pcs))}, 'FontSize', 13, 'Color', 'r')
hold off;
set(gca, 'FontSize', 13)

subplot(2,4,4)
imagesc(original_corr);
colorbar
title({['Identifiability','(Original)']})
set(gca, 'FontSize', 13)
ylabel('[Subjects x Sessions]')
xlabel('[Subjects x Sessions]')

print([fig_path, '/DI_values_panelA-modified_Idiff'],'-dpng');