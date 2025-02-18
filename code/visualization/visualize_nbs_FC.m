DI_FC_PCA_path = '/home/iesteves/FC/data/results/DI/DI_FC_PCA';
nbs_stats_path = '/home/iesteves/FC/data/results/nbs/stats/';

% colors
c_preictal =[223, 128, 58]/255;
c_ictal =[240, 0, 32]/255;
c_postictal =[225, 119, 168]/255;
c_interictal =[45, 149, 228]/255;
c_premenstrual=[0, 205, 133]/255;
c_midcycle =[0, 133, 50]/255;

colors_aux = [c_premenstrual; c_midcycle; c_preictal; c_ictal; c_postictal; c_interictal];
colors = colors_aux(end:-1:1,:);

%%
nr_patients = 10;
nr_controls = 14; 
nr_sessions = 6;
ncomp = 19;
contrast_type = '-1';
method = 'Extent';
threshold = 4;
    
load([DI_FC_PCA_path, '/FC_PCA_recon.mat'])
FC = FC_PCA_recon(:,:,:,ncomp);

comparisons = {'preic-interic', 'ic-interic', 'postic-interic', 'prem-mid', 'ic-prem', 'postic-prem', 'preic-ic'};
contrast_types = {'-1', '-1', '-1', '-1', '1', '1', '-1'};

face_alpha = 0.5;

figure('pos', [50 50 1800 1000]);
for c = 1:length(comparisons)
    comparison = comparisons{c};
    contrast_type = contrast_types{c};


    load([nbs_stats_path, 'stats-', num2str(ncomp),'PCs/nbs_', num2str(ncomp), 'PCs_', comparison, '_contrast', contrast_type, '_method-',method, '_cluster-', num2str(round(threshold))]);
    adj_matrix=full(nbs.NBS.con_mat{1}');


    m_val = NaN*ones(nr_controls, nr_sessions);

    for s = 1:6
        if s < 5
            ind_patient = s;
            area = FC(:,:,ind_patient:4:40);
            for k = 1:nr_patients
                area_patient = area(:,:,k);
                area_sig = area_patient(logical(adj_matrix));
                val = mean(area_sig);
                m_val(k,s) = val;      
            end

        else
           ind_control = s+36;
           area = FC(:,:,ind_control:2:68); 
            for k = 1:nr_controls
                area_control = area(:,:,k);
                area_sig = area_control(logical(tril(adj_matrix, 0)));
                val = mean(area_sig);
                m_val(k,s) = val;      
            end
        end
    end


    subplot(2, 4, c)
    b = boxplot([m_val(:,5:6), m_val(:,1:4)], 'Color', 'k');
    set(b,{'linew'},{1.5})
    hold on;
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',face_alpha);
    end
    plot_bp_pnts([m_val(:,5:6), m_val(:,1:4)], 15)
    line([2.5 2.5], [0 1], 'Color', 'k', 'LineStyle', ':', 'Linewidth', 1.5)
    set(gca, 'XTickLabel', [])
    set(gca, 'FontSize', 13)
    ylim([0 1])
    grid minor
    title(comparison)
    set(gca, 'XTickLabel', {'pm', 'postov', 'pre', 'ic', 'postic', 'interic'})
    set(gca, 'XTickLabelRotation', 30)
    hAx=gca;
    hAx.XAxis.TickLabelInterpreter='tex';      
    hAx.XTickLabel(1)={['\bf \color{black}HC-\color[rgb]{', num2str(c_premenstrual),'}pm']};
    hAx.XTickLabel(2)={['\bf \color{black}HC-\color[rgb]{', num2str(c_midcycle),'}postov']};
    hAx.XTickLabel(3)={['\bf \color{black}M-\color[rgb]{', num2str(c_preictal),'}pre']};
    hAx.XTickLabel(4)={['\bf \color{black}M-\color[rgb]{', num2str(c_ictal),'}ict']};
    hAx.XTickLabel(5)={['\bf \color{black}M-\color[rgb]{', num2str(c_postictal),'}post']};
    hAx.XTickLabel(6)={['\bf \color{black}M-\color[rgb]{', num2str(c_interictal),'}inter']};

end



%% Avg FC within and between-network

networks = {'FPN', 'DMN', 'DAN', 'LN', 'VAN', 'SMN', 'VN', 'SUB', 'CRB'};
nr_networks = length(networks);
nr_areas = 130;
areas = [1 13 1; 14 37 2; 38 52 3; 53 57 4; 58 69 5; 70 83 6; 84 100 7; 101 112 8; 113 nr_areas 9];

net_colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0 0 0];


m_val = NaN*ones(nr_networks, nr_networks, nr_controls, nr_sessions);

for s = 1:6
    if s < 5
        ind_patient = s;

        for a = 1:nr_networks
            for b = 1:nr_networks                
                area = FC(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_patient:4:40);              
                for k = 1:nr_patients
                    area_patient = area(:,:,k);
                    if a == b
                        aux = tril(ones(length(area_patient)),-1);
                        area_patient(aux == 0) = NaN;
                    end
                    m_val(a, b, k,s) = nanmean(area_patient(:));                    
                end
            end
        end


    else
       ind_control = s+36;
       for a = 1:nr_networks
            for b = 1:nr_networks                
                area = FC(areas(a,1):areas(a,2), areas(b,1):areas(b,2),ind_control:2:68);              
                for k = 1:nr_controls
                    area_control = area(:,:,k);
                    if a == b
                        aux = tril(ones(length(area_control)),-1);
                        area_control(aux == 0) = NaN;
                    end
                    m_val(a, b, k,s) = nanmean(area_control(:));                    
                end
            end
        end
    end
end

figure('pos', [50 50 1800 1000])
%h = suptitle({[corr_type, ' correlation - Avg network FC vs Cycle Day - No Thresh - ', num2str(nComp), ' PCs'] , '  ', '  '});
%h.Interpreter = 'none';
d = 1;
for a = 1:nr_networks
    for b = 1:nr_networks
        subplot(9,9,d)
        net_FC = squeeze(m_val(a, b,:,:));
        boxplot([net_FC(:,5:6), net_FC(:,1:4)], 'Color', 'k')
        hold on;
        h = findobj(gca,'Tag','Box');
        for j=1:length(h)
            patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',face_alpha);
        end
        plot_bp_pnts([net_FC(:,5:6), net_FC(:,1:4)], 5)
        line([2.5 2.5], [0 1], 'Color', 'k', 'LineStyle', ':', 'Linewidth', 1.5)
        set(gca, 'XTickLabel', [])

        ylim([0 0.9])
        d = d + 1;

        if a~= nr_networks
            set(gca, 'XTickLabel', [])
        end
        if b >1 
            set(gca, 'YTickLabel', [])
        end
        if a == 1 
            title(['\bf \color[rgb]{', num2str(net_colors(b,:)),'}', networks{b}]);
        end
        if b == 1
            ylabel(['\bf \color[rgb]{', num2str(net_colors(a,:)),'}', networks{a}]);
        end
    end
end




