function plot_networks_labels(areas)
% 
% plot_networks_labels - Plots network divisions for different brain
% networks, specified by the areas matrix
%
% INPUTS:
% > areas - A matrix where each row specifies a region with the format:
%           [start_point end_point network_label]
%           The start_point and end_point define the range of the region, and
%           the network_label is an integer that determines the color used for
%           that region.
%

colors_networks = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0.25 0.25 0.25];


    % horizontal
    for a =1:size(areas, 1)

        y1=-10; y2=0.5;
        x1= a-0.5; x2=a+0.5;
        v=[x1 y1; x1 y2; x2 y2; x2 y1];
        f=[1 2 3 4];
        patch('Faces',f,'Vertices',v,'FaceColor',colors_networks(areas(a,3),:), 'EdgeColor', colors_networks(areas(a,3),:));
        hold on;
        line([x2 x2], [0 9.5], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', ':')
    end
    ylim([0 9.5]) 

    % vertical
    for a =1:size(areas, 1)

        x1=-10; x2=0.5;
        y1= a-0.5; y2=a+0.5;
        v=[x1 y1; x1 y2; x2 y2; x2 y1];
        f=[1 2 3 4];
        patch('Faces',f,'Vertices',v,'FaceColor',colors_networks(areas(a,3),:), 'EdgeColor', colors_networks(areas(a,3),:));
        hold on;
        line([0 9.5], [y2 y2], 'Color', 'k', 'LineWidth', 1.4, 'LineStyle', ':')
    end
    ylim([0 9.5])
    xlim([0 9.5])
