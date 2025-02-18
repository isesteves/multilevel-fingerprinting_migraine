function plot_atlas_labels(atlas, xmax)
%
%   plot_atlas_labels generates a plot of brain regions
%   based on the specified atlas, by adding colored rectangles on the sides. The plot is scaled according to the 
%   maximum X value provided by xmax.
%
%   INPUTS:
%   > atlas - A string specifying the brain atlas to use. Valid options are:
%           'HO' - Harvard-Oxford atlas.
%           'AAL116' - Automated Anatomical Labeling atlas with 116 regions.
%           'Desikan' - Desikan-Killiany atlas.
%           'SchaeferSubCRB7100' - Schaefer atlas with 100 regions and subcortical/cerebellar areas.
%   > xmax - A numeric value specifying the maximum X value for scaling the plot.
%

switch atlas
    
    case 'HO'
        areas = [1 15 1; 31 37 2; 16 16 3; 19 20 3; 38 52 3; 23 30 4; 17 18 5; 21 22 5; 53 63 5];

    case  'AAL116'
        areas = [1 32 1; 33 34 3; 35 40 1; 41 42 5; 43 44 3; 45 46 5; ...
            47 60 4; 61 74 2; 75 82 5; 83 94 3; 95 116 6];

    case 'Desikan'
        areas = [1 26 1; 27 30 3; 31 38 4; 39 52 2; 53 66 3];
        
    case 'SchaeferSubCRB7100'
        areas = [1 13.5 1; 13.5 37.5 2; 37.5 52.5 3; 52.5 57.5 4; 57.5 69.5 5; 69.5 83.5 6; 83.5 100.5 7; 100.5 112.5 8; 112.5 130.5 9];
    
    otherwise
        disp('Unrecognized, please choose a different atlas')
end

regions = {'Frontal', 'Parietal', 'Temporal', 'Occipital', 'Subcortical', 'Cerebellar'};
networks = {'Cont', 'Default', 'DorsAttn', 'Limbic', 'SalVentAttn', 'SomMot', 'Vis', 'Subcort', 'Cereb'};
colors = [0.6350 0.0780 0.1840; 0.9290 0.6940 0.1250; 0.4660 0.6740 0.1880; ...
    0 0.4470 0.7410; 0.6 0.6 0.6; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330; 0.8500 0.3250 0.0980; 0.25 0.25 0.25];

% vertical
for a =1:size(areas, 1)
  
    x1=-10; x2=1;
    y1=areas(a, 1); y2=areas(a, 2);
    v=[x1 y1; x1 y2; x2 y2; x2 y1];
    f=[1 2 3 4];
    hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
    h(a) = hi(1);
end
xlim([-xmax*0.05 xmax])

% horizontal
for a =1:size(areas, 1)
  
    y1=-10; y2=1;
    x1=areas(a, 1); x2=areas(a, 2);
    v=[x1 y1; x1 y2; x2 y2; x2 y1];
    f=[1 2 3 4];
    hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(areas(a,3),:), 'EdgeColor', colors(areas(a,3),:));
    h(a) = hi(1);
end
ylim([-xmax*0.05 xmax])

% lines
val = 0;
for a=1:size(areas, 1)
    line([areas(a,1)-val areas(a,2)-val], [areas(a, 1)-val, areas(a, 1)-val], 'Color', colors(areas(a, 3), :), 'LineWidth', 1.2) % top
    line([areas(a,1)-val areas(a,2)-val], [areas(a, 2)-val, areas(a, 2)-val], 'Color', colors(areas(a, 3), :), 'LineWidth', 1.2) % bottom
    line([areas(a,1)-val areas(a,1)-val], [areas(a, 1)-val, areas(a, 2)-val], 'Color', colors(areas(a, 3), :), 'LineWidth', 1.2) % left
    line([areas(a,2)-val areas(a,2)-val], [areas(a, 1)-val, areas(a, 2)-val], 'Color', colors(areas(a, 3), :), 'LineWidth', 1.2) % right    
end

% if any(strcmp(atlas, {'HO', 'AAL116', 'Desikan'}))
%     L1 = plot(nan, nan, 'color', colors(1,:),'LineWidth', 2);
%     L2 = plot(nan, nan, 'color', colors(2,:),'LineWidth', 2);
%     L3 = plot(nan, nan, 'color', colors(3,:),'LineWidth', 2);
%     L4 = plot(nan, nan, 'color', colors(4,:),'LineWidth', 2);
%     L5 = plot(nan, nan, 'color', colors(5,:),'LineWidth', 2);
%     L6 = plot(nan, nan, 'color', colors(6,:),'LineWidth', 2);
%     hl = legend([L1, L2, L3, L4, L5, L6], regions, 'Orientation', 'Horizontal', 'Location', 'southoutside');
%     hl.Position = [0.5 0.01 0.01 0.01]; 
% else
%     L1 = plot(nan, nan, 'color', colors(1,:),'LineWidth', 2);
%     L2 = plot(nan, nan, 'color', colors(2,:),'LineWidth', 2);
%     L3 = plot(nan, nan, 'color', colors(3,:),'LineWidth', 2);
%     L4 = plot(nan, nan, 'color', colors(4,:),'LineWidth', 2);
%     L5 = plot(nan, nan, 'color', colors(5,:),'LineWidth', 2);
%     L6 = plot(nan, nan, 'color', colors(6,:),'LineWidth', 2);
%     L7 = plot(nan, nan, 'color', colors(7,:),'LineWidth', 2);
%     L8 = plot(nan, nan, 'color', colors(8,:),'LineWidth', 2);
%     L9 = plot(nan, nan, 'color', colors(9,:),'LineWidth', 2);
%     hl = legend([L1, L2, L3, L4, L5, L6, L7, L8, L9], networks, 'Orientation', 'Horizontal', 'Location', 'southoutside');
%     hl.Position = [0.5 0.01 0.01 0.01]; 
% end