function plot_carpet_general(data, labels, mo_metric, mo_metric_data, mo_data)
% data - 1 x m cell array for m ROIS
% data{1,R} - v x nr_vols array with the timecourse for ROI "R"


nr_vol = max(size(mo_metric_data));

[irow,~,~] = find(mo_data);

sizedata = zeros(1,size(data,2)+1);
sizedata(1) = 0;
maskdata = [];
for n = 1:size(data,2)
    maskdata = [maskdata; data{1,n}];
    sizedata(n+1) = size(maskdata,1);
    
end

% dvars 
subplot(4,1,1)
yrange = [min(mo_metric_data)-0.5*std(mo_metric_data) max(mo_metric_data)+0.5*std(mo_metric_data)];
plot(mo_metric_data, '-*')
hold on
go = line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r');
title(mo_metric)
xlim([1 nr_vol])
ylim(yrange)
xlabel('Volume number')

subplot(4,1,[2 3 4])
figure;
imagesc(maskdata);
set(gca,'ytick',[]);
colormap(gray);
c = colorbar;
c.LineWidth = 1.5;
c.Location = 'northoutside'; 
xlim([-3 nr_vol]);
title('Carpet plot');
ylabel('Voxel Intensity');
xlabel('Volume Number');

% [0.1 0.9 0.3; 0.4 0.8 0.7; 0.5 0 1]
colors = repmat(lines, 2,1);
for k =1:length(sizedata)-1
    %line([1 nr_vol], [sizedata(k+1) sizedata(k+1)], 'Color', colors(k,:), 'Linewidth', 0.8)
    
    x1=-3; x2=0;
    y1=sizedata(k); y2=sizedata(k+1);
    v=[x1 y1; x1 y2; x2 y2; x2 y1];
    f=[1 2 3 4];
    hi = patch('Faces',f,'Vertices',v,'FaceColor',colors(k,:));
    h(k) = hi(1);
end

legend(h,labels, 'Orientation', 'Horizontal', 'Location', 'southoutside')
