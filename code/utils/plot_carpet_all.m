function plot_carpet_all(gm_data, wm_data, csf_data, resp_data, resp_times, mo_metric, mo_metric_data, mo_data, TR)


nr_vol = max(size(mo_metric_data));

[irow,~,~] = find(mo_data);

alldata = double([gm_data(:); wm_data(:); csf_data(:)]);
clims = prctile(alldata, [1 99]);

% resp
subplot(5,1,1)
plot(resp_times, resp_data)
xlim([1 nr_vol*TR])
xlabel('Time (s)')
title('resp')

% dvars 
subplot(5,1,2)
yrange = [min(mo_metric_data)-0.5*std(mo_metric_data) max(mo_metric_data)+0.5*std(mo_metric_data)];
plot(mo_metric_data, '-*')
hold on
go = line(repmat(irow', 2, 1), repmat(yrange, length(irow), 1)', 'Color', 'r');
title(mo_metric)
xlim([1 nr_vol])
ylim(yrange)
xlabel('Volume number')

subplot(5,1,[3 4 5])
cim = imagesc([gm_data; wm_data; csf_data]);
set(gca, 'CLim', clims)
set(gca,'ytick',[]);
colormap(gray);
c = colorbar;
c.LineWidth = 1.5;
c.Location = 'northoutside';
line([1 nr_vol], [size(gm_data,1) size(gm_data,1)], 'Color','g', 'Linewidth', 0.8);
line([1 nr_vol], [size([gm_data; wm_data],1) size([gm_data; wm_data],1)], 'Color', 'b', 'Linewidth', 0.8);
line([1 nr_vol], [size([gm_data; wm_data; csf_data],1) size([gm_data; wm_data; csf_data],1)], 'Color', 'r', 'Linewidth', 0.8);

xlim([-3 nr_vol]);
title('Carpet plot');
ylabel('Voxel Intensity');
xlabel('Volume Number');

x1=-3; x2=0;
y1=0; y2=size(gm_data,1);
v=[x1 y1; x1 y2; x2 y2; x2 y1];
f=[1 2 3 4];
g1 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.9 0.3]);

x1=-3; x2=0;
y1=size(gm_data,1); y2=size([gm_data; wm_data],1);
v=[x1 y1; x1 y2; x2 y2; x2 y1];
f=[1 2 3 4];
g2 = patch('Faces',f,'Vertices',v,'FaceColor',[0.1 0.1 1]);

x1=-3; x2=0;
y1=size([gm_data; wm_data],1); y2=size([gm_data; wm_data; csf_data],1);
v=[x1 y1; x1 y2; x2 y2; x2 y1];
f=[1 2 3 4];
g3 = patch('Faces',f,'Vertices',v,'FaceColor',[1 0.1 0.1]);

g = [g1(1), g2(1), g3(1)];
legend(g, 'GM', 'WM', 'CSF', 'Orientation', 'Horizontal', 'Location', 'southoutside')
