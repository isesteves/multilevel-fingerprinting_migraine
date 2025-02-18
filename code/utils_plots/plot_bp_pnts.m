function plot_bp_pnts(data, mksize)


    x1=ones(size(data,1)).*(1+(rand(size(data,1))-0.5)/30);
    
    for k = 1:size(data,2)
        f1=scatter(x1(:,1)*k,data(:,k),mksize, 'k','filled');f1.MarkerFaceAlpha = 0.4;hold on
    end