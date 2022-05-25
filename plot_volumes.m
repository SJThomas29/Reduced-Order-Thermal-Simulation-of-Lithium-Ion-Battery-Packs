function plot_volumes(data, pos1, pos2, dir1, dir2, layer, barRange, cMap, dataName, dataUnits)
% Will plot the grid volumes at a colour based on data
% data - matrix with data for each volume
% pos1/2 - position data in first/second dimension
% dir1/2 - directions of plot, defines the plane
% layer - layer in the volume to be plotted
% barRange - 2 element vector with upper and lower limits for colour bar
% cmap - colourmap name, or leave as nan for default 
% dataName - name of the data being plotted for labels
% dataUnits - units of the data being plotted for axis

% Define position structure
pos.(dir1) = pos1; pos.(dir2) = pos2;

% Define dimensions
dims.(dir1) = length(pos.(dir1)) - 1; dims.(dir2) = length(pos.(dir2)) - 1;

% Find third dimension
if ~contains([dir1 dir2],'x')
    dir3 = 'x'; 
elseif ~contains([dir1 dir2],'y')
    dir3 = 'y'; 
else
    dir3 = 'z';
end

% Assign axis
xlabel([dir1 ' [m]']); ylabel([dir2 ' [m]']);

% Assign title
title(['Volume ' dataName ' shown in the ' dir1 dir2 ' plane at ' dir3 ' layer ' num2str(layer)]);

% Loops through dir1 and dir2 to plot lines
points1 = zeros(4,dims.(dir1)*dims.(dir2)); points2 = points1;
volColour = zeros(width(points1),1);
dataInd.(dir3) = layer;
n = 1; % Initialise counter
for ind1 = 1:dims.(dir1) % Loop through dir1 regions
    dataInd.(dir1) = ind1;
    for ind2 = 1:dims.(dir2) % Loop through dir2 regions
        dataInd.(dir2) = ind2;
        points1(:,n) = [pos.(dir1)(ind1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1)]';
        points2(:,n) = [pos.(dir2)(ind2) pos.(dir2)(ind2) pos.(dir2)(ind2 + 1) pos.(dir2)(ind2 + 1)]';
        volColour(n,1) = data(dataInd.x,dataInd.y,dataInd.z);
        n = n + 1;
    end
end
% Plot patches and apply colour map
patch(points1,points2,volColour);
if ~isnan(cMap)
    if strcmp(cMap,'thermal')
        colormap(cmocean('Thermal'));
    else
        colormap(cMap);
    end
else
    colormap();
end
bar = colorbar; bar.Label.String = [dataName ' [' dataUnits ']'];
bar.Label.FontSize = 12;
if ~isnan(barRange)
    caxis(barRange);
end
axis equal; 
axis tight;
end

