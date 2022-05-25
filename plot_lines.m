function plot_lines(data1, data2, pos1, pos2, dir1, dir2, lineWidth, layer, barRange, cMap, dataName, dataUnits)
% Will plot the grid lines at a colour based on data
% data1/2 - line data in first/second dimension
% pos1/2 - position data in first/second dimension
% dir1/2 - directions of plot, defines the plane
% lineWidth - width of plotted lines
% layer - layer in the volume to be plotted
% barRange - 2 element vector with upper and lower limits for colour bar
% cmap - colourmap name, or leave as nan for default 
% dataName - name of the data being plotted for labels
% dataUnits - units of the data being plotted for axis

% Define data structure
data.(dir1) = data1; data.(dir2) = data2;

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

% Define and plot outline for clarity
outline1 = [(pos.(dir1)(1) - 10*lineWidth) (pos.(dir1)(end) + 10*lineWidth) (pos.(dir1)(end) + 10*lineWidth) (pos.(dir1)(1) - 10*lineWidth)]';
outline2 = [(pos.(dir2)(1) - 10*lineWidth) (pos.(dir2)(1) - 10*lineWidth) (pos.(dir2)(end) + 10*lineWidth) (pos.(dir2)(end) + 10*lineWidth)]';
outlineC = [0 0 0];
patch(outline1, outline2, outlineC,'EdgeColor','none');

% Assign axis
xlabel([dir1 ' [m]']); ylabel([dir2 ' [m]']);

% Assign title
title(['Boundary ' dataName ' shown in the ' dir1 dir2 ' plane at ' dir3 ' layer ' num2str(layer)]);

% Loops through dir1 and dir2 to plot lines
points1 = zeros(4,dims.(dir1)*dims.(dir2)*4); points2 = points1;
volColour = zeros(width(points1),1);
dataInd.(dir3) = layer;
n = 1; % Initialise counter
for ind1 = 1:dims.(dir1) % Loop through dir1 regions
    dataInd.(dir1) = ind1;
    for ind2 = 1:dims.(dir2) % Loop through dir2 regions
        dataInd.(dir2) = ind2;
        % Lines in the dir1 direction will use transfer/conductivity in the
        % dir2 direction
        % First line in dir1 direction - using a rectangle for colour
        points1(:,n) = [pos.(dir1)(ind1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1)]';
        points2(:,n) = [(pos.(dir2)(ind2) - lineWidth/2) (pos.(dir2)(ind2) - lineWidth/2) (pos.(dir2)(ind2) + lineWidth/2)  (pos.(dir2)(ind2) + lineWidth/2)]';
        volColour(n,1) = data.(dir2)(dataInd.x,dataInd.y,dataInd.z);
        n = n + 1;
        % First line in dir2 direction - using a rectangle for colour
        points1(:,n) = [(pos.(dir1)(ind1) - lineWidth/2) (pos.(dir1)(ind1) - lineWidth/2) (pos.(dir1)(ind1) + lineWidth/2)  (pos.(dir1)(ind1) + lineWidth/2)]';
        points2(:,n) = [pos.(dir2)(ind2) pos.(dir2)(ind2 + 1) pos.(dir2)(ind2 + 1) pos.(dir2)(ind2)]'';
        volColour(n,1) = data.(dir1)(dataInd.x,dataInd.y,dataInd.z);
        n = n + 1;
        
        if ind1 == dims.(dir1) % Last in dir1 need the end line
            dataInd.(dir1) = dataInd.(dir1) + 1;
            points1(:,n) = [(pos.(dir1)(ind1 + 1) - lineWidth/2) (pos.(dir1)(ind1 + 1) - lineWidth/2) (pos.(dir1)(ind1 + 1) + lineWidth/2)  (pos.(dir1)(ind1 + 1) + lineWidth/2)]';
            points2(:,n) = [pos.(dir2)(ind2) pos.(dir2)(ind2 + 1) pos.(dir2)(ind2 + 1) pos.(dir2)(ind2)]'';
            volColour(n,1) = data.(dir1)(dataInd.x,dataInd.y,dataInd.z);
            n = n + 1;
            dataInd.(dir1) = dataInd.(dir1) - 1;            
        end
        
        if ind2 == dims.(dir2) % Last in dir2 need the end line
            dataInd.(dir2) = dataInd.(dir2) + 1;
            points1(:,n) = [pos.(dir1)(ind1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1 + 1) pos.(dir1)(ind1)]';
            points2(:,n) = [(pos.(dir2)(ind2 + 1) - lineWidth/2) (pos.(dir2)(ind2 + 1) - lineWidth/2) (pos.(dir2)(ind2 + 1) + lineWidth/2)  (pos.(dir2)(ind2 + 1) + lineWidth/2)]';
            volColour(n,1) = data.(dir2)(dataInd.x,dataInd.y,dataInd.z);
            n = n + 1;
            dataInd.(dir2) = dataInd.(dir2) - 1;            
        end
    end
end
% Plot patches and apply colour map
patch(points1,points2,volColour,'EdgeColor','none');
if ~isnan(cMap)
    colormap(cMap);
else
    colormap();
end
bar = colorbar; bar.Label.String = [dataName ' [' dataUnits ']'];
if ~isnan(barRange)
    caxis(barRange);
end
axis equal; 
axis tight;
end

