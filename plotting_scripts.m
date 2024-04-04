% Load example data
clear;
load('Examples\ExampleOutput.mat');

%% Volume plot example
% plot_volumes(data, pos1, pos2, dir1, dir2, layer, barRange, cMap, dataName, dataUnits)
close all;

% Specific heat capacity with default map and range
figure('WindowState', 'maximize');
plot_volumes(thermMatID, thermX, thermY, 'x', 'y', 3, nan, nan, 'material index', '#')

% Specific heat capacity with default map and range
figure('WindowState', 'maximize');
plot_volumes(thermCp, thermX, thermY, 'x', 'y', 3, nan, nan, 'specific heat capacity', 'J/kg*K')

% Temperature with hot map and specified range
figure('WindowState', 'maximize');
randTemps = 10 + 20*rand(thermDims);
plot_volumes(randTemps, thermX, thermY, 'x', 'y', 3, [10 30], 'hot', 'Temperature', [char(176) 'C'])

%% Line plot example
% plot_lines(data1, data2, pos1, pos2, dir1, dir2, lineWidth, layer, barRange, cMap, dataName, dataUnits)
close all;

% Thermal conductivity with default map and range
figure('WindowState', 'maximize');
plot_lines(thermKx, thermKy, thermX, thermY, 'x', 'y', 1e-3, 4, nan, nan, 'thermal conductivity', 'W/K')

%% Grid plot example
% plot_grid(x,y,z,figTitle,figView)
close all;
% Plot full pack
figure('WindowState', 'maximize');
plot_grid(thermX, thermY, thermZ, 'Pack', '3d'); % view can be, 3d, xy, xz, yz
axis equal;

%% Joined plot example
close all;
figure('WindowState', 'maximize');
range = [min(thermMatID,[],'all') max(thermMatID,[],'all')];
subplot(4,1,1:2);
plot_volumes(thermMatID, thermX, thermY, 'x', 'y', 3, range, nan, 'material index', '#')

subplot(4,1,3);
plot_volumes(thermMatID, thermX, thermZ, 'x', 'z', 2, range, nan, 'material index', '#')

subplot(4,1,4);
plot_volumes(thermMatID, thermY, thermZ, 'y', 'z', 3, range, nan, 'material index', '#')
% do this plot with velocity on the bottom instead, maybe for the report 

%% Plot of material indices showing the battery structure
close all;
f = figure('DefaultAxesFontSize',10,'Position', [100 100 1000 400]); %sgtitle('Radiation linearisation results','FontSize',14)
sgtitle('Structure and materials of a reduced size pack in the xy, xz and yz planes');
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
set(gcf,'color','white');

subplot(2,2,[1 3]);
set(gca, 'OuterPosition', [0, 0.01, 0.5, 0.9]);
plot_volumes(thermMatID, thermX, thermY, 'x', 'y', 3, [1 3], nan, 'material', '#');
title('');%'Reduced pack structure and materials in the xy plane');
colorbar off;
% cb = colorbar;
% cb.Ticks = [1 2 3];
% cb.TickLabels = {'Cell','Aluminium','Coolant'};

subplot(2,2,2);
set(gca, 'OuterPosition', [0.45, 0.52, 0.5, 0.4]);
plot_volumes(thermMatID, thermX, thermZ, 'x', 'z', 2, [1 3], nan, 'material', '#');
title('');%'Reduced pack structure and materials in the xz plane');
colorbar off;
% cb = colorbar;
% cb.Ticks = [1 2 3];
% cb.TickLabels = {'Cell','Aluminium','Coolant'};

subplot(2,2,4);
set(gca, 'OuterPosition', [0.45, 0.12, 0.5, 0.4]);
plot_volumes(thermMatID, thermY, thermZ, 'y', 'z', 2, [1 3], nan, 'material', '#');
title('');%'Reduced pack structure and materials in the yz plane');
colorbar off;


h = axes(f,'visible','off');
set(gca, 'OuterPosition', [0.45, -0.09, 0.54, 0.7]);
cb = colorbar(h,'southoutside');%'Position',[0.95 0.168 0.022 0.7]
cb.Ticks = [0 0.5 1];
cb.TickLabels = {'Cell','Aluminium','Coolant'};

%% Plotting the mapping from volumes to cell
close all;
f = figure('DefaultAxesFontSize',10,'Position', [100 100 500 400]);
sgtitle('Cell positioning for a reduced 4S3P pack in the xy plane');
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
set(gcf,'color','white');
range = [min(thermHeatInds,[],'all') max(thermHeatInds,[],'all')];

plot_volumes(thermHeatInds, thermX, thermY, 'x', 'y', 3, range, nan, 'Cell', '#');
title('');

%% Plotting two example cell stacks
close all;
load('Examples/ExampleOrientationPlotData.mat')

figure('DefaultAxesFontSize',10,'Position', [100 100 800 250]);
sgtitle('Cell orientation default and alternative example (length, width, thickness)')
setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
set(gcf,'color','white');

subplot(1,2,1);
plot_grid(cellStackYZX.x, cellStackYZX.y, cellStackYZX.z, 'Cell stack', '3d'); % view can be, 3d, xy, xz, yz
axis equal;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title('Default orientation (Z, Y, X)');
set(gca, 'OuterPosition', [0.05, 0.05, 0.45, 0.85]);

subplot(1,2,2);
plot_grid(cellStackXYZ.x, cellStackXYZ.y, cellStackXYZ.z, 'Cell stack', '3d'); % view can be, 3d, xy, xz, yz
axis equal;
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
title('Alternative orientation (X, Y, Z)');
set(gca, 'OuterPosition', [0.5, 0.05, 0.45, 0.85]);