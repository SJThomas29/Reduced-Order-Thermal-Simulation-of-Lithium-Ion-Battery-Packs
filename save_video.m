% Makes a video showing progression of model over time
close all;
% Get inputs for the video preferences - can use script inputs
if ~exist('script', 'var') % Skip this part if called by a script
%     input = ["y" "volTemp" "temperature" "K" "x" "y" "3" "3e-3" "hot" "300" "15" "-273.15" "y"];
    input = string(inputdlg({'Include velocity? [y/n]' 'Data name - from out [""]'...
        'Data name for label [""]' 'Data units for label [""]' 'Plot direction 1 [x,y,z]' 'Plot direction 2 [x,y,z]'...
        'Layer? [#]' 'Line width [m]' 'Colormap [""] - Thermal uses cmocean (https://uk.mathworks.com/matlabcentral/fileexchange/57773-cmocean-perceptually-uniform-colormaps)' 'Frames [#]' 'Framerate [#]' 'Include coolant? [y/n]'},...
        'Video config', [1 50], {'y' 'volTemp' 'temperature' [char (176) 'C'] 'x' 'y' '4' '3e-3' 'thermal' '450' '15' 'y'}));      
end

% Get data using input
data = out.(input(2)).Data;
time = out.(input(2)).Time;
if input(1) == "y"
   vel = out.velocity.Data; 
end
if input(12) == "y"
    coolOut = out.coolOutTemp.Data;
    coolIn = out.coolInTemp.Data;
    cellTemp = out.cellTemp.Data;
    meanTemp = squeeze(mean(cellTemp,[1 2]));
end

% Define position structure
pos.x = thermX; pos.y = thermY; pos.z = thermZ;

% Define values from inputs
dataLabel = char(input(3)); unitLabel = char(input(4));
dir1 = char(input(5)); dir2 = char(input(6));
layer = str2double(input(7)); 
lineWidth = str2double(input(8)); 
cMap = char(input(9));
frames = str2double(input(10)); frameRate = str2double(input(11)); 

% Define figure
figure('WindowState', 'maximize');

% Define frames for video and loop through
vidInds = round(linspace(1, length(time), frames));
cMin = min(data,[],'all'); cMax = max(data,[],'all'); % Min and max values for colorbar
barRange = [cMin cMax];
for frame = 1:frames % Step through each time value
    % Display frame progress
    disp(['Frame ' num2str(frame) ' of ' num2str(frames)]);
    
    % Plot data
    step = vidInds(frame);
    fullTitle = ['Time: ' num2str(time(step)) 's'];
    if input(1) == "y" % Include velocity subplot
        subplot(3,1,3);
        plot(time(1:step),vel(1:step));%,'k-');
        xlabel('time [s]'); 
        ylabel('velocity [m/s]');
        xlim([time(1) time(end)]);
        ylim([0 max(vel)+5]);
        grid on;

        if input(12) == "y" && (topCool || bottomCool)  % Include coolant on other axis, if used in model
            yyaxis right;
            plot(time(1:step),meanTemp(1:step));
            hold on;
            plot(time(1:step),coolIn(1:step),'b-');
            plot(time(1:step),coolOut(1:step),'m-');
            ylabel(['temperature [' char(176) 'C]']);
            hold off;
            legend('Vehicle velocity','Mean cell temp','Coolant in temp','Coolant out temp');
            title('Graph of vehicle velocity. coolant and cell temperature');
            ylim([floor(min([min(coolIn) min(coolOut)])) - 1 ceil(max([max(coolIn) max(coolOut)])) + 1]);
            yyaxis left;
        else
            title('Graph of vehicle velocity');
        end
        subplot(3,1,[1 2]); 
    end
    
    % Volume plot
    plot_volumes(data(:,:,:,step), pos.(dir1), pos.(dir2), dir1, dir2, layer, barRange, cMap, dataLabel, unitLabel);
    
    if input(1) ~= "y" % append time to single title if needed 
        figTitle = get(gca,'title');
        titleText = get(figTitle,'String');
        set(figTitle,'String',[titleText ' - ' fullTitle]);
    end
            
    % Assign frame
    if frame == 1
        F = getframe(gcf);
    else
        F(frame) = getframe(gcf);
    end
end

% Define video
video = VideoWriter('ModelVid','Motion JPEG AVI');
video.FrameRate = frameRate;
open(video);
writeVideo(video,F);
close(video);
disp('Video saved');

close all;
