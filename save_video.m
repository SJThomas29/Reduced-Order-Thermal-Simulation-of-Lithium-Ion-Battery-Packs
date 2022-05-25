% Makes a video showing progression of model over time
close all;
% Get inputs for the video preferences - can use script inputs
if ~exist('script', 'var') % Skip this part if called by a script
%     input = ["v" "y" "volTemp" "temperature" "K" "x" "y" "3" "3e-3" "hot" "300" "15" "-273.15" "y"];
    input = string(inputdlg({'Plot type - volumes or line [v/l]' 'Include velocity? [y/n]' 'Data name - from out [""]'...
        'Data name for label [""]' 'Data units for label [""]' 'Plot direction 1 [x,y,z]' 'Plot direction 2 [x,y,z]'...
        'Layer? [#]' 'Line width [m]' 'Colormap [""]' 'Frames [#]' 'Framerate [#]' 'Include coolant? [y/n]'},...
        'Video config', [1 50], {'v' 'y' 'volTemp' 'temperature' [char (176) 'C'] 'x' 'y' '4' '3e-3' 'thermal' '450' '15' 'y'}));      
end

% Get data using input
data = out.(input(3)).Data;
time = out.(input(3)).Time;
if input(2) == "y"
   vel = out.velocity.Data; 
end
if input(13) == "y"
    coolOut = out.coolOutTemp.Data;
    coolIn = out.coolInTemp.Data;
    cellTemp = out.cellTemp.Data;
    meanTemp = squeeze(mean(cellTemp,[1 2]));
end

% Define position structure
pos.x = thermX; pos.y = thermY; pos.z = thermZ;

% Define values from inputs
dataLabel = char(input(4)); unitLabel = char(input(5));
dir1 = char(input(6)); dir2 = char(input(7));
layer = str2double(input(8)); 
lineWidth = str2double(input(9)); 
cMap = char(input(10));
frames = str2double(input(11)); frameRate = str2double(input(12)); 

% Define figure
figure('WindowState', 'maximize');

% Define frames for video and loop through
vidInds = round(linspace(1, length(time), frames));
cMin = min(data,[],'all'); cMax = max(data,[],'all'); % Min and max values for colorbar
barRange = [cMin cMax];
for frame = 1:frames % Step through each time value
    % Display frame progress
    disp(['Frame ' num2str(frame) ' of ' num2str(frames)]);
    
    % Define title
%     if input(14) == "y"; fullTitle = [fullTitle ', Coolant output temp: ' num2str(round(coolOut(step),1)) ' [' char(176) 'C]']; end
    
    % Plot data
    step = vidInds(frame);
    fullTitle = ['Time: ' num2str(time(step)) 's'];
    if input(2) == "y" % Include velocity subplot
        subplot(3,1,3);
        plot(time(1:step),vel(1:step));%,'k-');
        xlabel('time [s]'); 
        ylabel('velocity [m/s]');
        xlim([time(1) time(end)]);
        ylim([0 max(vel)+5]);
        
        if input(13) == "y" % Include coolant on other axis
            yyaxis right;
            ylabel(['temperature [' char(176) 'C]']);
            plot(time(1:step),meanTemp(1:step));
            hold on;
            plot(time(1:step),coolIn(1:step));%,'b-');
            plot(time(1:step),coolOut(1:step));%,'r-');
            grid on;
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
    
%     if input(1) == "v" % Volume plot
        plot_volumes(data(:,:,:,step), pos.(dir1), pos.(dir2), dir1, dir2, layer, barRange, cMap, dataLabel, unitLabel);
%     else
%         
%         
%     end
    
    if input(2) ~= "y" % append time to single title if needed 
        figTitle = get(gca,'title');
        titleText = get(figTitle,'String');
        set(figTitle,'String',[titleText ' - ' fullTitle]);
    end
    
    % Add title
        
    % Assign frame    
    F(frame) = getframe(gcf);
%     cla; % Clear axis for next frame
end

% Define video
video = VideoWriter('ModelVid','Motion JPEG AVI');
video.FrameRate = frameRate;
open(video);
writeVideo(video,F);
close(video);
disp('Video saved');

close all;
