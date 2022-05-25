function plot_grid(x,y,z,figTitle,figView)
    % Plots a grid with a single line with x,y,z values from FVM
    [X1, Y1, Z1] = meshgrid(x([1 end]),y,z); % Meshgrid with start and end x values
    X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = permute(Z1,[2 1 3]); % Rearrange the dimensions to correct order
    X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN; % Stops connection back to the start

    [X2, Y2, Z2] = meshgrid(x,y([1 end]),z); % Meshgrid with start and end y values - do not need to rearrange 
    X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN; % Stops connection back to the start

    [X3, Y3, Z3] = meshgrid(x,y,z([1 end]));% Meshgrid with start and end z values
    X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = permute(Z3,[3 1 2]); % Rearrange the dimensions to correct order
    X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN; % Stops connection back to the start
  
    line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title(figTitle);
    switch figView % Rotate cell stack
        case '3d'
            view(3)
        case 'xy'
            view(0,90)
        case 'xz'
            view(0,0)
        case 'yz'
            view(90,0)
        otherwise
            disp('Invalid view input - leaving default')
    end
end

