function [] = animateGraph(x, y, t, u, 'Parent', axesHandle2)

%Function to animate the 2d graphs to plot a contour plot of width of tile
%against tile thickness against temperature in that section of the tile
%with varying time
%
% Input Arguments:
% x - Tile thickness (m)
% y - Width of tile (m)
% t - time(s)
% u - temperature at the section of the tile (C)

    [X, Y] = meshgrid(x, y);

    % Create a figure
    figure;

    % Loop for animation
    for i = 1:length(t)
        % Update the surface plot
        if i == 1
            U = u(:,:,i);
            % If first iteration, create the surface plot
            surf(Y,X,U, 'Parent', axesHandle2);  
            xlim([0 0.05])   % Limit graph axes
            colormap('jet');  
            colorbar;  
            xlabel('Thickness of Tile (m)');
            ylabel('Width of Tile (m)');
            zlabel('Temperature(C)');
            title(sprintf('Surface Plot at t = %.2f', t(i)));  % Example: set title with current time
            % Customize the plot appearance as needed
        else
            U = u(:,:,i);
            % For subsequent iterations, update the surface plot
            h = surf(X, Y, U);  
            xlim([0 0.05])        % Limit graph axes
            title(sprintf('Surface Plot at t = %.2f', t(i)));  
            xlabel('Thickness of Tile (m)');
            ylabel('Width of Tile (m)');
            zlabel('Temperature(C)');

            shading interp
            % set y axis range
            zlim([-0 1300])
            % set colour map range
            caxis ([0 1300]);
            colorbar
         end
    
        % Update the figure
        drawnow;
    
        % Pause for a short duration to control animation speed
        pause(0.001);  % Adjust as needed
    end