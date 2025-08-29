function main()

% Function to load a GUI enabling user to select and simulate the heat flow and damage through tiles at different locations as they burn up. 

    % Create a figure
    fig = figure('Position', [500, 400, 800, 600], 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', 'Name', 'Temperature variation through a space shuttle tile 2D');

    % List of initial temperature conditions and images for each tile
    % location
    options = {'Tile 468', 'Tile 480', 'Tile 502', 'Tile 590', 'Tile 597', 'Tile 711', 'Tile 730', 'Tile 850'};
    images = {'temp468R.jpg', 'temp480R.jpg', 'temp502R.jpg', 'temp590R.jpg', 'temp597R.jpg', 'temp711R.jpg','temp730R.jpg', 'temp850R.jpg'}; % Paths to image files

    % Create listbox for tile locations
    listbox = uicontrol('Style', 'listbox', 'Position', [350, 445, 100, 110], 'String', options, 'Callback', @listboxCallback);

    % Create axes for displaying images
    axesHandle = axes('Position', [0.02, 0.50, 0.4, 0.6], 'Title', 'Boundary Conditions');
    axesHandle2 = axes('Position', [0.09, 0.075, 0.85, 0.55]);

    % List of methods 
    optionsM = {'forward', 'dufort-frankel', 'backward', 'crank-nicholson' };
    
    % Create listbox for methods
    listboxM = uicontrol('Style', 'listbox', 'Position', [470, 445, 100, 110], 'String', optionsM, 'Callback', @listboxCallbackM);

    % Simulation button
    button = uicontrol('Style', 'pushbutton', 'String', 'Simulate', 'Position', [600, 480, 100, 30], 'Callback', @buttonCallback)

    % Load 468 tile as inital tile & display parameters
    img = imread(images{1});
    [x,y,t,u] = simulation('temp468R.jpg', 'forward');
    imshow(img, 'Parent', axesHandle);
  
    % Create lables for parameters
    textLabelBC = uicontrol(fig, 'Style', 'text', 'Position', [32, 550, 300, 30], 'String', 'External boundary conditions', 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabelFG = uicontrol(fig, 'Style', 'text', 'Position', [64, 350, 100, 30], 'String', 'Simulation:', 'HorizontalAlignment', 'center', 'FontSize', 12)
    

    % Callback function for listbox
    function listboxCallback(~, ~)
        % Get selected index
        selectedIndex = listbox.Value;
        % Get selected option
        selectedOption = options{selectedIndex};
        % Display selected option
        disp(['Selected Option: ', selectedOption]);
        % Load and display corresponding image
        img = imread(images{selectedIndex});
        imshow(img, 'Parent', axesHandle);
    end

% Callback function for listboxM
    function listboxCallbackM(~, ~)
        % Get selected index
        selectedIndex = listboxM.Value;
        % Get selected option
        method = optionsM{selectedIndex};
        % Display selected option
        disp(['Selected Option: ', method]);
    end


    % Call back function for button
    function buttonCallback(~, ~)
        selectedIndex = listboxM.Value;
        method = optionsM{selectedIndex};
        selectedIndex = listbox.Value;
        selectedOption = images{selectedIndex};
        [x, y, t, u] = simulation(selectedOption, method);
        [X, Y] = meshgrid(x, y);

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
            h = surf(X, Y, U, 'Parent', axesHandle2);  
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
        
    end

    

end
