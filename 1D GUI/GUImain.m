function GUI()

    % Create a figure
    fig = figure('Position', [500, 400, 800, 600], 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', 'Name', 'Temperature variation through a space shuttle tile 1D');

    % Define the list of options and corresponding images
    options = {'Tile 468', 'Tile 480', 'Tile 502', 'Tile 590', 'Tile 597', 'Tile 711', 'Tile 730', 'Tile 850'};
    images = {'temp468R.jpg', 'temp480R.jpg', 'temp502R.jpg', 'temp590R.jpg', 'temp597R.jpg', 'temp711R.jpg','temp730R.jpg', 'temp850R.jpg'}; % Paths to image files

    % Create listbox for options
    listbox = uicontrol('Style', 'listbox', 'Position', [350, 445, 100, 110], 'String', options, 'Callback', @listboxCallback);

    % Create axes for displaying images
    axesHandle = axes('Position', [0.02, 0.50, 0.4, 0.6], 'Title', 'Boundary Conditions');
    axesHandle2 = axes('Position', [0.09, 0.055, 0.85, 0.6]);

    % Insert a plot button
    button = uicontrol('Style', 'pushbutton', 'String', 'Simulate', 'Position', [350, 410, 100, 30], 'Callback', @buttonCallback)

    % Load 468 tile as inital tile
    img = imread(images{1});
    [x,t,u,nt, nx, n, bestMethod] = simulation('temp468R.jpg');
    imshow(img, 'Parent', axesHandle);
    surf(x,t,u,'Parent', axesHandle2)
    xlabel(axesHandle2, 'Distance(m)')
    ylabel(axesHandle2, 'Time(s)')
    zlabel(axesHandle2, 'Temperature(c)')
  

    textLabelBC = uicontrol(fig, 'Style', 'text', 'Position', [32, 550, 300, 30], 'String', 'External boundary conditions', 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabelFG = uicontrol(fig, 'Style', 'text', 'Position', [64, 350, 100, 30], 'String', 'Simulation:', 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabelnx = uicontrol(fig, 'Style', 'text', 'Position', [470, 525, 150, 30], 'String', 'no. Spatial steps:', 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabelnt = uicontrol(fig, 'Style', 'text', 'Position', [476, 495, 150, 30], 'String', 'no. Time steps:', 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabeln = uicontrol(fig, 'Style', 'text', 'Position', [480, 465, 150, 30], 'String', 'Tile thickness:' , 'HorizontalAlignment', 'center', 'FontSize', 12)
    textLabelBM = uicontrol(fig, 'Style', 'text', 'Position', [483, 435, 150, 30], 'String', 'Best method:', 'HorizontalAlignment', 'center', 'FontSize', 12)

    outputFieldnx = uicontrol(fig, 'Style', 'edit', 'Position', [610, 535, 80, 20], 'HorizontalAlignment', 'left', 'Enable', 'inactive');
    outputFieldnt = uicontrol(fig, 'Style', 'edit', 'Position', [610, 505, 80, 20], 'HorizontalAlignment', 'left', 'Enable', 'inactive');
    outputFieldn = uicontrol(fig, 'Style', 'edit', 'Position', [610, 475, 80, 20], 'HorizontalAlignment', 'left', 'Enable', 'inactive');
    outputFieldBM = uicontrol(fig, 'Style', 'edit', 'Position', [610, 445, 80, 20], 'HorizontalAlignment', 'left', 'Enable', 'inactive');

    outputFieldnx.String = sprintf(num2str(nx))
    outputFieldnt.String = sprintf(num2str(nt))
    outputFieldn.String = sprintf(num2str(n))
    outputFieldBM.String = sprintf(num2str(bestMethod))

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

    % Call back function for button
    function buttonCallback(~, ~)
        selectedIndex = listbox.Value;
        selectedOption = images{selectedIndex};
        [x, t, u, nt, nx, n, bestMethod] = simulation(selectedOption);
        surf(x,t,u,'Parent', axesHandle2)
        xlabel(axesHandle2, 'Distance(m)')
        ylabel(axesHandle2, 'Time(s)')
        zlabel(axesHandle2, 'Temperature(c)')

        % Output variables ot main window
        outputFieldnx.String = sprintf(num2str(nx))
        outputFieldnt.String = sprintf(num2str(nt))
        outputFieldn.String = sprintf(num2str(n))
        outputFieldBM.String = sprintf(num2str(bestMethod))
    end

    

end
