function [origin, TRC, time, temperature] = ImgScan(image_name)

% ImgScan - Scans an image of the inital boundary conidtions for a space
% shuttle tile recording the pixel locations of all the data in the image
% along with the origin and top right corner for scaling.
%
% Input arguments:
% image_name       - Name of image for inital conditions i.e. 'temp597.jpg'
%
% Ouput arguments:
% origin           - 1 x 2 matrix, graph origin
% TRC              - 1 x 2 matrix, graph Top Right Corner
% time             - Time matrix containing all read time values
% temperature      - Temperature matrix containing all read time values


    % Import image for processing and turn it into a grayscale matrix
    image = imread(image_name);
    GreyImage = image(:,:,3);

    % Preallocate matrix size for speed
    [x, y] = size(GreyImage);
    p = zeros(x, y); % Used to find axes
    q = zeros(x, y); % Used to find data

    for n = 1:x
        for m = 1:y
            % Use colour filters to identify white and non white points
            if  GreyImage(n,m) < 100 
                p(n,m) = 0; % Marks non-white points with a 0
            else
                p(n,m) = 1; % Marks white points with a 1
            end

            % Use colour filters to identify red points
            if GreyImage(n,m) < 100 &&  20 < GreyImage(n,m)
               q(n,m) = 0;  % Marks red points with a 0
            else
               q(n,m) = 1;  % Marks non-red points with a 1
            end
        end
    end


    % Finding Origin
    % Sum rows and columns
    rowsum = sum(p, 2);
    colsum = sum(p, 1);

    % Lowest sums will be the Axes
    [row_max_0, loc_row_max_0] = min(rowsum);  
    [col_max_0, loc_col_max_0] = min(colsum);

    % Use this to determine the origin
    origin = [loc_col_max_0, loc_row_max_0]; 

    % Finding Top Right Corner (TRC) for scaling

    % column and row data for analysis
    columnData = p(:, origin(1));
    rowData = p(origin(2),:);

    % Initialise variables
    maxZeroLengthR = 0;
    currentZeroLength = 0;

    % Analysing rows
    for i = 1:height(columnData)
        % If data is zero, increment currentZeroLength count
        if columnData(i) == 0
            currentZeroLength = currentZeroLength + 1;
        else
            % Non-zero element, update maxZeroLength if required
            if currentZeroLength > maxZeroLengthR
                maxZeroLengthR = currentZeroLength;
            end
            % Reset currentZeroLength 
            currentZeroLength = 0;
        end
    end

    % Update maxZeroLength if the longest sequence is at the end
    if currentZeroLength > maxZeroLengthR
        maxZeroLengthR = currentZeroLength;
    end

    % Reset variables
    maxZeroLengthC = 0; 
    currentZeroLength = 0;

    % Analysing columns
    for j = 1:size(rowData, 2)
        % If data is zero, increment currentZeroLength count
        if rowData(j) == 0
            currentZeroLength = currentZeroLength + 1;
        else
            % Non-zero element, update maxZeroLength if required
            if currentZeroLength > maxZeroLengthC
                maxZeroLengthC = currentZeroLength;
            end
            % Reset currentZeroLength
            currentZeroLength = 0;
        end
    end

    % Update maxZeroLengthC if the longest sequence is at the end
    if currentZeroLength > maxZeroLengthC
        maxZeroLengthC = currentZeroLength;
    end

    % Using origin and length of axes to find top right corner (TRC) of graph axes
    TRC = [origin(1) + maxZeroLengthC, origin(2) - maxZeroLengthR];


    % Reformatting q using p to only include data within the axes

    % Adjust columns
    q(:, 1:origin(1)+1) = 1;
    % Adjust rows
    q(origin(2):x, :) = 1;
    q(1:TRC(2), :) = 1;

    % Removing double data locations, causes errors
    % Largest data point, worst case (highest temperature)
    for m = 1:length(q)
        for n = 1:height(q)
            if q(n,m) == 0
                MaxZeroRow = n;      % Record which row the first zero occurs in
                q(:,m) = 1;          % Set all other values in the column to 1
                q(MaxZeroRow,m) = 0; % Reset location of highest zero
            end
        end
    end

    % Storing row and column data of zeros
    % Preallocate time and temperature arrays
    time = [];
    temperature = [];

    % Loop through q array
    for i = 1:height(q)
        for j = 1:length(q)
            % Identify 0 locations and record rows and columns as time and temperature
            if q(i,j) == 0
                time = [time, j];
                temperature = [temperature, i];
            end
        end
    end
end
