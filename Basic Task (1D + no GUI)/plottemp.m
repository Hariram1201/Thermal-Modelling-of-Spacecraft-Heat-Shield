% Script to plot image of measured temperature, and trace it using the mouse.
%
% Image from http://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-aiaa-2001-0352.pdf
%
% D N Johnston 05/01/24

%Retrieves the image with graph of temperature against time 
name = 'temp597';
img=imread([name '.jpg']);

image(img);
title("Click left button to set data points, right button to end");

%Creates the arrays needed to store the data points to be used
timeData = [];
tempData = [];

hold on
% Gets origin location and places a blue cross at this point
[xOrigin, yOrigin, button] = ginput(1);
plot(xOrigin, yOrigin, '+b')

%Gets the top right location and places a red cross at this point
[xEdge, yEdge, button] = ginput(1);
plot(xEdge, yEdge, '+r')

%Determines the multiplier to determine the actual coordinates of the
%points to be used
xDiff = xEdge - xOrigin;
yDiff = yEdge - yOrigin;

xMultiplier = 2000 / xDiff;
yMultiplier = 2000 / yDiff;

while 1 % infinite loop
    [x, y, button] = ginput(1); % get one point using mouse
    if button ~= 1 % break if anything except the left mouse button is pressed
        break
    end

    %Determines the actual coordinates of the points to be used
    xPoint = ((x - xOrigin) * xMultiplier);
    yPoint = ((y - yOrigin) * yMultiplier);
    
    % Add data point to vectors. Note that x and y are pixel coordinates.
    % You need to locate the pixel coordinates of the axes, interactively
    % or otherwise, and scale the values into time (s) and temperature (F, C or K).
    %Only adds the data points to the array if they lie within the
    %constraints of the graph
    if (xPoint > 0 && xPoint < 2000) && (yPoint > 0 && yPoint < 2000)

        %Converts  the Temperature from Fahrenheits to Celsius
        yPoint = (yPoint - 32 ) / 1.8;

        %Plots a green circle on the graph for the selected data points
        plot(x, y, 'og')

        %Adds the new points into the array
        timeData = [timeData, xPoint];
        tempData = [tempData, yPoint];
    end
end
hold off

%save data to .mat file with same name as image file
save(name, 'timeData', 'tempData')

