function [finalSpaceStep] = investMethodsSpace(tempData, timeData)
%Function to Determine the Stability and Accuracy of the Four Methods
%through Determining the Inner Surface Temperature at Time = 4000s 
%against the Number of Spatial Steps

%Input Arguments:
%tempData: Array containing initial temperature data
%timeData: Array containing initial time data

%Return Arguments:
%finalSpaceStep - Array containing the Smallest Number of Space Step for 
% each of the Four Methods before the Desired Accuracy is Compromised


    %Initialises the Variables 
    nt = 1001;
    thick = 0.05;
    tmax = 4000;

    % numbers of spatial steps
    nx = 2:1:21;

    % caculate all spatial step sizes
    dx = thick./(nx-1);

    % preallocate result vectors for efficiency
    uf = zeros(size(nx));
    ub = zeros(size(nx));
    ud = zeros(size(nx));
    uc = zeros(size(nx));

    % Run simulations for each nx value for all four methods
    for i = 1:length(nx)

        disp (['nx = ' num2str(nx(i)) ', dx = ' num2str(dx(i)) ' m'])

        %Determines the Temperature at the Inner Surface Boundary for the
        %Forward Differencing Method
        [~, ~, u] = calctemp(tmax, nt, thick, nx(i), 'forward', timeData, tempData);
        uf(i) = u(nt,end);

        %Determines the Temperature at the Inner Surface Boundary for the
        %Backward Differencing Method
        [~, ~, u] = calctemp(tmax, nt, thick, nx(i), 'backward', timeData, tempData);
        ub(i) = u(nt, end);

        %Determines the Temperature at the Inner Surface Boundary for the
        %Dufort-Frankel Method
        [~, ~, u] = calctemp(tmax, nt, thick, nx(i), 'dufort-frankel', timeData, tempData);
        ud(i) = u(nt, end);
        
        %Determines the Temperature at the Inner Surface Boundary for the
        %Crank Nicholson Method    
        [~, ~, u] = calctemp(tmax, nt, thick, nx(i), 'crank-nicholson', timeData, tempData);
        uc(i) = u(nt, end);
    end


    %Compiles the Final Distances for the Four Methods and Assigns to a New
    %Variable called 'space' and also extracts the temperature at the points 
    %where the Number of Space Step is at its' Highest Test Value
    space = [uf;ub;ud;uc];
    startSpace = [uf(end), ub(end), ud(end), uc(end)];

    %Sets Accuracy Desired for Retrieving a Suitable Number of Space Steps
    accuracy = 0.001;

    %Determines the Range in which the Temperature Values are Allowed to
    %Lie within
    for i = 1:4
        plus5 = startSpace(i) + accuracy * startSpace(i);
        minus5 = startSpace(i) - accuracy * startSpace(i);

        plusMinus5Space(2*i - 1) = plus5;
        plusMinus5Space(2*i) = minus5;
    end

    %Initialises the Array to Store the Smallest Possible Number of Space 
    %Step for each Method before the Accuracy of the Solution is Compromised
    finalSpaceStep = [21,21,21,21];

    %Determines the Smallest Number of Space Steps allowed for Each Method 
    %before the Accuracy of the Solution is Compromised
    for n = 1:4
        for x = 1:length(space(1,:))
            testCol = space(:,x);
            if testCol(n) > plusMinus5Space(2*n - 1) || testCol(n) < plusMinus5Space(2*n)
                finalSpaceStep(n) = x+1;
            end
        end
    end
end