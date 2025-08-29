function [ finalTimeStep ] = investMethodsTime (tempData,timeData)
%Function to Determine the Stability and Accuracy of the Four Methods
%through Determining the Inner Surface Temperature at Time = 4000s 
%against the Timestep

%Return Arguments:
%finalTimeStep - Array containing the Largest Possible Timestep for each of
%the Four Methods before the Desired Accuracy is Compromised


    %Retrieves the Temperature and Time Data to be used and Stored in the
    %Arrays 'tempData' and 'timeData'
    load temp597.mat

    %Initialises the Variables 
    nx = 21;
    thick = 0.05;
    tmax = 4000;

    % numbers of timesteps
    nt = 41:20:1001;

    % caculate all timestep sizes
    dt = tmax./(nt-1);

    % preallocate result vectors for efficiency
    uf = zeros(size(nt));
    ub = zeros(size(nt));
    ud = zeros(size(nt));
    uc = zeros(size(nt));

    % Run simulations for each nt value for all four methods
    for i = 1:length(nt)
        disp (['nt = ' num2str(nt(i)) ', dt = ' num2str(dt(i)) ' s'])

        %Determines the Temperature at the Inner Surface Boundary for the
        %Forward Differencing Method
        [~, ~, u] = calctemp(tmax, nt(i), thick, nx, 'forward', timeData, tempData);
        uf(i) = u(end, nx);

        %Determines the Temperature at the Inner Surface Boundary for the
        %Backward Differencing Method
        [~, ~, u] = calctemp(tmax, nt(i), thick, nx, 'backward', timeData, tempData);
        ub(i) = u(end, nx);

        %Determines the Temperature at the Inner Surface Boundary for the
        %Dufort-Frankel Method
        [~, ~, u] = calctemp(tmax, nt(i), thick, nx, 'dufort-frankel', timeData, tempData);
        ud(i) = u(end, nx);

        %Determines the Temperature at the Inner Surface Boundary for the
        %Crank Nicholson Method
        [~, ~, u] = calctemp(tmax, nt(i), thick, nx, 'crank-nicholson', timeData, tempData);
        uc(i) = u(end, nx);
    end

    %Plots a figure to determine accuracy and stability of each method
    figure(2)
    plot(dt, [uf; ub; ud; uc])
    ylim([0 200])
    title('Investigation in Stability and Accuracy of Different Timesteps')
    xlabel('Time Step')
    ylabel('Temperature(C) on the Inner Surface at t = 4000s')
    legend ('Forward', 'Backward', 'Dufort-Frankel', 'Crank-Nicholson')

    %Compiles the Final Times for the Four Methods and Assigns to a New
    %Variable called 'time' and also extracts the time at the points where
    %the Time Step is at its' Lowest Test Value
    time = [uf;ub;ud;uc];
    startTime = [uf(end), ub(end), ud(end), uc(end)];

    %Sets Accuracy Desired for Retrieving a Suitable Time Step
    accuracy = 0.05;

    %Determines the Range in which the Temperature Values are Allowed to
    %Lie within
    for i = 1:4
        plus5 = startTime(i) + accuracy * startTime(i);
        minus5 = startTime(i) - accuracy * startTime(i);

        plusMinus5Time(2*i - 1) = plus5;
        plusMinus5Time(2*i) = minus5;
    end

    %Initialises the Array to Store the Largest Possible Step for each
    %Method before the Accuracy of the Solution is Compromised
    finalTimeStep = [0,0,0,0];

    %Determines the Time Step allowed for Each Method before the Accuracy
    %of the Solution is Compromised
    for n = 1:4
        for x = 1:length(time(1,:))
            testCol = time(:,x);
            if testCol(n) > plusMinus5Time(2*n - 1) || testCol(n) < plusMinus5Time(2*n)
                finalTimeStep(n) = x;
            end
        end
        finalTimeStep(n) = 41 + 20*finalTimeStep(n);
        finalTimeStep(n) = tmax./(finalTimeStep(n)-1);
    end
