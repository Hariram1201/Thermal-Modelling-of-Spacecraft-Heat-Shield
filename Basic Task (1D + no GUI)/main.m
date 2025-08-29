% Retrieves the time and temperature data points
load('temp597.mat');
  
% timeData = [0, 60, 500, 1000, 1500, 1750, 4000]; % s
% tempData = [16, 16, 820, 760, 440, 16, 16]; % degrees C

%Determines the temperature at each distance for each time step for the 4
%different methods
[x1, t1, u1] = calctemp(4000, 1001, 0.05, 21, 'forward', timeData, tempData);
[x2, t2, u2] = calctemp(4000, 501, 0.05, 21, 'dufort-frankel', timeData, tempData);
[x3, t3, u3] = calctemp(4000, 501, 0.05, 21, 'backward', timeData, tempData);
[x4, t4, u4] = calctemp(4000, 501, 0.05, 21, 'crank-nicholson', timeData, tempData);

%Plots a subplot for forward differencing
subplot(2,2,1)
surf(x1, t1, u1)
title('Forward Differencing')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')
shading interp

%Plots a subplot for Dufort-Frankel
subplot(2,2,2)
surf(x2, t2, u2)
title('Dufort-Frankel')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')
shading interp

%Plots a subplot for backward differencing
subplot(2,2,3)
surf(x3, t3, u3)
title('Backward Differencing')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')
shading interp

%Plots a subplot for Crank Nicholson
subplot(2,2,4)
surf(x4, t4, u4)
title('Crank Nicholson')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')
shading interp

%Investigates the Stability and Accuracy of the Four Methods
[finalTimeStep] = investMethodsTime (tempData,timeData);
[finalSpaceStep] = investMethodsSpace (tempData,timeData);

%Determines which of the Four Methods is the Most Accurate
[bestMethod, index] = detBestMethod(finalTimeStep, finalSpaceStep);

%Determines the Values of dt and nt that are within the Required Accuracy
%whilst also Saving Computational Power
nx = round(0.05 / finalSpaceStep(index) + 1);
dt = finalTimeStep(index);
nt = round((4000 / dt) + 1);

%Determines the Final Thickness of the Tile that will not affect the Tile
n = detTileThickness(nt, nx,bestMethod, timeData, tempData);

%Plots the Final Graph
[x, t, u] = calctemp(4000, nt, n, nx, bestMethod, timeData, tempData);
figure(5)
surf(x,t,u)
title('Final Graph')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')