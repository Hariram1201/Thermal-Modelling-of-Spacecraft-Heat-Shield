function[x, t, u, nt, nx, n, bestMethod] = simulation(selectedOption)

% Function - simulation: Performs numerical analysis to simulate heat flow
% through a space shuttle tile. Determines appropriate spatial and time
% steps, tile thickness, and method
%
% Input arguments:
% selectedOption - Selected image containing tile inital locations
%
% Ouput arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
% nt     - number of timesteps
% nx     - number of spatial steps
% n      - Tile thickness (m)
% bestMethod - Best numerical method for simulation




% Retrieves the time and temperature data points  
[origin, TRC, time, temperature] = ImgScan(selectedOption);
[tempData, timeData] = DataScale(origin, TRC, time, temperature);

%Determines the temperature at each distance for each time step for the 4
%different methods
[x1, t1, u1] = calctemp(4000, 1001, 0.05, 21, 'forward', timeData, tempData);
[x2, t2, u2] = calctemp(4000, 501, 0.05, 21, 'dufort-frankel', timeData, tempData);
[x3, t3, u3] = calctemp(4000, 501, 0.05, 21, 'backward', timeData, tempData);
[x4, t4, u4] = calctemp(4000, 501, 0.05, 21, 'crank-nicholson', timeData, tempData);

%Investigates the Stability and Accuracy of the Four Methods
[finalTimeStep] = investMethodsTime(tempData, timeData);
[finalSpaceStep] = investMethodsSpace(tempData, timeData);

%Determines which of the Four Methods is the Most Accurate
[bestMethod, index] = detBestMethod(finalTimeStep, finalSpaceStep);

%Determines the Values of dt and nt that are within the Required Accuracy
%whilst also Saving Computational Power
nx = finalSpaceStep(index);
dt = finalTimeStep(index);
nt = round((4000 / dt) + 1);

%Determines the Final Thickness of the Tile that will not affect the Tile
n = detTileThickness(nt, nx,bestMethod, timeData, tempData);

%Plots the Final Graph
[x, t, u] = calctemp(4000, nt, n, nx, bestMethod, timeData, tempData);
figure(5)
surf(x,t,u)
title('Final Graph, Interactive')
xlabel('Distance(m)')
ylabel('Time(s)')
zlabel('Temperature(C)')
end