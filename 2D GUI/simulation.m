function[x,y,t,u] = simulation(selectedOption, method)

% Function - simulation: Simulates the damage across a tile using different
% numerical methods
%
% Input arguments:
% method - The method to be used in the numerical analysis
%
% Output arguments
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
% y      - width of tile(m)



% Retrieves the time and temperature data points
[origin, TRC, time, temperature] = ImgScan(selectedOption);
[tempData, timeData] = DataScale(origin, TRC, time, temperature);

%Initialises Tile Width and Number of Spatial Steps across the Width
tileWidth = 0.2;
ny = 20;

%Simultes Forward Method
[x, y, t, u] = calctemp2d(4000, 2001, 0.05, 21, method, timeData, tempData, tileWidth, ny);


end