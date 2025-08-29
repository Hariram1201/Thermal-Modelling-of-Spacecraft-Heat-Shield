function [tempData,timeData] = DataScale(origin, TRC, time, temperature)
% DataScale - Scales the data extracted from ImgScan using the origin and
% TRC as scaling factors 
%
% Input arguments:
% origin           - 1 x 2 matrix, graph origin
% TRC              - 1 x 2 matrix, graph Top Right Corner
% time             - Time matrix containing all read time values
% temperature      - Temperature matrix containing all read time values
%
% Output arguments
% tempData         - Scaled temperature data 
% timeData         - Scaled time data

% Determine multiplier to scale pixel co-ordinates

xDiff = TRC(1) - origin(1);
yDiff = TRC(2) - origin(2);

xMultiplier = 2000 / xDiff;
yMultiplier = 2000 / yDiff;

% Determine actual locations of the points to be used

for n = 1:length(time)

    % Scaling time data
    timeData(n) = ((time(n) - origin(1)) * xMultiplier);

    % Scaling and converting temperature data F -> C
    tempData(n) = (((temperature(n) - origin(2)) * yMultiplier) - 32) /1.8;

end