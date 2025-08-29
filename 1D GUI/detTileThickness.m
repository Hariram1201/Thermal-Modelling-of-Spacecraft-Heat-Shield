function [n] = detTileThickness(nt, nx,bestMethod, timeData, tempData)
    
%Function for determining the smallest tile thickness without any
%damage occurring to the tile

%Input Arguments:
%nt - number of timesteps
%nx - number of spatial steps
%bestMethod - stores the best method found from the analysis
%timeData - time points from the data provided (s)
%tempData - temperature points from the data provided

%Return Arguments:
%n - final tile thickness (m)

    %Initialises the Required Thickness as 0
    tileThickness = 0;

    %Traverses through the Tile Thicknesses to be Tested and Determines a Value
    % that is Appropriate
    for n = 0.12: -0.01:0.02
        % Plots a Graph of Temperature at the Inner Boundary against Time for 10
        % Different Tile Thicknesses
        [x , t , u ] = calctemp (4000 , nt , n , nx , bestMethod , timeData , tempData );
        figure (4)
        hold on
        plot (t , u (: , end ) )
        title ( ['Plot of Temperature at Inner Surface against Time with Suitable ' ...
            'TimeStep and Spatial Step '])
        xlabel ( ' Time ( s ) ')
        ylabel ( ' Temperature ( C ) on the Inner Surface ')
        legend ( ' 0.12 ' , ' 0.11 ' , ' 0.10 ' , ' 0.09 ' , ' 0.08 ' , ' 0.07 ' , ' 0.06 ' , ' 0.05 ' , ...
            ' 0.04 ' , ' 0.03 ' , ' 0.02 ')
        % Determines the First Tile Thickness at which the Maximum Temperature
        % Experienced by the Inner Surface is less than the Maximum Temperature
        % the Material of the Space Shuttle can Experience
        if max ( u (: , end ) ) > 150 && tileThickness == 0
            tileThickness = n + 0.01;
        end
    end
    n = tileThickness ;
end