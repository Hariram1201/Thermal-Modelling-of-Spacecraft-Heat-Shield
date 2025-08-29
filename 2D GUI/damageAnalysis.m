function [damageArr,u] = damageAnalysis(ny,nx,nt,u)

%Function to model the damage in the tile 
%Input Arguments:
% ny - number of spatial steps in the thickness of the tile
% nx - number of spatial steps in the width of the tile
% nt - number of time steps
% u - temperatre at each subsection of the tile

%Return arguements:
% damageArr - Details the areas in which the tile is damaged, used to
% animate the graph
% u - stores the temperature values at each section of the tile

    damageArr = zeros(ny,nx,nt);

    %Assigns the temperates at time = 0 into the damage array
    damageArr(:,:,1) = u(:,:,1);
    %Determining whether the section of the tile is damaged by checking
    %whether it exceeds the maximum allowable temperature or if the tile
    %subsection has already been damaged
    for n = 2:nt
        for a = 1:nx
            for b = 1:ny
                damageArr(b,a,n) = u(b,a,n);
                if u(b,a,n) > 660 || isnan(damageArr(b,a,n-1))
                    u(b,a,n) = u(b,1,n);
                    damageArr(b,a,n) = NaN;
                end
            end
        end
     end