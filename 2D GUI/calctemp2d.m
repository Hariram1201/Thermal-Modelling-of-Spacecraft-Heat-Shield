function [x, y, t, damageArr, u] = calctemp2d(tmax, nt, xmax, nx, method, timeData, tempData, tileWidth, ny)
% Function for modelling temperature in a space shuttle tile
% D N Johnston  14/02/24
%
% Input arguments:
% tmax   - maximum time (s)
% nt     - number of timesteps
% xmax   - total thickness (m)
% nx     - number of spatial steps
% method - solution method ('forward', 'backward' etc)
% timeData - time vector for surface temperatures (s)
% tempData - surface temperature vector (C or K)
% tileWidth - width of tile (m)
% ny - numer of spatial steps in the width of the tile
%
% Return arguments:
% x      - distance vector (m)
% t      - time vector (s)
% u      - temperature matrix (C or K)
%
% For example, to perform a  simulation with 501 time steps
%   [x, t, u] = calctemp(4000, 501, 0.05, 21, 'forward', timeData, tempData);
%

% Set material properties and derived values (LI-900)
% Obtained from NASA document: Structures and Materials: Space Shuttle Tiles, Grades 5-8 - NASA
% Note that we're assuming constant properties.
thermCon = 0.0484; % W/m K; 0.028 BTU/ft/hr/F at 70F and 1 atm
density  = 144;    % kg/m^3; 9 lb/ft^3
specHeat = 628;    % J/kg/K; 0.15 Btu/lb/F

% Initialise everything.
dt = tmax / (nt-1);
t = (0:nt-1) * dt;
dx = xmax / (nx-1);
x = (0:nx-1) * dx;
dy = tileWidth / (ny-1);
y = (0:ny-1) * dy;
u = zeros(ny, nx, nt);

%Determine value of p
thermDiff = thermCon / (density * specHeat);
p = (thermDiff * dt) / (dx ^ 2);

% Use interpolation to get outside temperature at time vector t 
% and store it as left-hand boundary vector L.
L = interp1(timeData, tempData, t, "linear", "extrap");

% % set initial conditions equal to boundary temperature at t=0.
u(:,:,1) = L(1);

temp2000 = 0;

% Select method and run simulation.
switch method
    %Simulates the Forward Differencing Method
    case 'forward'

        % Outside boundary condition
        for q =1:ny
            u(q,1,:) = L;
        end

        % set up index vectors
        i  = 2:nx-1;
        im = 1:nx-2;
        ip = 3:nx;
        j  = 2:ny-1;
        jm = 1:ny-2;
        jp = 3:ny;
        %now loop through time
        for n=1:nt-1   
            if (n * dt) < 2000 
                temp2000 = u(:,1,n+1);
            else               
                u(:,1,n+1) = temp2000;
            end

            % calculate internal values using forward differencing
            u(j, i, n+1) = (1 - 4 * p) * u(j, i, n) +  p * (u(j, im, n) + u(j, ip, n) + u(jm, i, n) + u(jp, i, n));

            %Neumann Boundary for north boundary (j=1)
            u(1,i,n+1) = (1-4*p)*u(1,i,n) + p*(u(1,im,n) + u(1,ip,n) + 2*u(2,i,n));
            %Neumann Boundary for south boundary (j=ny)
            u(ny,i,n+1) = (1 - 4*p) * u(ny,im,n) + p*(u(ny,im,n) + u(ny,ip,n) + 2*u(ny-1,i,n));
            %Neumann Boundary for west boundary (i=nx)
            u(j,nx,n+1) = (1-4*p) * u(j,nx,n) + p*(2*u(j,nx-1,n) + u(jm,nx,n) + u(jp,nx,n));
            %Neumann Boundary for north-west boundary (j=1 & i=nx)
            u(1,nx,n+1) = (1-4*p) * u(1,nx,n) + p*(2*u(2,nx,n) + 2*u(1,nx-1,n));
            %Neumann Boundary for south-west boundary (j=ny & i =nx)
            u(ny,nx,n+1) = (1-4*p) * u(ny,nx,n) + p*(2*u(ny-1,nx,n) + 2*u(ny,nx-1,n));  
        end

        %Analysis of Damage in the Tile
        [damageArr,u] = damageAnalysis(ny,nx,nt,u);

    %Simulates the Dufort-Frankel Method
    case 'dufort-frankel' 

        % Outside boundary condition
        for q =1:ny
            u(q,1,:) = L;
        end

        % set up index vectors
        i  = 2:nx-1;
        im = 1:nx-2;
        ip = 3:nx;
        j  = 2:ny-1;
        jm = 1:ny-2;
        jp = 3:ny;

        %now loop through time
        for n=1:nt-1   
            %Data only goes up to 2000s, so only Calculates Values for this Range only
            if (n * dt) < 2000 
                temp2000 = u(:,1,n+1);
            else               
                u(:,1,n+1) = temp2000;
            end

            %Slightly Different Equation is used to Prevent it from Failing
            if n ==1
                % calculate internal values using forward differencing
                u(j, i, n+1) = ((1 - 4 * p) * u(j, i, n) +  2*p * (u(j, im, n) + u(j, ip, n) + u(jm, i, n) + u(jp, i, n))) / (1+4*p);
                %Neumann Boundary for north boundary (j=1)
                u(1, i, n+1) = ((1 - 4 * p) * u(1, i, n) +  2*p * (u(1, im, n) + u(1, ip, n) + 2*u(2, i, n))) / (1+4*p);  
                %Neumann Boundary for south boundary (j=ny)
                u(ny, i, n+1) = ((1 - 4 * p) * u(ny, i, n) +  2*p * (u(ny, im, n) + u(ny, ip, n) + 2*u(ny-1, i, n))) / (1+4*p);
                %Neumann Boundary for west boundary (i=nx)
                u(j, nx, n+1) = ((1 - 4 * p) * u(j, nx, n) +  2*p * (2*u(j, nx-1, n) + u(jm, nx, n) + u(jp, nx, n))) / (1+4*p);
                %Neumann Boundary for north-west boundary (j=1 & i=nx)
                u(1, nx, n+1) = ((1 - 4 * p) * u(1, nx, n) +  2*p * (2*u(1, nx-1, n) + 2*u(2, nx, n))) / (1+4*p);
                %Neumann Boundary for south-west boundary (j=ny & i =nx)
                u(ny, nx, n+1) = ((1 - 4 * p) * u(ny, nx, n) +  2*p * (2*u(ny, nx-1, n) + 2*u(ny-1, nx, n))) / (1+4*p);
            else
                % calculate internal values using forward differencing
                u(j, i, n+1) = ((1 - 4 * p) * u(j, i, n-1) +  2*p * (u(j, im, n) + u(j, ip, n) + u(jm, i, n) + u(jp, i, n))) / (1+4*p);
                %Neumann Boundary for north boundary (j=1)
                u(1, i, n+1) = ((1 - 4 * p) * u(1, i, n-1) +  2*p * (u(1, im, n) + u(1, ip, n) + 2*u(2, i, n))) / (1+4*p);
                %Neumann Boundary for south boundary (j=ny)
                u(ny, i, n+1) = ((1 - 4 * p) * u(ny, i, n-1) +  2*p * (u(ny, im, n) + u(ny, ip, n) + 2*u(ny-1, i, n))) / (1+4*p);
                %Neumann Boundary for west boundary (i=nx)
                u(j, nx, n+1) = ((1 - 4 * p) * u(j, nx, n-1) +  2*p * (2*u(j, nx-1, n) + u(jm, nx, n) + u(jp, nx, n))) / (1+4*p);
                %Neumann Boundary for north-west boundary (j=1 & i=nx)
                u(1, nx, n+1) = ((1 - 4 * p) * u(1, nx, n-1) +  2*p * (2*u(1, nx-1, n) + 2*u(2, nx, n))) / (1+4*p);
                %Neumann Boundary for south-west boundary (j=ny & i =nx)
                u(ny, nx, n+1) = ((1 - 4 * p) * u(ny, nx, n-1) +  2*p * (2*u(ny, nx-1, n) + 2*u(ny-1, nx, n))) / (1+4*p);
            end
        end

        %Analysis of Damage in the Tile
        [damageArr,u] = damageAnalysis(ny,nx,nt,u);

    %Simulates the Backwards Differencing Method
    case 'backward' 

        % Outside boundary condition
        for q =1:ny
            u(q,1,:) = L;
        end

        maxiterations = 100;
        tolerance = 1.e-4;

        %now loop through time
        for n=1:nt-1
            %Data only goes up to 2000s, so only Calculates Values for this Range only
            if (n * dt) < 2000 
                temp2000 = u(:,1,n+1);
            else               
                u(:,1,n+1) = temp2000;
            end
    
            % calculate internal values iteratively using Gauss-Seidel
            % Starting values are equal to old values
            u(2:ny-1, 2:nx-1, n+1) = u(2:ny-1, 2:nx-1, n);
    
            for iteration = 1:maxiterations
                change = 0;
                for i=2:nx
                   
                    for j=1:ny
                        %Neumann Boundary for north boundary (j=1)
                        if j == 1
                            jm = 2;
                        else
                            jm = j -1;
                        end

                        %Neumann Boundary for south boundary (j=ny)
                        if j == ny
                            jp = ny - 1;
                        else
                            jp = j + 1;
                        end

                        %Neumann Boundary for west boundary (i=nx)
                        if i == nx
                            ip = nx - 1;
                        else
                            ip = i + 1;
                        end

                        uold = u(j, i, n+1);
                        % calculate internal values using forward differencing
                        u(j, i, n+1) = ((u(j, i, n) + p * (u(jm, i, n+1) + u(jp, i, n+1) + u(j, i-1, n+1) + u(j, ip, n+1)))/(1+4*p));
                        change = change + abs(u(j, i, n+1) - uold);
                    end
                end

                if change < tolerance
                    break
                end
            end
        end

        %Analysis of Damage in the Tile
        [damageArr,u] = damageAnalysis(ny,nx,nt,u);

    %Simulates the Crank-Nicholson Method
    case 'crank-nicholson'

        % Outside boundary condition
        for q =1:ny
            u(q,1,:) = L;
        end

        maxiterations = 100;
        tolerance = 1.e-4;

        %now loop through time
        for n=1:nt-1
            %Data only goes up to 2000s, so only Calculates Values for this Range only
            if (n * dt) < 2000 
                temp2000 = u(:,1,n+1);
            else               
                u(:,1,n+1) = temp2000;
            end
    
            % calculate internal values iteratively using Gauss-Seidel
            % Starting values are equal to old values
            u(2:ny-1, 2:nx-1, n+1) = u(2:ny-1, 2:nx-1, n);

            for iteration = 1:maxiterations
                change = 0;
                for i=2:nx
                   
                    for j=1:ny
                        %Neumann Boundary for north boundary (j=1)
                        if j == 1
                            jm = 2;
                        else
                            jm = j -1;
                        end

                        %Neumann Boundary for south boundary (j=ny)
                        if j == ny
                            jp = ny - 1;
                        else
                            jp = j + 1;
                        end

                        %Neumann Boundary for west boundary (i=nx)
                        if i == nx
                            ip = nx - 1;
                        else
                            ip = i + 1;
                        end

                        uold = u(j, i, n+1);
                        % calculate internal values using forward differencing
                        u(j, i, n+1) = ((1-2*p) * u(j,i,n) + (p/2) * (u(j,i-1,n) + u(j,ip,n) + u(jm,i,n) + u(jp,i,n) + u(j,i-1,n+1) + u(j,ip,n+1) + u(jm,i,n+1) + u(jp,i,n+1))) / (1+2*p);
                        change = change + abs(u(j, i, n+1) - uold);
                    end
                end

                if change < tolerance
                    break
                end
            end
        end

        %Analysis of Damage in the Tile
        [damageArr,u] = damageAnalysis(ny,nx,nt,u);


    %Returns an Error Message if an Appropriate Method is not Entered
    otherwise
        error (['Undefined method: ' method]);
end

    