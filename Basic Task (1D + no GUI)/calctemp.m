function [x, t, u] = calctemp(tmax, nt, xmax, nx, method, timeData, tempData)
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
u = zeros(nt, nx);

iVec = 2:nx-1;
a = 1:nx;

%Determine value of p
thermDiff = thermCon / (density * specHeat);
p = (thermDiff * dt) / (dx ^ 2);

% Use interpolation to get outside temperature at time vector t 
% and store it as left-hand boundary vector L.
L = interp1(timeData, tempData, t, "linear", "extrap");

% set initial conditions equal to boundary temperature at t=0.
u(1, :) = L(1);

temp2000 = 0;

% Select method and run simulation.
switch method
    %Simulates the Forward Differencing Method
    case 'forward'

        % Outside boundary condition
        u(:, 1) = L;

        %Traverses through each Time Step (Rows)
        for n=1:nt-1
            %Data only goes up to 2000s, so only Calculates Values for this
            %Range only
            if (n * dt) < 2000 

                temp2000 = u(n+1,1);

            else
                
                u(n+1,1) = temp2000;
            end

            %Determines the Value for the Next Temperature
            u(n+1,iVec) = (1 - 2 * p) * u(n,iVec) + p * (u(n,iVec-1) + u(n,iVec+1));

            %Determines the Values for the Internal Boundary
            u(n+1,nx) = (1-2*p) * u(n, nx) + 2*p*u(n, nx-1);
        end

    %Simulates the Dufort-Frankel Method
    case 'dufort-frankel' 

        % Outside boundary condition
        u(:, 1) = L; 

        %Traverses through each Time Step (Rows)
        for n = 1:nt-1
            %Data only goes up to 2000s, so only Calculates Values for this
            %Range only
            if (n * dt) < 2000 

                temp2000 = u(n+1,1);
            else

                u(n+1,1) = temp2000;
            end

            %Determines whether the Method is on the Second Row where a
            %Slightly Different Equation is used to Prevent it from
            %Failing
            if (n + 1) == 2
                %Determines the Value for the Next Temperature
                u(n+1,iVec) = ((1-2*p) * u(n,iVec) + 2*p*(u(n,iVec-1) + u(n,iVec+1))) / (1+2*p);

                %Determines the Values for the Internal Boundary
                u(n+1,nx) = ((1 - 2*p) * u(n,nx) + 4*p*u(n,nx-1)) / (1+2*p); 
            else  
                %Determines the Value for the Next Temperature
                u(n+1,iVec) = ((1-2*p) * u(n-1,iVec) + 2*p*(u(n,iVec-1) + u(n,iVec+1))) / (1+2*p);
                    
                %Determines the Values for the Internal Boundary
                u(n+1,nx) = ((1 - 2*p) * u(n-1,nx) + 4*p*u(n,nx-1)) / (1+2*p);
            end
        end

    %Simulates the Backwards Differencing Method
    case 'backward' 

        %Traverses through each Time Step (Rows)
        for n=1:nt-1
            if (n*dt) < 2000

                %Sets Values to be used for the Tri-Diagonal Matrix Method
                b(1)    = 1;
                c(1)    = 0;
                d(1)    = L(n+1);
                a(iVec) = -p;
                b(iVec) = 1 + 2*p;
                c(iVec) = -p;
                d(iVec) = u(n,iVec);
                a(nx)   = -2*p;
                b(nx)   = 1 + 2*p;
                d(nx)   = u(n,nx);

                %Updates the Next Column with the Calculated Temperature
                %Values
                u(n+1,:) = tdm(a,b,c,d);

                temp2000 = L(n+1);

            else
                %Sets Values to be used for the Tri-Diagonal Matrix Method
                b(1)    = 1;
                c(1)    = 0;
                d(1)    = temp2000;
                a(iVec) = -p;
                b(iVec) = 1 + 2*p;
                c(iVec) = -p;
                d(iVec) = u(n,iVec);
                a(nx)   = -2*p;
                b(nx)   = 1 + 2*p;
                d(nx)   = u(n,nx);

                %Updates the Next Column with the Calculated Temperature
                %Values
                u(n+1,:) = tdm(a,b,c,d);
            end

        end

    %Simulates the Crank-Nicholson Method
    case 'crank-nicholson'

        %Traverses through each Time Step (Rows)
        for n=1:nt-1
 
            if (n*dt) < 2000
                %Sets Values to be used for the Tri-Diagonal Matrix Method
                b(1)    = 1;
                c(1)    = 0;
                d(1)    = L(n+1);
                a(iVec) = -p/2;
                b(iVec) = 1 + p;
                c(iVec) = -p/2;
                d(iVec) = (p/2) * u(n,iVec - 1) + (1-p) * u(n,iVec) + (p/2) * u(n,iVec + 1);
                a(nx)   = -p;
                b(nx)   = 1 + p;
                d(nx)   = p * u(n, nx-1) + (1-p)*u(n, nx);

                %Updates the Next Column with the Calculated Temperature
                %Values
                u(n+1,:) = tdm(a,b,c,d); 

                temp2000 = L(n+1);
            else
                %Sets Values to be used for the Tri-Diagonal Matrix Method
                b(1)    = 1;
                c(1)    = 0;
                d(1)    = temp2000;
                a(iVec) = -p/2;
                b(iVec) = 1 + p;
                c(iVec) = -p/2;
                d(iVec) = (p/2) * u(n,iVec - 1) + (1-p) * u(n,iVec) + (p/2) * u(n,iVec + 1);
                a(nx)   = -p;
                b(nx)   = 1 + p;
                d(nx)   = p * u(n, nx-1) + (1-p)*u(n, nx);

                %Updates the Next Column with the Calculated Temperature
                %Values
                u(n+1,:) = tdm(a,b,c,d); 
            end
        end

    %Returns an Error Message if an Appropriate Method is not Entered
    otherwise
        error (['Undefined method: ' method]);
end

    