# Spacecraft Heat Shield Thermal Modelling – MATLAB
This repository contains the MATLAB implementation of a thermal analysis project for spacecraft heat shield tiles, simulating heat transfer during atmospheric re-entry. The project covers both 1D and 2D heat conduction modelling, using multiple numerical methods to solve the heat equation and evaluate thermal performance under realistic NASA re-entry conditions.

# Features

-**1D & 2D Heat Equation Solvers:** Forward Differencing, DuFort-Frankel, Backward Differencing, Crank-Nicolson.
-**Solver Techniques:** Tri-diagonal matrix algorithm for 1D implicit methods and Gauss-Seidel iteration for 2D implicit solutions.
-**Boundary Conditions:** Supports Neumann boundaries for zero-heat-flux inner surfaces.
-**Discretisation Optimisation:** Automated determination of optimal time and space steps to balance accuracy (<5% error) and computational efficiency.
-**Tile Thickness Analysis:** Calculates minimum safe heat shield thickness to protect the internal structure.
-**Damage Analysis:** Maps tile failure zones when local temperature exceeds material limits (660 °C).
-**Interactive GUI:** 1D and 2D visualisations of temperature profiles, numerical method performance, and damage progression.

# Usage

Run the main.m script to launch the GUI and simulate heat transfer for various tile locations. Select numerical methods, time/space discretization, and view interactive 1D or 2D visualizations of heat flow and damage.
