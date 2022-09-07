%% rpOML Examples - Plotting in the CR3BP
clear; clc; close all; format long g; rpOMLstart();


%% Overview
% This script shows the various plotting capabilities in the CR3BP.


%% Loading Constants and CR3BP Systems
c = constants();
c3= constantscr3bp();

% Define a system of interest
c3sys = c3.EarthMoon;









%% Simple CR3BP System Plot
cr3bp_plotsystemNoState(c3sys)









%% CR3BP System Plot with Specific Jacobi Constant (J) Zero-Velocity Curves
J = 3.173;
cr3bp_plotsystemNoState(c3sys,J)









%% Integrating an Initial Condition and Plotting Results
% Integrate IC below and show a full plot of the system. The ND Rot. plot 
% will include the zero-velocity curve and the Lagrange points.
% Additionally, an inertial plot is also available.

x0      = [0.12; 0; 0; 0; 3.45; 0];                 % ND Rot. ICs.
tf      = 25;                                       % ND Integ. Time
mu      = c3sys.mu;                                 % Characteristic Mass

dt      = 0.001;                                    % Output Time Step                      
options = odeset('reltol',1e-8, 'abstol',1e-8);     % Integ. Tol. Vals.

% Calling Integration (ODE113)
[t,x] = ode113(@(t,Y) eomCR3BP(t,Y,mu), 0:dt:tf, x0, options);


% Plotting Results
cr3bp_plotsystem(c3sys,t,x, ...
    'plotLagrangePts',true, ...
    'plotZVC'        ,true, ...
    'plotInertial'   ,true)









%% Showing Body Radii (ND) and a Trajectory
% Good for very very low energy trajectories that a closely bounded by one
% of the bodies.

x0      = [0.98; 0; 0; 0; 1.2; 0];                  % ND Rot. ICs.
tf      = 08;                                       % ND Integ. Time
mu      = c3sys.mu;                                 % Characteristic Mass


% Calling Integration (ODE113)
[t,x] = ode113(@(t,Y) eomCR3BP(t,Y,mu), 0:dt:tf, x0, options);


% Plotting Results
cr3bp_plotsystem(c3sys,t,x, ...
    'plotLagrangePts',true, ...
    'plotZVC'        ,false, ...
    'plotBodyRadius', true);                % <--- Plotting Lunar ND Radius
xlim([0.80 1.20]);
zlim([-0.01 0.01])
view([45 45])

















