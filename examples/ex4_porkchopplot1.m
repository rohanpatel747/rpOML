%% rpOML Examples - Porkchop Plot Between Earth and Mars (1)
clear; clc; format long g; rpOMLstart();

%% Overview
% This script is intended to show how to use the porkchop generation and
% plotting functions in the library.

% This example will be from Earth to Mars


%% Problem Setup
% Get Planetary Constants Data
pcd = constants();
mu  = pcd.Sun.mu;


%% Porkchop Plot Setup
% The function pkchp_porkchop() takes in a structure input with the
% following entries:

in             = struct;            % Initialize Structure
in.mu          = pcd.Sun.mu;        % Center Body Gravitational Constant
in.depBody     = 3;                 % Departure Planet (3 = Earth)
in.arrBody     = 4;                 % Arrival   Planet (4 = Mars)
in.depDate     = '04-Jun-2005';     % Earliest Departure Date
in.arrDate     = '01-Dec-2005';     % Earliest Arrival   Date
in.daysPastDep = 140;               % Latest Departure (Days since Early)
in.daysPastArr = 450;               % Latest Arrival   (Days since Early)
in.dayStep = 2;                     % Search Step Size (Days)

% Plotting Bounds
in.ctr1 = [8, 9, 10, 12, 16, 20, 25, 30, 48, 72, 100];  % Dep. C3 Contours
in.ctr2 = [1:0.5:5];                                    % Arr. Vinf Contours
in.ntof = 15;                                           % # of TOF Contours


%% Run Porkchop Plotting Algorithm and Plot Results
out  = pkchp_porkchop(in);
plt1 = pkchp_plt(out);
