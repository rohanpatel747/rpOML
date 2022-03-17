%% rpOML Examples - Broad Search Example
clear; close all; clc; format long g; rpOMLstart();


%% Overview
% This example shows how the broad search functions work to generate
% initial guess trajectories for an optimization process. The broad search
% is designed to connect Lambert arcs between planets to create
% multi-flyby trajectories. 

% For the following problem, we are interested in finding a low cost (C3)
% trajectory to Jupiter via multiple flybys of Earth and Venus. The
% sequence used for this example is: EVEEJ




%% Problem Setup
% Maximum Allowable Launch C3 (km2/s2)
c3max = 30;

% Setting Up Days to Search Within (per Encounter)
e1 = juliandate('01-feb-2023','dd-mmm-yyyy');
e2 = juliandate('01-sep-2023','dd-mmm-yyyy');
e3 = juliandate('01-aug-2024','dd-mmm-yyyy');
e4 = juliandate('01-nov-2027','dd-mmm-yyyy');
e5 = juliandate('01-dec-2030','dd-mmm-yyyy') - 300;

% Maximum Jupiter Arrival V-Infinity (km/s)
maxJVinf = 16;




%% Inputs Setup
% We start by creating a structure that the broadsearch() function accepts.

in = struct;

% The first field to create is 'sequence'. This tells the broad search
% which bodies to flyby and the lower/upper encounter dates (the lower
% dates were listed above in lines 22-26). It is ordered as such:
%           PlanetIDNumber, EncounterLowerJulian, EncounterUpperJulian; 

in.sequence = [
    3, e1, e1+150;
    2, e2, e2+150;
    3, e3, e3+150;
    3, e4, e4+150;
    5, e5, e5+600;
];

% Technically, we are done now. But additional inputs can be specified.
% Additional Inputs Include:

% V-Infinity Constraints
% We won't constrain all the encounter V-Infinity magnitudes. We only want
% the departure (due to C3) and arrival at Jupiter. Therefore, we can do
% the following:

in.constrainVinf   = [sqrt(c3max); maxJVinf];

% Let's search the grid of possible dates between our upper/lower bounds by
% a 10 day spacing. To do this, we add the following input: 

in.spacing = [10;10;10;10;10];

% If spacing is not specified, broadsearch() has a default spacing integer
% that is specified by planet (ex: Earth=10, Jupiter=50, Pluto=200 days)

% Finally we want to see what the broadsearch is computing, so turn the 
% verbose mode on:

in.verbose = true;




%% Calling broadsearch()
out = broadsearch(in);

% In the command window, the number of Lambert arcs to compute and patching
% of these legs will show up. 


%% Plotting Results
broadsearch_plot(out)

% We can now plot the results to see the various trajectories it found. A
% cost plot is also included to see how much V-Infinity discontinuity
% exists in each trajectory.



























