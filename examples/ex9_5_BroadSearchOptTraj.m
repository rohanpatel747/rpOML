%% rpOML Examples - Broad Search Optimization 
clear; close all; clc; format long g;  beep off;

%% Overview
% This example shows you how to take certain trajectories from the
% BroadSearch algoirthm set in this library and convert them into to
% useable interplanetary trajectories in the patched conics 2-Body model.
% By usable, I mean where the incoming and outgoing V-inf. discontinuities
% are minimized and flyby properties are valid.



% For the following problem, consider a simple Earth(3)-Mars(4)-Jupiter(5)
% trajectory.




%% Inputs Setup
% The optimizer function needs to know the sequence, lower bound of
% enconter dates, guess for encounter dates, and upper bound of enconter
% dates. This is seen in the setup variables "seq", "e1..e3".

% The optimizer function broadsearch_optTraj() takes a structure input with
% the fields:
%   "sequence"  for an array of enconters and guesses,
%   "c3max"     for the maximum launch C3 (km3/s2)
%   "rminFB"    for the minimum flyby altitude for every encounter (km)
%   "optType"   for the optimization cost function.
%                   'minTOF' = minimum flight time
%                   'minC3'  = minimum launch C3
%   "ephemType" for the ephemeris type
%                   'meeus'  = Meeus Algorithm Solution
%                   'ephem'  = CSPICE DE4XX Ephemeris Solution

seq = [3; 4; 5]; 
e1  = 2462868;
e2  = e1+150;
e3  = e2+700;

in = struct;
in.sequence = [
    seq(1), e1-500, e1, e1+500;
    seq(2), e2-500, e2, e2+500;
    seq(3), e3-1000, e3, e3+1000;
];

in.c3max     = 100;
in.rminFb    = 600;
in.optType   = 'minTOF';
in.ephemType = 'ephem';



%% Optimize Enconters 
% We can specify if we want parallel computing to be used as shown in the
% function call optional input 'useParallel'.

% If we are using SPICE data and want to specify custom .bsp files, put
% them in the cspiceSPK cell array. We won't for this example so its empty
% (or we don't even need to specify it as an input to the function).

in  = broadsearch_optTraj(in, 'useParallel', false, 'cspiceSPK', {});



%% Results and Plotting
% Now that our trajectory is optimized, we can view the specific results
% with the following function. The input structure we created is appended
% with results from optimization (from the previous function), so now we
% pass this to the results function. Encounter data, flyby properties, 
% launch stats, and optimization results are all printed to the console.

out = broadsearch_optTrajResults(in);


















