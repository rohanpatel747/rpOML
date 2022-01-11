%% Examples of Using rpOML and CSPICE Mice
clear; close all; clc; format long; rpOMLstart();


% Get Planetary Constants Data
pcd = planetaryconstants();





% Example: create_state()
r_xyz = [7000; 7000; 0];
v_xyz = [0; 0; 8];
xi    = [r_xyz;v_xyz];

fX1 = create_state(xi, pcd.muEarth, 'rv');

fX2 = fX1;
fX2.ta = pi/4;
fX2.changed = 'aeiowta';
fX2 = create_state(fX2, pcd.muEarth, 'fullX', true, false);





% Example: Integrate EOMs & Plot Result
[t,stateX] = prop2bp(fX2, [0 3600], pcd.rEarth, 'given');
fX3 = create_state(stateX(end,:).', pcd.muEarth, 'rv');





% Example: Propagate using Kepler's Eqn. of a Ellipse
[E,ta] = propKepElip(fX2, 3600);

disp(fX3.ta *(pi/180))  % <-- TA from numerical integration
disp(ta)                % <-- TA from Kepler's Eqn.
disp(' ')
