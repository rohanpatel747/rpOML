% rpOML - broadsearch testing script (to be removed eventually)
clear; close all; clc; format long g; rpOMLstart();



trajCase = 1;


% _________________________________________________________________________
% Inputs
if trajCase == 1
% EJS Trajectory
ei = getJulianDate('01-08-2037');
ji = getJulianDate('01-02-2039');
si = getJulianDate('01-09-2040');
in = struct;
in.verbose = true;
%in.spacing = [10; 10];
in.sequence = [
    3, ei, ei+60;
    5, ji, ji+200;
    6, si, si+400;
];
in.constrainVinf   = [sqrt(125); 16];





elseif trajCase == 2
% EVEEJ Trajectory
e1 = juliandate('01-feb-2023','dd-mmm-yyyy');
e2 = juliandate('01-sep-2023','dd-mmm-yyyy');
e3 = juliandate('01-aug-2024','dd-mmm-yyyy');
e4 = juliandate('01-nov-2027','dd-mmm-yyyy');
e5 = juliandate('01-dec-2030','dd-mmm-yyyy') - 300;
c3 = 30;
in = struct;
in.verbose = true;
in.spacing = [10;10;10;10;10];
in.sequence = [
    3, e1, e1+150;
    2, e2, e2+150;
    3, e3, e3+150;
    3, e4, e4+150;
    5, e5, e5+600;
];
in.constrainVinf   = [sqrt(c3); 16];

end


out = broadsearch(in);
broadsearch_plot(out);






















