%% rpOML Examples - Multi-Rev. Lambert Solution
clear; close all; clc; format long g; rpOMLstart();


%% Overview
% This script is intended to show an example of a multi-rev. Lambert
% solution between Earth and Venus




%% Problem Setup
% Get Planetary Constants Data
pcd = constants();
mu  = pcd.Sun.mu;

% Specify Departure/Arrival Bodies, associated States, and Flight Times
depbdy = 3;
arrbdy = 2;

ti = 2460545;
tf = 2460919;
dt = (tf-ti)*86400;

sv1 = getStatePlanet(depbdy, ti).x;
sv2 = getStatePlanet(arrbdy, tf).x;




%% Compute Lambert
% Specify Type III Transfer (1 Revolution, Left Branch Solution)
type = 3;


% Compute Lambert Arc
l  = lambertNrev(sv1, sv2, dt, mu, type)



%% Results
% The computed results must match the following data
% Data Provided by Professor Davis (ASEN6008)
%{

Launch (JD)     2460545					
Arrival (JD)    2460919


Planet 1 at Departure [x y z; vx vy vz]:	
130423562.1     -76679031.85	3624.816561	
14.61294123     25.56747613     -0.001503446

Planet 2 at arrival [x y z; vx vy vz]:	
19195371.67     106029328.4     348953.802	
-34.57913611	6.064190776     2.078550651



Departure Velocity [vx vy vz]:
12.76771134     22.79158874     0.090338826

Arrival Velocity [vx vy vz]:
-37.30072389    -0.176853447    -0.066693083			



Outgoing V-infinity	[vx vy vz]: 
-1.845229883	-2.77588739     0.091842272			

Incoming V-infinity [vx vy vz]:
-2.721587781	-6.241044223	-2.145243734			



C3 (km2/s2):                    11.11885913					
Magnitude V-infinity (km/s):	7.138609369	

%}


% Plot Results
figure()
rev = 1;
[psiU,psiL] = lambert_getPsiNRev(rev);
hold on
yline(psiU,'b','linewidth',1.5)
yline(psiL,'b','linewidth',1.5)
yline(l.psimin,'r','linewidth',2)
plot(l.psiV,'linewidth',1.5)
hold off
xlabel('Itterations')
ylabel('Universal Variable     \Psi (rad.)')
legend({'\Psi Upper Bound','\Psi Lower Bound', ...
    '\Psi Branch Split (min. TOF)', '\Psi @ Each Itteration'}, ...
    'fontsize',12,'location','southeast')
grid on; box on; set(gcf,'color','w');
title('\Psi versus Itterations','fontsize',14);

figure()
hold on
plot(l.tofV./86400)
hold off
xlabel('Itterations')
ylabel('Time of Flight (days)')
grid on; box on; set(gcf,'color','w');
ylim([0 500])
title('Calculated TOF versus Itterations','fontsize',14);














