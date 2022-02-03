% rpOML Examples - NBody Propagation

clear; close all; clc; format long g; rpOMLstart();









% TODO:
%{
    1. figure out why 2bp doesnt land exactly on mars
    2. make the backoffset more elegant


    3. Constants() add meeus data and do an "if requested for only one body
    kinda of input and output"
%}












% Get Planetary Constants Data
pcd = constants();
% Override Constants
pcd.muSun   = 1.327e11;
pcd.muEarth = 3.968e5;
pcd.muMars  = 4.305e4;
pcd.rMars = 3397.2;


rE = pcd.Earth.r;
aE = pcd.Earth.sma;
rM = pcd.Mars.r;
aM = pcd.Mars.r;
mu = pcd.Sun.mu;
muE= pcd.Earth.mu;
muM= pcd.Mars.mu;
rS = 696347.055;






%function [t, state] = nbp(xi, ti, dt, mu, [3,4], cbname)
% NBP Propagates N-Body EOMs
%
%   Assumptions/Warnings:
%   	1. 
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. xi [6x1] Initial Inertial State [x y z vx vy vz] in km and km/s
%       2. ti [1x1] Initial Julian Date
%       3. dt [1x1] Propagation Time (seconds)
%       4. mu [1x1] Central Body Gravitational Parameter (km3/s2)
%       5. bd [nx1] Additional Gravitational Bodies
%
%
%   TODO:
%       1. ability to use SPK data to find body pos.
%               


%if useSPK
%    ti = cspice_str2et(['JD',jdate])
%
%
%end


ti = 'Jul 01, 2020';
eti = cspice_str2et({ti});


%{

USING NAIF SPK DATA
    -> Write sep. NBP EOMs function
    -> Use ET for time calculation (bc its in seconds)
       mice_spkezr(naifID, et, 'J2000', 'NONE', cbname);

       where:
        1. naifID = spkID of planet
        2. et = et_i + t    ;where:  t    = integration dt time
                                     et_i = ephemeris time of initial epoch
        3. cbname = Central Body Name 'SOLAR SYSTEM BARYCENTER'
%}





% Lambert Sample Problem
depBdy = '3';
arrBdy = '4';
etp    = cspice_str2et( {'Jul 30, 2020', 'Feb 18, 2021'} );

[ctr_bdy] = mice_bodc2n(0);
tof = etp(2) - etp(1);
sv1 = mice_spkezr(depBdy, etp(1), 'J2000', 'NONE', ctr_bdy.name ).state;
sv2 = mice_spkezr(arrBdy, etp(2), 'J2000', 'NONE', ctr_bdy.name ).state;
[vito, vfto] = glambert(mu, sv1, sv2, tof, 0);









% Verification via 2BP & NBP Integration
noPts = 1000;
xi = [sv1(1:3); vito];

% Get Planet Positions Over Time Span
tspan = (0:noPts-1) * ( etp(2) - etp(1) )/noPts + etp(1);
for i=1:length(tspan)
   xb1(i,1:3) =  mice_spkezr(depBdy, tspan(i), 'J2000', 'NONE', ctr_bdy.name ).state(1:3).';
   xb2(i,1:3) =  mice_spkezr(arrBdy, tspan(i), 'J2000', 'NONE', ctr_bdy.name ).state(1:3).';
end



% Numerical Integration with ODE45
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12) ;

% 2-Body EOM Integration
[t2bp,s2bp] = ode45(@(t,Y) eom2BP(t,Y,mu), [0:noPts:tof], xi, options);

% N-Body EOM Integration
backOffset  = 86400*1.5;

pList       = [3; 4; 5; 6];
muList      = [pcd.Earth.mu; pcd.Mars.mu; pcd.Jupiter.mu; pcd.Saturn.mu];
eti         = etp(1);
cbname      = ctr_bdy.name;

[tnbp,snbp] = ode45(@(t,Y) eomNBP(t,Y, mu, pList, muList, eti, cbname, backOffset), ...
    [0:noPts:tof], xi, options);










% Plotting
if true
    
figure()
hold on
plot3(xb1(:,1), xb1(:,2), xb1(:,3), 'k')
plot3(xb2(:,1), xb2(:,2), xb2(:,3), 'k')
plot3(s2bp(:,1), s2bp(:,2), s2bp(:,3), 'r')
plot3(snbp(:,1), snbp(:,2), snbp(:,3), 'b')

hold off
legend({'Earth','Mars', '2BP Traj.', 'NBP Traj.'}, ...
    'fontsize', 12, 'location', 'southwest');
xlabel('X (km)','fontsize', 12);
ylabel('Y (km)','fontsize', 12);
zlabel('Z (km)','fontsize', 12);
%title('Comparing 2BP versus NBP w/ Earth and Mars Gravitational Effects', ...
%    'fontsize', 12)
grid on; box on; axis equal;
set(gcf,'color','w');



% Error Plot
t = tnbp./86400;    % Time in days (/86400)

for i=1:length(t)
   norm_err(i,1) = abs(norm(snbp(i,1:3)) - norm(s2bp(i,1:3)));
   norm_err(i,2) = abs(norm(snbp(i,4:6)) - norm(s2bp(i,4:6)));
end

figure()
f1t = tiledlayout(2,1,'Tilespacing','compact','Padding','compact');
title(f1t,'Difference (NBP-2BP) ','fontsize',14,'fontweight','bold');
set(gcf,'color','w');

nexttile
hold on
plot(t,snbp(:,1)-s2bp(:,1),'r');
plot(t,snbp(:,2)-s2bp(:,2),'b');
plot(t,snbp(:,3)-s2bp(:,3),'g');
plot(t,norm_err(:,1), 'c')
hold off
legend({'X-Position','Y-Position','Z-Position','Abs. Norm. Error'}, ...
    'fontsize', 12, 'location', 'southwest');
ylabel('Position','fontsize', 12);
grid on; box on;
set(gcf,'color','w');
xlim([0 210]);

nexttile
hold on
plot(t,snbp(:,4)-s2bp(:,4),'r');
plot(t,snbp(:,5)-s2bp(:,5),'b');
plot(t,snbp(:,6)-s2bp(:,6),'g');
plot(t,norm_err(:,2), 'c')
hold off
legend({'X-Velocity','Y-Velocity','Z-Velocity','Abs. Norm. Error'}, ...
    'fontsize', 12, 'location', 'southwest');
xlabel('Time (days)','fontsize', 12);
ylabel('Velocity (km/s)','fontsize', 12);
grid on; box on;
set(gcf,'color','w');
xlim([0 210]);
end



















% 2BP EOMs
function dYdt = TBPrates(t,Y)
    mu   = 1.327000000000000e+11;
    rvec = Y(1:3);
    vvec = Y(4:6);
    r    = sqrt(rvec(1)^2+rvec(2)^2+rvec(3)^2) ;
    Dx   = vvec;
    D2x  = -mu/r^3*rvec;
    dYdt = [Dx; D2x];
end





















