

clear; clc; close all; format long g; rpOMLstart();
pcd = constants();
mu = pcd.Sun.mu;



encs = [2460523.50000000; 2461710.50000000];
dt = (encs(2)-encs(1))*86400;


xd = getStatePlanet(3, encs(1), 'meeus').x.';
xa = getStatePlanet(3, encs(2), 'meeus').x.';

arc1 = lambertNrev(xd, xa, dt, mu, 3);




xi1 = xd;
xi1(4:6) = arc1.vi;




options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

[~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], xd, options);
[~,state2] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], xa, options);
[~,state3] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], xi1, options);

figure()
hold on
plot3(state1(:,1),state1(:,2),state1(:,3),'k');
plot3(state2(:,1),state2(:,2),state2(:,3),'k');
plot3(state3(:,1),state3(:,2),state3(:,3),'b');
hold off
grid on; box on; axis equal;




% datetime(ans,'convertfrom','juliandate')