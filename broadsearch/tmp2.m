

clear; clc; close all; format long g; rpOMLstart();
pcd = constants();
mu = pcd.Sun.mu;






load EVEEJ_tmpoutput.mat;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;


figure()
hold on

% Plot Planets in Sequence
dtPlanets = (out.encs(end) - out.encs(1)) * 86400;
for i=1:length(out.seq)
    x = getStatePlanet(out.seq(i), out.encs(1), 'meeus').x.';
    xf= getStatePlanet(out.seq(i), out.encs(i), 'meeus').x.';
    
    [~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dtPlanets], x, options);
    plot3(state1(:,1),state1(:,2),state1(:,3),'k');
    scatter3(xf(1), xf(2), xf(3),'r','filled');
    clear state1;
end


% Plot Trajectory in Sequence
for i=1:length(out.seq)-1
    
    x1 = getStatePlanet(out.seq(i)  , out.encs(i)  , 'meeus').x.';
    x2 = getStatePlanet(out.seq(i+1), out.encs(i+1), 'meeus').x.';
    dt = (out.encs(i+1) - out.encs(i)) * 86400;

    if (out.seq(i) == out.seq(i+1)) && (dt>pcd.(getPlanetName(out.seq(i))).t)
        arc = lambertNrev(x1, x2, dt, mu, 3);
    else
        arc = lambert0rev(x1, x2, dt, mu);
    end
    
    x1 = [x1(1:3); arc.vi];
    [~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], x1, options);
    plot3(state1(:,1),state1(:,2),state1(:,3),'r');

end







hold off
grid on; box on; axis equal;
set(gcf,'color','w');




%{
encs = [2460523.50000000; 2461710.50000000];
dt = (encs(2)-encs(1))*86400;


xd = getStatePlanet(3, encs(1), 'meeus').x.';
xa = getStatePlanet(3, encs(2), 'meeus').x.';

arc1 = lambertNrev(xd, xa, dt, mu, 3);


xi1 = xd;
xi1(4:6) = arc1.vi;


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
%}