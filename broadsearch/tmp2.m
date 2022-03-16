

clear; clc; close all; format long g; rpOMLstart();
global mu pcd
pcd = constants();
mu = pcd.Sun.mu;






load EVEEJ_tmpoutput2.mat;




figure()
hold on
for i=1:height(out.encs)
    broadsearch_plotTraj(out.seq,out.encs(i,:))
end
hold off
grid on; box on; axis equal;
set(gcf,'color','w');








function broadsearch_plotTraj(seq,encs)

    global mu pcd

    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

    % Plot Planets in Sequence
    dtPlanets = (encs(end) - encs(1)) * 86400;
    for i=1:length(seq)
        x = getStatePlanet(seq(i), encs(1), 'meeus').x.';
        xf= getStatePlanet(seq(i), encs(i), 'meeus').x.';

        [~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dtPlanets], x, options);
        plot3(state1(:,1),state1(:,2),state1(:,3),'k');
        scatter3(xf(1), xf(2), xf(3),'r','filled');
        clear state1;
    end


    % Plot Trajectory in Sequence
    for i=1:length(seq)-1

        x1 = getStatePlanet(seq(i)  , encs(i)  , 'meeus').x.';
        x2 = getStatePlanet(seq(i+1), encs(i+1), 'meeus').x.';
        dt = (encs(i+1) - encs(i)) * 86400;

        if (seq(i) == seq(i+1)) && (dt>pcd.(getPlanetName(seq(i))).t)
            arc = lambertNrev(x1, x2, dt, mu, 3);
        else
            arc = lambert0rev(x1, x2, dt, mu);
        end

        x1 = [x1(1:3); arc.vi];
        [~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], x1, options);
        plot3(state1(:,1),state1(:,2),state1(:,3),'r');

    end
end