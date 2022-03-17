

clear; clc; close all; format long g; rpOMLstart();
global mu pcd
pcd = constants();
mu = pcd.Sun.mu;





%load EJP_NH.mat;
%load EVEEJ_tmpoutput3.mat;
load EH_EJS_Lab6.mat;
Jmax = 5;

out.encs = out.encs(3,:);


gaseqname = broadsearch_sequencename(out.seq);

if true
figure()
hold on
for i=1:height(out.encs)
    broadsearch_plotTraj(out.seq,out.encs(i,:))
end
hold off
title([gaseqname, ' Trajectory'],'fontsize',12,'interpreter','latex');
xlabel('X EC','fontsize',12,'interpreter','latex');
ylabel('Y EC','fontsize',12,'interpreter','latex');
zlabel('Z EC','fontsize',12,'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
grid on; box on;
set(gcf,'color','w'); axis equal;
end

if true
figure()
hold on
scatter([1:height(out.cost)],out.costsum,'r','filled');
lgnd = {'Total'};

for i=1:width(out.cost)
    scatter([1:height(out.cost)],out.cost(:,i),'x');
    lgnd{1,i+1} = [gaseqname(i+1),'GA']; %[getPlanetName(out.seq(i)),' GA'];
end
hold off
if max(out.costsum)+1 < Jmax
    ylim([-0.0 max(out.costsum)+1]);
else
    ylim([-0.0 Jmax]);
end
xticks([1:1:height(out.cost)]);
title([gaseqname, ' Sequence Cost'],'fontsize',12,'interpreter','latex');
xlabel('Trajectory Number','fontsize',12,'interpreter','latex');
ylabel('Velocity Discontinuity at Flyby (Cost) (km/s)','fontsize',12,'interpreter','latex');
legend(lgnd,'fontsize',12,'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
grid on; box on;
set(gcf,'color','w');
end












function b = broadsearch_sequencename(seq)

    b = [];
    for i=1:length(seq)
        a = getPlanetName(seq(i));
        a = a(1);
        b = [b,a];
    end
end



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