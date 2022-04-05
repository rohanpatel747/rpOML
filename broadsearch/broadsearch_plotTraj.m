function broadsearch_plotTraj(seq,encs,mu,pcd,plotColor)
%BROADSEARCH_PLOTTRAJ Intermediate Integration + Plotting of Trajectories
%
%   Assumptions/Warnings:
%   	[NOTE: This function is not meant to be standalone]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. broadsearch_plotResults()
%       2. getStatePlanet()
%       3. getPlanetName()
%       4. bplanefromVi1Vi2()
%       5. eom2BP()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    

    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

    % Plot Planets in Sequence
    dtPlanets = (encs(end) - encs(1)) * 86400;
    for i=1:length(seq)
        x = getStatePlanet(seq(i), encs(1), 'meeus').x.';
        xf= getStatePlanet(seq(i), encs(i), 'meeus').x.';

        [~,state1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dtPlanets], x, options);
        plot3(state1(:,1),state1(:,2),state1(:,3),'k');
        scatter3(xf(1), xf(2), xf(3),plotColor,'filled');
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
        plot3(state1(:,1),state1(:,2),state1(:,3),plotColor);

    end
end