function cr3bp_plotsystem(sys, time, state, opts)
%CR3BP_PLOTSYSTEM Plot N.D. and Inertial System
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sys      [struct]    Structure Containing System L,V,T,mu vals.
%       2. time     [nx1]       Array of Non-Dim. Time Values
%       3. state    [nx6]       Array of Non-Dim. States
%       4. optional
%               plotBodyRadius [T/F] N.D. Radius of Primary/Secondary Bdys.
%               plotInertial   [T/F] Plot Inertial Coordinates Transform
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (none)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
       sys;
       time;
       state;
       opts.plotBodyRadius {mustBeNumericOrLogical} = false;
       opts.plotInertial   {mustBeNumericOrLogical} = true;
    end
    
    if opts.plotInertial == true
        tiledFigTotal = 3;
    else
        tiledFigTotal = 1;
    end
    
    
    J = cr3bp_computeJacobiConstant(sys,state(1,:).');
    
    

    pltBdyRad = opts.plotBodyRadius;
    if pltBdyRad
        [x,y,z] = sphere(30);
        xr1 = x*sys.rb1;    xr2 = x*sys.rb2;
        yr1 = y*sys.rb1;    yr2 = y*sys.rb2;
        zr1 = z*sys.rb1;    zr2 = z*sys.rb2;
    end
    

    for i=1:length(time)
       state1(i,1:6) = sys.b1.';
       state2(i,1:6) = sys.b2.';
    end

    out  = cr3bp_convertRotNDToInertial(sys, time, state);
    out1 = cr3bp_convertRotNDToInertial(sys, time, state1);
    out2 = cr3bp_convertRotNDToInertial(sys, time, state2);

    figure()
    ft = tiledlayout(tiledFigTotal,1,'Tilespacing','compact','Padding','compact');
    title(ft, [sys.name,' System  |  Jacobi Constant (C) = ', num2str(J)], ...
        'fontsize',14,'interpreter','latex');
    set(gcf,'color','w');

    nexttile
    hold on
    plot3(state(:,1),state(:,2),state(:,3));
    
    if pltBdyRad
        s1 = surf(xr1+sys.b1(1),yr1+sys.b1(2),zr1+sys.b1(3),'FaceAlpha',0.5);
        s2 = surf(xr2+sys.b2(1),yr2+sys.b2(2),zr2+sys.b2(3),'FaceAlpha',0.5);
        s1.EdgeColor = 'none';
        s2.EdgeColor = 'none';
        s1.FaceColor = 'r';
        s2.FaceColor = 'b';
    else
        scatter3(sys.b1(1),sys.b1(2),sys.b1(3),100,'r','filled');
        scatter3(sys.b2(1),sys.b2(2),sys.b2(3),'k','filled');
    end
    
    
    hold off
    grid on; box on; axis equal;
    legend({'Spacecraft','Primary Body','Secondary Body'},'fontsize',12,'interpreter','Latex');
    xlabel('X [ND]','fontsize',12,'interpreter','Latex');
    ylabel('Y [ND]','fontsize',12,'interpreter','Latex');
    zlabel('Z [ND]','fontsize',12,'interpreter','Latex');
    
    
    if opts.plotInertial
        
        nexttile
        hold on
        plot3(out(:,2), out(:,3), out(:,4));
        plot3(out1(:,2), out1(:,3), out1(:,4),'r');
        plot3(out2(:,2), out2(:,3), out2(:,4),'k');
        hold off
        grid on; box on; axis equal;
        legend({'Spacecraft','Primary Body','Secondary Body'},'fontsize',12,'interpreter','Latex');
        xlabel('X Inertial (km)','fontsize',12,'interpreter','Latex');
        ylabel('Y Inertial (km)','fontsize',12,'interpreter','Latex');
        zlabel('Z Inertial (km)','fontsize',12,'interpreter','Latex');

        nexttile
        hold on
        plot3(out(:,2), out(:,3), out(:,4));
        plot3(out1(:,2), out1(:,3), out1(:,4),'r');
        plot3(out2(:,2), out2(:,3), out2(:,4),'k');
        hold off
        view(45,45)
        grid on; box on; axis equal;
        legend({'Spacecraft','Primary Body','Secondary Body'},'fontsize',12,'interpreter','Latex');
        xlabel('X Inertial (km)','fontsize',12,'interpreter','Latex');
        ylabel('Y Inertial (km)','fontsize',12,'interpreter','Latex');
        zlabel('Z Inertial (km)','fontsize',12,'interpreter','Latex');
        
    end
    
end