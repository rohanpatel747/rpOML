function cr3bp_plotsystem(sys, time, state, opts)
%CR3BP_PLOTSYSTEM Plot N.D. and Inertial System
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sys      [struct]    Structure Containing System L,V,T,mu vals.
%       2. time     [nx1]       Array of Non-Dim. Time Values
%       3. state    [nx6]       Array of Non-Dim. States
%       4. optional
%               plotBodyRadius  [T/F] N.D. Radius of Primary/Secondary Bdys.
%               plotInertial    [T/F] Plot Inertial Coordinates Transform
%               plotLagrangePts [T/F] Plots L1-5 from System
%               plotZVC         [T/F] Plots Zero-Velocity Contour
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (none)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
       sys;
       time;
       state;
       opts.plotLagrangePts {mustBeNumericOrLogical} = true;
       opts.plotZVC         {mustBeNumericOrLogical} = false;
       opts.plotBodyRadius  {mustBeNumericOrLogical} = false;
       opts.plotInertial    {mustBeNumericOrLogical} = false;
    end
    
    if opts.plotInertial == true
        tiledFigTotal = 2;
    else
        tiledFigTotal = 1;
    end
    
    
    % Compute System Jacobi Constant From Initial State
    J = cr3bp_computeJacobiConstant(sys,state(1,:).');
    

    % Plots ND Rot Radius of Primary and Secondary Bodies
    if opts.plotBodyRadius
        [x,y,z] = sphere(30);
        xr1 = x*sys.rb1;    xr2 = x*sys.rb2;
        yr1 = y*sys.rb1;    yr2 = y*sys.rb2;
        zr1 = z*sys.rb1;    zr2 = z*sys.rb2;
    end

    % Find Dim. Inertial Coordinates and Time 
    if opts.plotInertial
        for i=1:length(time)
           state1(i,1:6) = sys.b1.';
           state2(i,1:6) = sys.b2.';
        end
    
        out  = cr3bp_convertRotNDToInertial(sys, time, state);
        out1 = cr3bp_convertRotNDToInertial(sys, time, state1);
        out2 = cr3bp_convertRotNDToInertial(sys, time, state2);
    end






    figure()
    ft = tiledlayout(tiledFigTotal,1,'Tilespacing','compact','Padding','compact');
    title(ft, [sys.name,' System  |  Jacobi Constant (C) = ', num2str(J,8)], ...
        'fontsize',14,'interpreter','latex');
    set(gcf,'color','w');

    nexttile
    lgndNum = 1;
    hold on


    % Plot Rot. Spacecraft State
    plot3(state(:,1),state(:,2),state(:,3));
    lgnd = {'Spacecraft'};
    

    % Plots ND Rot Radius of Primary and Secondary Bodies
    if opts.plotBodyRadius
        s1 = surf(xr1+sys.b1(1),yr1+sys.b1(2),zr1+sys.b1(3),'FaceAlpha',0.5);
        s2 = surf(xr2+sys.b2(1),yr2+sys.b2(2),zr2+sys.b2(3),'FaceAlpha',0.5);
        s1.EdgeColor = 'none';
        s2.EdgeColor = 'none';
        s1.FaceColor = [0.4660 0.6740 0.1880];
        s2.FaceColor = 'k';
    else
        scatter3(sys.b1(1),sys.b1(2),sys.b1(3),100,[0.4660 0.6740 0.1880],'filled');
        scatter3(sys.b2(1),sys.b2(2),sys.b2(3),'k','filled');
    end
    lgnd{1,lgndNum+1} = 'Primary Body';
    lgnd{1,lgndNum+2} = 'Secondary Body';
    lgndNum = lgndNum+2;


    % Plots Zero Velocity Curves for Initial State Jacobi Constant
    if opts.plotZVC
        zvs = cr3bp_jacobiZVC(sys,J,'onlyplotcurve',true,'simplePlot',true);
        lgnd{1,lgndNum+1} = 'ZVC';
        lgndNum = lgndNum + 1;
    end


    % Plots Lagrange Points of the System
    if opts.plotLagrangePts
        sys = cr3bp_computeLagrangePoints(sys);
        Lpts  = [sys.L1.';
            sys.L2.';
            sys.L3.';
            sys.L4.';
            sys.L5.'];
        hold on
        scatter(Lpts(:,1),Lpts(:,2),'d','r','filled');
        hold off
        lgnd{1,lgndNum+1} = 'Lagrange Points';
        lgndNum = lgndNum + 1;
    end

    
    hold off
    grid on; box on; axis equal;
    legend(lgnd,'fontsize',12,'interpreter','latex');
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
        %view(45,45)
        grid on; box on; axis equal;
        legend({'Spacecraft','Primary Body','Secondary Body'},'fontsize',12,'interpreter','Latex');
        xlabel('X Inertial (km)','fontsize',12,'interpreter','Latex');
        ylabel('Y Inertial (km)','fontsize',12,'interpreter','Latex');
        zlabel('Z Inertial (km)','fontsize',12,'interpreter','Latex');
        
    end
    
end