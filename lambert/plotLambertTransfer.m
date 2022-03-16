function [dt, stateX3] = plotLambertTransfer(mu, x1_, x2_, dv1_, dt)

    x1_t = x1_;
    x1_t(4) = x1_(4) + dv1_(1);
    x1_t(5) = x1_(5) + dv1_(2);
    x1_t(6) = x1_(6) + dv1_(3);
    

    fX1 = create_state(x1_, mu, 'rv', false, false, true);          % Initial
    fX2 = create_state(x2_, mu, 'rv', false, false, true);          % Final
    fX3 = create_state(x1_t, mu, 'rv', false, false, true);         % Transfer
    [~,stateX1] = plotTraj(fX1, [0 dt], 1000, 'none');              % Integ. Initial  [0 dt]
    [~,stateX2] = plotTraj(fX2, [dt 0], 1000, 'none');              % Integ. Final    [dt 0]
    [~,stateX11] = plotTraj(fX1, [0 fX1.T], 1000, 'none');          % Integ. Initial  [0 T]
    [~,stateX22] = plotTraj(fX2, [fX2.T 0], 1000, 'none');          % Integ. Final    [0 T]
    [~,stateX3] = plotTraj(fX3, [0 dt], 1000, 'none');              % Integ. Transfer [0 dt]
    [~,stateX33] = plotTraj(fX3, [0 fX3.T], 1000, 'none');          % Integ. Transfer [0 T]

    
    %figure
    hold on
    % Original/Final Orbit
    plot3(stateX11(:,1),stateX11(:,2),stateX11(:,3),'-.k')
    plot3(stateX22(:,1),stateX22(:,2),stateX22(:,3),'-.k')
    plot3(stateX33(:,1),stateX33(:,2),stateX33(:,3),'-.k')

    % Original/Final Traveled Orbit
    plot3(stateX1(:,1),stateX1(:,2),stateX1(:,3),'b','linewidth',1.5)
    plot3(stateX2(:,1),stateX2(:,2),stateX2(:,3),'r','linewidth',1.5)

    % r1_ and r2_
    quiver3(0,0,0,x1_(1),x1_(2),x1_(3),'autoscale','off','color','b')
    quiver3(0,0,0,x2_(1),x2_(2),x2_(3),'autoscale','off','color','r')
    
    %quiver3(0,0,0,fX3.e_(1)*1e9,fX3.e_(2)*1e9,fX3.e_(3)*1e9,'autoscale','off','color','g')
    %quiver3(0,0,0,-fX3.e_(1)*1e9,-fX3.e_(2)*1e9,-fX3.e_(3)*1e9,'autoscale','off','color','g')

    % Transfer Orbit
    plot3(stateX3(:,1),stateX3(:,2),stateX3(:,3),'g','linewidth',1.5)

    hold off
    grid on; box on; axis equal;
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    set(gcf,'color','w');

    
    
    
    
    
    
    
    
    
    
    
    
    function [t,stateX] = plotTraj(fX, dt, rcb, pltType)
    %PLOTTRAJ Plot Trajectory Integrated from ti to tf
    %   Detailed explanation goes here
    tic

        % Auxilary Inputs (Not Required)
        auxIn = struct;
        auxIn.cbColor = 'c';
        auxIn.cbAlpha = 0.20;
        auxIn.cbLinesColor = [0.3, 0.3, 0.3,];



        % Integrate Trajectory
        mu = fX.mu;
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;
        [t,stateX] = ode45(@rates, dt, fX.xi, options);

        % Post Processing
        rmag = zeros(length(stateX),1);
        vmag = zeros(length(stateX),1);
        for i=1:length(stateX)
            rmag(i,1) = norm(stateX(i,1:3));
            vmag(i,1) = norm(stateX(i,4:6));
        end




        % TODO
        % 1. add orbital elements to plot
        % 2. varargin auxilary inputs structure




        % Create Figure
        if strcmp(pltType,'none')

        elseif strcmp(pltType,'given')
            hold on
            [X,Y] = meshgrid(-10000:1000:10000);
            Z = zeros(21);
            [a, b, c] = sphere;
            cbState = [0 0 0 rcb];
            cbColor = auxIn.cbLinesColor; 
            cbOBJ = mesh(a*cbState(1,4),b*cbState(1,4),c*cbState(1,4),'facecolor', auxIn.cbColor,'facealpha',auxIn.cbAlpha);
            colormap(cbColor);

            plot3(stateX(:,1),stateX(:,2),stateX(:,3),'linewidth',2.5)
            scale = fX.rp+100;
            if strcmp(pltType,'elements')
                e__ = fX.e_/fX.e;
                quiver3(0,0,0,e__(1)*scale,e__(2)*scale,e__(3)*scale,'Autoscale','off','linewidth',2)
                quiver3(0,0,0,fX.N_(1),fX.N_(2),fX.N_(3),'Autoscale','off','linewidth',2,'color','g')
                legend({'Central Body','Orbit','Ecc. Vector','Nodes Vector'},'fontsize',14,'location','southeast')
            end
            h = drawCircle(0,0,rcb);
            hold off
            xlabel('x (km)')
            ylabel('y (km)')
            axis equal
            grid on
            setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
            set(gcf,'color','w')
        else
            figure
            [X,Y] = meshgrid(-10000:1000:10000);
            Z = zeros(21);
            [a, b, c] = sphere;
            cbState = [0 0 0 rcb];
            cbColor = auxIn.cbLinesColor; 
            cbOBJ = mesh(a*cbState(1,4),b*cbState(1,4),c*cbState(1,4),'facecolor', auxIn.cbColor,'facealpha',auxIn.cbAlpha);
            colormap(cbColor);
            hold on
            plot3(stateX(:,1),stateX(:,2),stateX(:,3),'linewidth',2.5)
            scale = fX.rp+100;
            if strcmp(pltType,'elements')
                e__ = fX.e_/fX.e;
                quiver3(0,0,0,e__(1)*scale,e__(2)*scale,e__(3)*scale,'Autoscale','off','linewidth',2)
                quiver3(0,0,0,fX.N_(1),fX.N_(2),fX.N_(3),'Autoscale','off','linewidth',2,'color','g')
                legend({'Central Body','Orbit','Ecc. Vector','Nodes Vector'},'fontsize',14,'location','southeast')
            end
            h = drawCircle(0,0,rcb);
            hold off
            xlabel('x (km)')
            ylabel('y (km)')
            axis equal
            grid on
            setAxes3DPanAndZoomStyle(zoom(gca),gca,'camera');
            set(gcf,'color','w')
        end

    toc

        % 2BP EOMs
        function dYdt = rates(t,Y)
            rvec    = Y(1:3);
            vvec    = Y(4:6);
            r = sqrt(rvec(1)^2+rvec(2)^2+rvec(3)^2) ;
            Dx   = vvec;
            D2x  = -mu/r^3*rvec;
            dYdt = [Dx; D2x];
        end

        % Draw Circle
        function h = drawCircle(x,y,r)
            hold on
            th = 0:pi/50:2*pi;
            xunit = r * cos(th) + x;
            yunit = r * sin(th) + y;
            h = plot(xunit, yunit);
            hold off
        end

    end

    
end
