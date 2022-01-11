function [t,stateX] = propTraj(fX, dt, rcb, pltType)
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

