function fig = plotBplane(mu, rbody, vinf)
%PLOTBPLANE Creates B-Plane Plot w/ Impact Disc Given Inputs
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]      Flyby Body Gravitational Parameter (km3/s2)
%       2. rbody  [1x1]      Flyby Body Radius (km)
%       3. vinf   [1x1]      Flyby Vinfinity Magnitude (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. fig    [figure]   MATLAB Figure to Plot B-Plane Coords. on
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    vi = vinf;
    rp = rbody;

    Bi = (mu/(vi^2))*((((1 + (((vi^2)*rp)/mu))^2) - 1)^(0.5));
    j=1;
    for i=0:0.05:2*pi+0.05
        bimp(j,1) = Bi*cos(i);
        bimp(j,2) = Bi*sin(i);
        j=j+1;
    end






    
    fig = figure();
    hold on;
    %scatter(targets(2,1),targets(2,2),'d','filled');
    %text(targets(2,1),targets(2,2),'FKE','FontSize',14)

    plot(bimp(:,1), bimp(:,2),'b','linewidth',1.5)
    
    hold off;
    grid on; box on; set(gcf,'color','w');
    set(gca, 'YDir','reverse')
    xlabel('B*T (km)','fontsize',16);
    ylabel('B*R (km)','fontsize',16);
    ax = gca; ax.FontSize = 16; axis equal;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    title(['B-Plane Targets | Vinf=', num2str(vi,4),'km/s'],'fontsize',18)

    figureDefaults(gca(),"figType",-2);
end