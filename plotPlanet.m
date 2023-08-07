function plotPlanet(rcb, opts)
%PLOTPLANET Plots a sphere of radius r (km) on figure
    
    arguments
        rcb;
        opts.cbColor = 'c';
        opts.cbAlpha = 0.20;
        opts.cbLinesColor = [0.3, 0.3, 0.3,];
    end


    

    [X,Y] = meshgrid(-10000:1000:10000);
    Z = zeros(21);
    [a, b, c] = sphere;
    cbState = [0 0 0 rcb];
    cbColor = opts.cbLinesColor; 
    cbOBJ = mesh(a*cbState(1,4),b*cbState(1,4),c*cbState(1,4),'facecolor', opts.cbColor,'facealpha',opts.cbAlpha);
    colormap(cbColor);
    h = drawCircle(0,0,rcb);

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