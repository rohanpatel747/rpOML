function zvs = cr3bp_jacobiZVC(c3sys,J,opts) 
%CR3BP_jacobiZVC Creates In-Plane Zero Velocity Contour given Jacobi Const.
%
%   Assumptions/Warnings:
%       *** NOTE: DESIGNED FOR Z=0 (X/Y PLANAR) ZVCs ONLY ***
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure Containing System L,V,T,mu vals.
%       2. J          [nx1]      Jacobi Constant Value
%       3. opts                  Optional Plot Specific Options Containing:
%           'onlyplotcurve'[T/F]    Add curve to exisiting figure
%           'lim'          [1x1]    X/Y Axis Limits
%           'Joffset'      [1x1]    Plot till Lower Jacobi Bound (higher E)
%           'cbarfactor'   [1x1]    Multiplying Factor for Number of Legend Vals.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (figure)
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
        c3sys;
        J;
        opts.onlyplotcurve {mustBeNumericOrLogical} = false;
        opts.lim           {mustBeNumeric} = 1.400;
        opts.levels        {mustBeNumeric} = 0.001;
        opts.Joffset       {mustBeNumeric} = 0.100;
        opts.cbarfactor    {mustBeNumeric} = 1.000;
    end

    onlyplotcurve = opts.onlyplotcurve;
    lim           = opts.lim;
    levels        = opts.levels;
    Joffset       = opts.Joffset;
    cbarfactor    = opts.cbarfactor;

    x_ = -lim:0.005:lim;
    y_ = x_;
    z_ = zeros(length(x_));

    mu= c3sys.mu;

    x1 = c3sys.b1(1);   x2 = c3sys.b2(1);
    y1 = c3sys.b1(2);   y2 = c3sys.b2(2);
    z1 = c3sys.b1(3);   z2 = c3sys.b2(3);
    
    %x1 = ones(length(x_)).*x1;
    %y1 = ones(length(x_)).*y1;
    %z1 = ones(length(x_)).*z1;
    
    %x2 = ones(length(x_)).*x2;
    %y2 = ones(length(x_)).*y2;
    %z2 = ones(length(x_)).*z2;

    
    zvs = zeros(length(x_),length(x_));
    for i=1:length(x_)
        for j=1:length(y_)
            for k=1:length(z_)

                x  = x_(i);
                y  = y_(j);
                z  = 0;

                r1 = sqrt((x-x1)^2 + (y-y1)^2 + (z-z1)^2);
                r2 = sqrt((x-x2)^2 + (y-y2)^2 + (z-z2)^2);

                C  = (x^2 + y^2) + ((2*(1-mu))/r1) + ((2*mu)/r2);

                if C<=J
                    zvs(j,i) = C;
                else
                    zvs(j,i) = NaN;
                end

            end
        end
    end
    
    
    %r1 = ((x_-x1).^2 + (y_-y1).^2 + (z_-z1).^2).^(1/2);
    %r2 = ((x_-x2).^2 + (y_-y2).^2 + (z_-z2).^2).^(1/2);
    %C  = (x_.^2 + y_.^2) + ((2*(1-mu))./r1) + ((2*mu)./r2);

    if onlyplotcurve

        hold on
        contourf(x_,y_,zvs,J-Joffset:levels:J,'edgecolor','none');
        cbr = colorbar;
        set(cbr,'YTick', ...
            [J-Joffset,J-Joffset+levels:(levels*cbarfactor):J-levels,J]);
        hold off

    else

        figure()
        ft = tiledlayout(1,1,'Tilespacing','compact','Padding','compact');
        title(ft,[c3sys.name,' System  |  Jacobi Constant (C) = ', num2str(J)], ...
            'fontsize',14,'interpreter','latex');
        set(gcf,'color','w');

        nexttile
        hold on
        scatter3(c3sys.b1(1),c3sys.b1(2),c3sys.b1(3),'r','filled');
        scatter3(c3sys.b2(1),c3sys.b2(2),c3sys.b2(3),'k','filled');


        contourf(x_,y_,zvs,J-Joffset:levels:J,'edgecolor','none');
        cbr = colorbar; 
        set(cbr,'YTick', ...
            [J-Joffset,J-Joffset+levels:(levels*cbarfactor):J-levels,J]);

        hold off
        grid on; box on; axis equal;
        legend({'Primary Body','Secondary Body'},'fontsize',12,'interpreter','Latex');
        xlabel('X [ND]','fontsize',12,'interpreter','Latex');
        ylabel('Y [ND]','fontsize',12,'interpreter','Latex');
        zlabel('Z [ND]','fontsize',12,'interpreter','Latex');

        xlim([-lim lim]);
        ylim([-lim lim]);

    end

end