function flybyPlot2D(in)%vp_, v1_, vinf1_, vinf2_, v2_)
%FLYBYPLOT2D Plot Flyby Velocity Vector Diagram
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. vp_    [2x1]     Planet Heliocentric Velocity Vector (km/s)
%       2. v1_    [2x1]     Spacecraft Incoming Velocity Vector (km/s)
%       3. vinf1_ [2x1]     Incoming Vinfinity Vector (km/s)
%       3. vinf2_ [2x1]     Outgoing Vinfinity Vector (km/s)
%       4. v2_    [2x1]     Spacecraft Outgoing Velocity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (Figure)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    vp_    = in.vp_;
    v1_    = in.v1_;
    vinf1_ = in.vinf1_;
    vinf2_ = in.vinf2_;
    v2_    = in.v2_;

    figure()
    hold on
    quiver(0     , 0     , vp_(1)   , vp_(2)   , 'k', 'autoscale', false)
    quiver(0     , 0     , v1_(1)  , v1_(2)  , 'b', 'autoscale', false)
    quiver(vp_(1), vp_(2), vinf1_(1), vinf1_(2), 'r', 'autoscale', false)
    quiver(vp_(1), vp_(2), vinf2_(1), vinf2_(2), 'r', 'autoscale', false)
    quiver(0     , 0     , v2_(1)   , v2_(2)   , 'g', 'autoscale', false)
    hold off
    grid on; box on; set(gcf,'color','w');
    legend({'V_p','V_s_c_1','V_\infty_1','V_\infty_2','V_s_c_2'}, ...
        'location','southeast','fontsize',12)
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    xlabel('V_x');
    ylabel('V_y');

end
