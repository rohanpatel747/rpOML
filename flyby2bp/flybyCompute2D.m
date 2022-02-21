function out = flybyCompute2D(mu, vp_, vinf_, rp, opt)
%FLYBYCOMPUTE2D Planar Post-Flyby Vinf and V2 Computation
%
%   Assumptions/Warnings:
%   	1. Turn angle direction not checked yet
%       2. Default Flyby Type: LEADING ('leading',true) in input args. 
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]     Central Body Gravitational Parameter (km3/s2)
%       2. vp_    [2x1]     Planet Heliocentric Velocity Vector (km/s)
%       3. vi1_   [2x1]     Incoming Vinfinity Vector to Planet (km/s)
%       4. rp     [1x1]     Close Approach RADIUS (km)
%       5. 'leading', True  Default Flyby Type (set to False for Trailing)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. e      [1x1]     Eccentricty
%       2. d      [1x1]     Turn Angle (rad.)
%       3. md     [1x1]     Miss-Distance (B*T) (km)
%       4. vperi  [1x1]     Close Approach Velocity (km/s)
%       5. h      [1x1]     Flyby Traj. Angular Momentum (km2/s)
%       6. v1_    [2x1]     Incoming Heliocentric Velocity Vector (km/s)
%       7. v2_    [2x1]     Outgoing Heliocentric Velocity Vector (km/s)
%       8. vi1_   [2x1]     Incoming Vinfinity Vector to Planet (km/s)
%       9. vi2_   [2x1]     Outgoing Vinfinity Vector w.r.t Planet (km/s)
%      10. vp_    [2x1]     Planet Heliocentric Velocity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
        mu;
        vp_;
        vinf_;
        rp;
        opt.leading {mustBeNumericOrLogical} = true;
    end

    vinf = norm(vinf_);
    e      = 1 + ((rp*(norm(vinf)^2)) / mu);
    d      = 2*asin(1/e);
    
    if opt.leading==false; d=-d; end
    
    md     = rp*sqrt(1+((2*mu)/(rp*(vinf^2))));
    vperi  = sqrt((vinf^2) + ((2*mu)/rp));
    h      = rp*vperi;
    %v2    = (vinf^2) + (vp^2) - 2*vinf*vp*cos(b2);
    rotMat = [cos(-d) -sin(-d); sin(-d)  cos(-d)];
    vinf2_ = rotMat*vinf_;
    v2_    = vp_ + vinf2_;
    
    
    out        = struct;
    out.e      = e;
    out.d      = d;
    out.md     = md;
    out.vperi  = vperi;
    out.h      = h;
    out.v1_    = vp_ + vinf_;
    out.v2_    = v2_;
    out.vi1_   = vinf_;
    out.vi2_   = vinf2_;
    out.vp_    = vp_;
end

