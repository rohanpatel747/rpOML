function out = bplanefromVi1Vi2(mu, vi1_, vi2_)
%BPLANEFROMVI1ANDVI2 Gets the B-Plane Parameters from In/Out Vinfinity
%
%   Assumptions/Warnings:
%   	1. This method and the one from bplanefromRV() will vary very
%   	   slightly due to this method ignorning the radius in the energy
%   	   equation (at infinity along a hyperbola).
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]     Central Body Gravitational Parameter (km3/s2)
%       2. vi1_   [3x1]     Incoming Vinfinity Vector (km/s)
%       3. vi2_   [3x1]     Outgoing Vinfinity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. b      [1x1]     B-Vector Magnitude (km)
%       2. b_     [3x1]     B-Vector
%       3. BR     [1x1]     B*R Component (km)
%       4. BT     [1x1]     B*T Component (km)
%       5. ta     [1x1]     Theta (Angle Between T and B) (rad.)
%       6. rp     [1x1]     Close Approach Radius (km)
%       7. d      [1x1]     Flyby Turning Angle (rad.)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. bplane_getSTR()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % Get STR Frame
    str = bplane_getSTR(vi1_);
    S_  = str.S_;
    T_  = str.T_;
    R_  = str.R_;
    vi1 = norm(vi1_);
    vi2 = norm(vi2_);
    
    % Angular Momentum Unit Vector
    H_= cross(vi1_,vi2_)/norm(cross(vi1_,vi2_));

    % Turning Angle
    d = atan2(norm(cross(vi1_,vi2_)),dot(vi1_,vi2_));
    % Alternate: d = acos(dot(vi1_,vi2_)/(vi1*vi2));
    
    % Radius of Periapsis
    rp = (mu/(vi1^2))*( (1/cos((pi-d)/2)) - 1);

    % B-Vector
    b  = rp*sqrt(1+((2*mu)/(rp*(vi1^2))));   % Miss-Distance Equation
    % Alternate: b  = (mu/(vi1^2))*((1 + ((vi1^2)*(rp/mu)))^2 - 1)^(1/2);
    B_ = cross(S_,H_);
    b_ = b*B_;
    
    
    % B-Plane Coordiantes
    BR = dot(b_,R_);
    BT = dot(b_,T_);
    ta = acos(dot(T_,B_)); if BR<0; ta=(2*pi)-ta; end
    
    
    % Outputs
    out    = struct;
    out.b  = b;
    out.b_ = b_;
    out.BR = BR;
    out.BT = BT;
    out.ta = ta;
    out.rp = rp;
    out.d  = d;

end

