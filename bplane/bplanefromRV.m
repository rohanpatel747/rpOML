function out = bplanefromRV(mu, x_)
%BPLANEFROMVI1ANDVI2 Gets the B-Plane Parameters from State Vector
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]     Central Body Gravitational Parameter (km3/s2)
%       2. x_     [6x1]     State Vector w.r.t. Flyby Body (km and km/s)
%                           [x y; z; vx; vy; vz]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. S_     [3x1]     S_hat Unit Vector
%       2. T_     [3x1]     T_hat Unit Vector
%       3. R_     [3x1]     R_hat Unit Vector
%       4. b      [1x1]     B-Vector Magnitude (km)
%       5. b_     [3x1]     B-Vector
%       6. BR     [1x1]     B*R Component (km)
%       7. BT     [1x1]     B*T Component (km)
%       8. ta     [1x1]     Theta (Angle Between T and B) (rad.)
%       9. rp     [1x1]     Close Approach Radius (km)
%      10. d      [1x1]     Flyby Turning Angle (rad.)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % State r and v
    r_ = x_(1:3); r = norm(r);
    v_ = x_(4:6); v = norm(v);

    % Angular Momentum and Eccentricity Unit Vectors:
    h_ = cross(r_,v_)/norm(cross(r_,v_));
    e_ = (1/mu)*(((v^2)-(mu/r))*r_ - dot(r_,v_)*v_);
    e  = norm(e_);
    
    % Half Angle Cosine and Sine
    cosp = 1/e;
    sinp = sqrt(1-(cosp^2));
    
    % STR Frame
    S_   = cosp*(e_/e) + sinp*(cross(h_,e_)/norm(cross(h_,e_)));
    T_   = cross(S_,[0;0;1])/norm(cross(S_,[0;0;1]));
    R_   = cross(S_,T_);
        
    % Energy and SMA
    E = ((v^2)/2) - (mu/r);
    a = -mu/(2*E);
      
    % B-Vector
    b = abs(a)*sqrt((e^2) - 1);
    B_= cross(S_,h_);
    b_= b*B_;
    
    % B-Plane Coordiantes
    BR = dot(b_,R_);
    BT = dot(b_,T_);
    ta = acos(dot(T_,B_)); if BR<0; ta=(2*pi)-ta; end

    
    % Rp, Vinf, and Turn Angle
    vi = sqrt(mu/abs(a));
    rp = (mu*(e-1))/(vi^2);
    d  = 2*asin(1/e);
    
    
    % Outputs
    out = struct;
    out.S_ = S_;
    out.T_ = T_;
    out.R_ = R_;
    out.b  = b;
    out.b_ = b_;
    out.BR = BR;
    out.BT = BT;
    out.ta = ta;
    out.rp = rp;
    out.d  = d;

end

