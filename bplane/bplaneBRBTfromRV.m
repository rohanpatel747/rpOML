function [BR,BT] = bplaneBRBTfromRV(mu, x_)
%BPLANEBRBTFROMRV Returns the B-Plane BR and BT Values Only Given and RV
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]     Central Body Gravitational Parameter (km3/s2)
%       2. x_     [6x1]     State Vector w.r.t. Flyby Body (km and km/s)
%                           [x y; z; vx; vy; vz]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: [BR, BT] = 
%       1. BR     [1x1]     B*R Component (km)
%       2. BT     [1x1]     B*T Component (km)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % State r and v
    r_ = x_(1:3); r = norm(r_);
    v_ = x_(4:6); v = norm(v_);

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

end

