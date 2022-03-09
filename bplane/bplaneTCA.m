function TCA = bplaneTCA(mu,x_)
%BPLANETCA Returns Linearized Time of Flight Given State and Grav. Param.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]      Central Body Gravitational Parameter (km3/s2)
%       2. x_     [6x1]      State Vector w.r.t. Flyby Body (km and km/s)
%                            [x y; z; vx; vy; vz]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. TCA    [1x1]      Time to go to Close Approach (s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % Find Orbital Elements and Anomalies (Ta=deg)
    s2e = conv_state2ele(x_,mu,false);
    an  = getAnomalyandDt(s2e, false);
    
    % Hyperbolic Eccentric and Mean Anomaly (rad)
    e  = s2e.e;
    a  = s2e.a;
    H  = an.ea;

    % Time to Closest Approach (TCA)
    TCA = (H - e*sinh(H)) / (sqrt(-mu/(a^3)));

end

