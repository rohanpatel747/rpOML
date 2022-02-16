function [E,ta] = kepElipAnomaly(M, e)
%PROPKEPELIP Newton Iterative Process to find F and theta given Mh and e
%   Inputs:
%       M  = Mean Anomaly (rads)
%       e  = Eccentricity (n/a)
%   Outputs:
%       E  = Eccentric Anomaly (rads)
%       ta = True Anomaly (rads)
%
%   Info:
%       > Calling function: propKepElip(fX,dt) but can be standalone.
%       > Algorithm Reference: Curtis Orbital Mechanics for Engineers
%                              Section 3.4 Elliptical Orbits (e<1)

    err = 1e-12;                                         % Error tollerance
    
    if M<pi
        E = M + e/2;                                    % Initial Guess E for M < 180 deg.
    else
        E = M - e/2;                                    % Initial Guess E for M > 180 deg.
    end
    
    ratio = 1;
    while abs(ratio) > err
        ratio = (E - e*sin(E) - M)/(1 - e*cos(E));      % Curtis Eq3.17
        E = E - ratio;
    end
    
    ta = sqrt((1+e)/(1-e))*tan(E/2);
    ta = 2*atan(ta);
    
end

