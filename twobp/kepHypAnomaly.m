function [F,ta] = kepHypAnomaly(Mh, e, tol)
%KEPHYPANOMALY Itterative Process to find F and theta given Mh and e
%   Inputs
%       Mh = Mean Hyperbolic Anomaly
%       e  = Eccentricity
%   Outputs:
%       F  = Hyperbolic Eccentric Anomaly
%       ta = True Anomaly (radians)

    if nargin < 3
        tol = 1e-8;
    end

    F = Mh;
    ratio = 1;
    while abs(ratio) > tol
        f  = e*sinh(F) - F - Mh;
        df = e*cosh(F) -1;
        ratio = f/df;
        F = F - ratio;    
    end

    ta = sqrt((e+1)/(e-1))*tanh(F/2);
    ta = 2*atan(ta);

end
