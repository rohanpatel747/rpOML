function [E,ta] = propKepElip(fX, dt)
%PROPKEPELIP Newton Itterative Process to find M and Ta given Me and e
%   Inputs
%       fX : Full State Structure Using:
%               T   = Period (s)
%               e   = Eccentricity (n/a)
%               dti = Initial State's Time Past Periapsis (s)
%       dt : Propagation time from initial time past periapsis (s)
%   Outputs:
%       E  : Eccentric Anomaly after dt (rad)
%       ta : True Anomaly after dt (rad)
%
%
%   Info:
%       > Calls: kepElipAnomaly(Me, e) which solves using Newton's Method.
%       > Algorithm Reference: Curtis Orbital Mechanics for Engineers
%                              Section 3.4 Elliptical Orbits (e<1)
%

    T  = fX.T;                       % Initial orbit's period (s)
    e  = fX.e;                       % Initial orbit's eccentricity (n/a)
    dti= fX.dt;                      % Initial state's time since periapsis (s)

    if e >= 1
        disp('Cannot use this function. propKepElip is only for eliptical/circular trajectories')
        E = NaN; ta = NaN;
    else
        dt = dt + dti;                   % New time past periapsis (s)

        Me = (2*pi)*(dt/T);              % Mean Anomaly at new dt (rad)

        [E,ta] = kepElipAnomaly(Me, e);   % Itterative method to find new E and ta
    end

end

