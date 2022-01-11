function pK = propKepHyp(fX,dt)
%PROPKEP HYPERBOLA Given initial state and some time duration, find final state.
%   (ONLY WORKS FOR HYPERBOLIC TRAJECTORIES)
%
%   Input:
%       1. fX  - Full State Structure (from create_state)
%       2. dt  - Propagation Time duration in seconds
%   Output Structure (pK) with fields:
%       1. ta  - True Anomaly (rad.)
%       2. taf - True Anomaly Inf. (rad.)
%       3. xf  - (6x1) State in km and km/s

    % Initial State & Orbit Properties
    x    = fX.xi;
    mu   = fX.mu;
    a    = fX.a;
    e    = fX.e;
    h    = fX.h;
    inc  = fX.i * (pi/180);
    ta   = fX.ta * (pi/180);
    raan = fX.o * (pi/180);
    aop  = fX.w * (pi/180);


    if e <= 1
        disp('Cannot use this function. propKep is only for hyperbolic trajectories')
        pk = NaN;
    else
        % Find True Anomaly Inf.
        taf = acos(-1/e);

        % Hyperbolic Eccentric Anomaly given Initial Theta (Need for dt correction)
        F = log((sqrt(e+1) + sqrt(e-1) * tan(ta/2))/(sqrt(e+1) - sqrt(e-1) * tan(ta/2)));
        Mh = e*sinh(F) - F;
        ti = (h^3/(mu^2))*(1/((e^2 - 1)^(3/2)))*Mh; % b/c state isn't at TA=0, need to find TA=0 -> TA=input state vec.
        dt = dt+ti;

        % Find Hyperbolic Anomaly and True Anomaly after dt
        Mh = ((mu^2)/(h^3))*((e^2 - 1)^(3/2))*dt;
        [~, ta] = kepHypAnomaly(Mh, e);

        % Convert from Keplerian Elements to Perifocal Frame
        rperi = (((h^2)/mu)*(1/(1+e*cos(ta)))).*[cos(ta);sin(ta);0];
        vperi = (mu/h)*[-sin(ta);e+cos(ta);0];

        % Convert from Perifocal Frame to Initial Frame (R313 DCM)
        QxGP = [cos(aop) sin(aop) 0; -sin(aop) cos(aop) 0; 0 0 1]* ...
            [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)]* ...
            [cos(raan) sin(raan) 0; -sin(raan) cos(raan) 0; 0 0 1];

        % Final State Vector after dt
        rMJ2KEC = QxGP'*rperi;
        vMJ2KEC = QxGP'*vperi;

        % Output
        pK = struct;
        pK.ta = ta;
        pK.taf = taf;
        pK.xf = [rMJ2KEC;vMJ2KEC];
    end
end