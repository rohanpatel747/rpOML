function dYdt = eom2BP(t,Y, mu)
%EOM2BP 2-Body EOM [6x1]
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1 Int.]  dt time step from ODE45
%       2. Y          [6x1 Dbl.]  State to pass into integ. instant
%       3. mu         [1x1 Dbl.]  Central Body Gravitational Param. (km3/s2)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [6x1 Dbl.]  State Space Vector for Integration
%

    rvec = Y(1:3);
    vvec = Y(4:6);
    r    = sqrt(rvec(1)^2+rvec(2)^2+rvec(3)^2) ;
    Dx   = vvec;
    D2x  = -mu/r^3*rvec;
    dYdt = [Dx; D2x];
end

