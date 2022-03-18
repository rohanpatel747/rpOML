function dYdt = eomdual2BP(t,Y, mu)
%EOMDUAL2BP 2-Body EOM [12x1] Cartesian State for 2 Bodies for given Julian Date(s) Vector
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1  Int.]  dt time step from ODE45
%       2. Y          [12x1 Dbl.]  State 1 and 2 to pass into integ. instant
%       3. mu         [1x1  Dbl.]  Central Body Gravitational Param. (km3/s2)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [12x1 Dbl.]  State Space Vector for Integration
%

    rvec1 = Y(1:3);
    vvec1 = Y(4:6);
    r    = sqrt(rvec1(1)^2+rvec1(2)^2+rvec1(3)^2) ;
    Dx   = vvec1;
    D2x  = -mu/r^3*rvec1;
    
    rvec2 = Y(7:9);
    vvec2 = Y(10:12);
    r2    = sqrt(rvec2(1)^2+rvec2(2)^2+rvec2(3)^2) ;
    Dx2   = vvec2;
    D2x2  = -mu/r2^3*rvec2;
    
    dYdt = [Dx; D2x; Dx2; D2x2];
end