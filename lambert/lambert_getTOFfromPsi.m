function tof = lambert_getTOFfromPsi(psi, sv1, sv2, mu)
%LAMBERT_GETTOFFROMPSI Computes unitterated TOF given a Psi
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. psi [1x1]        Transfer Universal Variable (radians)
%       2. sv1 [6x1]or[3x1] Initial State or Position Vector (km)
%       3. sv2 [6x1]or[3x1] Final   State or Position Vector (km)
%       4. mu  [1x1]        Central Body Gravitational Parameter (km3/s2)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. tof [1x1]        Time of Flight (s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Lambert's Problem" Professor Kate Davis
%       "Lambert-Handout-3.pdf"
%       2. Bate, Muller, and White "Fundamentals of Astrodynamics"
%       3. Vallado "Fundamentals of Astrodynamics and Applications"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    r0_ = sv1(1:3); r0 = norm(r0_);
    rf_ = sv2(1:3); rf = norm(rf_);

    nu0 = atan2(r0_(2), r0_(1));
    nuf = atan2(rf_(2), rf_(1));
    dnu = nuf-nu0;
    if dnu<0; dnu = (2*pi) + dnu; end

    DM=1;
    if dnu > pi; DM=-1; end

    cdnu = dot(r0_,rf_)/(r0*rf);

    A = DM*sqrt(r0*rf*(1+cdnu));
    
    [c2, c3] = lambert_getc2c3(psi);
    
    y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));
    xi = sqrt(y/c2);
    tof = ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));

end

