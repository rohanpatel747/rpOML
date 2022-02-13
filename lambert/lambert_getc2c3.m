function [c2,c3] = lambert_getc2c3(psi)
%LAMBERTC2C3 Compute the C2 and C3 Values for the Uni. Var. Lambert Alg.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. psi [1x1]        Transfer Universal Variable (radians)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. c2  [1x1]        C2 Parameter
%       2. c3  [1x1]        C3 Parameter
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Lambert's Problem" Professor Kate Davis
%       "Lambert-Handout-3.pdf"
%       2. Bate, Muller, and White "Fundamentals of Astrodynamics"
%       3. Vallado "Fundamentals of Astrodynamics and Applications"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    if psi > 1e-6       % Elliptical Case
        c2 = (1.0-cos(sqrt(psi)))/psi;
        c3 = (sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
    elseif psi < -1e-6  % Hyperbolic Case
        c2 = (1.0-cosh(sqrt(-psi)))/psi;
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/(sqrt((-psi)^3));
    else                % Parabolic Case
        c2 = 1/2;
        c3 = 1/6;
    end
end

