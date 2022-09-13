function out = cr3bp_eigValsLinSys(c3sys, xyz)
%CR3BP_EIGVALSLINSYS Finds the Linearized Dynamics' Eigen-Values.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure containing at least:  c3sys.mu
%                                where mu is the characteristic mass ratio
%       2. xyz        [3x1]      Position Vector [x; y; z]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: structure 'out' containing fields:
%
%       1. .eigVals   [6x1]      Set of Eigen Values:
%                                   [1:4,1] are the In-Plane Modes
%                                   [5:6,1] are the Out-of-Plane Modes
%       2. .U         [1x1]      Pseudo-Potential Function
%       3. .jU        [3x1]      1st Partial of Pseudo-Potential Function
%       4. .jUjU      [6x6]      2nd Partial of Pseudo-Potential Function
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % jUjU @ x,y,z
    out  = cr3bp_pseudoPotential(c3sys,xyz);
    Uxx = out.jUjU(1,1);     Uxy = out.jUjU(1,2);     Uxz = out.jUjU(1,3);
    Uyx = out.jUjU(2,1);     Uyy = out.jUjU(2,2);     Uyz = out.jUjU(2,3);
    Uzx = out.jUjU(3,1);     Uzy = out.jUjU(3,2);     Uzz = out.jUjU(3,3);
    
    % In-Plane Modes:
    a = 1;
    b = (4-Uxx-Uyy);
    c = (Uxx*Uyy - (Uxy^2));
    A1 = (-b+sqrt((b^2) - (4*a*c)))/(2*a);
    A2 = (-b-sqrt((b^2) - (4*a*c)))/(2*a);

    % Out-of-Plane Modes:
    oopm1 = 1i*sqrt(abs(Uzz));
    oopm2 = -oopm1;

    % System In/Out-Plane Modes (Eigen Values)
    out.eigVals = [sqrt(A1); -sqrt(A1); sqrt(A2); -sqrt(A2); oopm1; oopm2];

end