function J = cr3bp_computeJacobiConstant(c3sys, x_)
%CR3BP_COMPUTEJACOBICONSTANT Computes the Jacobi Constant at a State (x_)
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure Containing System L,V,T,mu vals.
%       2. x_         [6x1]      N.D. Rot. State (x;y;z;vx;vy;vz)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. J          [1x1]      Jacobi Constant Value
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    mu = c3sys.mu;

    x1 = c3sys.b1(1);   x2 = c3sys.b2(1);
    y1 = c3sys.b1(2);   y2 = c3sys.b2(2);
    z1 = c3sys.b1(3);   z2 = c3sys.b2(3);

    x  = x_(1);
    y  = x_(2);
    z  = x_(3);
    dx = x_(4);
    dy = x_(5);
    dz = x_(6);

    r1 = sqrt((x-x1)^2 + (y-y1)^2 + (z-z1)^2);
    r2 = sqrt((x-x2)^2 + (y-y2)^2 + (z-z2)^2);

    J  = (x^2 + y^2) + ((2*(1-mu))/r1) + ((2*mu)/r2) - dx^2 - dy^2 - dz^2;

end

