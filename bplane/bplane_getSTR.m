function out = bplane_getSTR(vi_)
%BPLANE_GETSTR Gets the STR Coordinate Frame Represented in Inertial Coords
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. vi_    [3x1]     Incoming Vinfinity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. S_     [3x1]     S_hat Unit Vector
%       2. T_     [3x1]     T_hat Unit Vector
%       3. R_     [3x1]     R_hat Unit Vector
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    out = struct;
    out.S_ = vi_/norm(vi_);
    out.T_ = cross(out.S_,[0;0;1])/norm(cross(out.S_,[0;0;1]));
    out.R_ = cross(out.S_,out.T_)/norm(cross(out.S_,out.T_));

end

