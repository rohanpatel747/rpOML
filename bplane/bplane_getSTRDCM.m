function in = bplane_getSTRDCM(in)
%BPLANE_GETSTRDCM Gets the STR->Inertial Direction Cosine Matrix
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%   	1. in     [struct]  Input Structure Containing Fields:
%          in.S_  [3x1]     S Unit Vector (km)
%          in.T_  [3x1]     T Unit Vector (km)
%          in.R_  [3x1]     R Unit Vector (km)
%
%       < Input 1 can be solved using:    bplane_getSTR(vinf_in) >
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: in structure containing fields:
%       1. S_     [3x1]     S_hat Unit Vector
%       2. T_     [3x1]     T_hat Unit Vector
%       3. R_     [3x1]     R_hat Unit Vector
%       4. STRDCM [3x3]     Direction Cosine Matrix of STR to Inertial
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Example:
%       vi1_str = dcm*vi1_; = [vinfmag; 0; 0]
%

    S = in.S_;  X = [1;0;0];
    T = in.T_;  Y = [0;1;0];
    R = in.R_;  Z = [0;0;1];

    in.DCM = [dot(S,X) dot(S,Y) dot(S,Z);  
              dot(T,X) dot(T,Y) dot(T,Z);
              dot(R,X) dot(R,Y) dot(R,Z)];

end

