function rev = lambert_getRevfromType(type)
%LAMBERT0REV Compute 0 Revolution Lambert Transfer
%
%   Assumptions/Warnings:
%   	1. MAX REV = 4 (Type 9 and 10 Trajectories)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. type [1x1] Trajectory Type (i.e. III (3), or IV (4))
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. rev [1x1] Associated Orbit Revolutions
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    if     (type==3) || (type==4)
        rev = 1;
    elseif (type==5) || (type==6)
        rev = 2;
    elseif (type==7) || (type==8)
        rev = 3;
    elseif (type==9) || (type==10)
        rev = 4;
    end
end

