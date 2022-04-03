function out = cr3bp_convertRotNDToInertial(sys, time, state, n)
%CR3BP_CONVERTTOINERTIAL Converts State and Time Arrays to Inertial Coords.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sys      [struct]    Structure Containing System L,V,T,mu vals.
%       2. time     [nx1]       Array of Non-Dim. Time Values
%       3. state    [nx6]       Array of Non-Dim. States
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: Array 'out' Containing Columns:
%       1. out(:,1) - Dim. Time (s)
%       2. out(:.2:6) - [x, y, z, vx, vy, vz]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Description:
%         Conversion Rules from Non-Dim. to Dim. System:
%             Distance: di = L*d
%             Velocity: si = V*s
%             Time: ti = t*(T/(2*pi))
%             (X,Y,Z) = Inertial
%             (x,y,z) = Rotational
%             Position: ( X, Y, Z)' = A * (x_, y_, z_)'
%             Position: ( x, y, z)' = L * (x_, y_, z_)'  
%             Velocity: (dX_, dY_, dZ_)' = A*(dx-y, dy+x, dz)'
%             Velocity: (dX,dY,dZ)' = V * (dx_, dy_, dz_)'
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. Equations of Motion (Lo Section 2.3, p.28_bk)
%

    L = sys.L;
    V = sys.V;
    T = sys.T;
   

    for i=1:height(state)
        p = state(i,1:3).';
        v = state(i,4:6).';
        t = time(i);

        A = [cos(t), -sin(t),   0;
             sin(t),  cos(t),   0;
                  0,       0,   1];

        p_i = A*p;
        v_i_= [v(1)-p(2);  v(2)+p(1);  v(3)];
        v_i = A*v_i_;
        

        

        out(i,1) = t * (T/(2*pi));
        out(i,2) = L * p_i(1);
        out(i,3) = L * p_i(2);
        out(i,4) = L * p_i(3);
        out(i,5) = V * v_i(1);
        out(i,6) = V * v_i(2);
        out(i,7) = V * v_i(3);

    end

end