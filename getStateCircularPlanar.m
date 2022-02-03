function out = getStateCircularPlanar(dt, xi, mu)
% GETSTATECIRCULAR Returns Cartesian State w.r.t Time (CIRCULAR-PLANAR ORBITS ONLY)
%
%   Assumptions/Warnings:
%   	1. Circular Planar Orbits ONLY
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. dt [1x1] Time since given initial state (xi)
%       2. xi [6x1] Initial State [x y z vx vy vz]
%       3. mu [1x1] Central Body Gravitational Parameter
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Outputs:
%       1. out [6x1] Inertial State [x y z vx vy vz]
%

    % Geometry
    x_ = [1; 0; 0];
    r_ = xi(1:3);
    v_ = xi(4:6);
    r  = norm(r_);
    v  = norm(v_);

    % Initial Angle
    ai= atan2(norm(cross(r_,x_)),dot(r_,x_));
    if r_(2)<0
        ai = (2*pi)-ai;
    end

    % Orbit Period
    T = 2*pi*sqrt((r^3)/mu);

    % Angle Traveled
    dA = (2*pi*dt)/T;

    % New Angle
    a = ai+dA;

    % New State
    rnew_ = [r*cos(a); r*sin(a); 0.00];
    vnew_ = [v*cos(a); v*sin(a); 0.00];
    xnew_ = [rnew_;vnew_];
    out = xnew_;

 end
 
 