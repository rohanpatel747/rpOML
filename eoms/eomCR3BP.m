function dYdt = eomCR3BP(t, Y, mu)
%EOMCR3BP CR3BP EOMs [6x1] 
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1 Int.]  dt time step from ODE45
%       2. Y          [6x1 Dbl.]  State to pass into integ. instant
%       3. mu         [1x1 Dbl.]  Secondary Body Gravitational Param.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [6x1 Dbl.]  State Space Vector for Integration
%

x = Y(1);
y = Y(2);
z = Y(3);

dx= Y(4);
dy= Y(5);
dz= Y(6);

r1 = sqrt(  ((x+mu))^2 + (y)^2 + (z)^2);
r2 = sqrt((x-(1-mu))^2 + (y)^2 + (z)^2);

d2x = x + 2*dy -(1-mu)*((x+mu)/(r1^3)) - mu*((x-(1-mu))/(r2^3));
d2y = y - 2*dx -(1-mu)*((y)/(r1^3))    - mu*(y/(r2^3));
d2z =          -(1-mu)*((z)/(r1^3))    - mu*(z/(r2^3));

dYdt = [dx; dy; dz; d2x; d2y; d2z];

end

