function dYdt = eomCR3BPwSTM6060(t, Y, mu)
%EOMCR3BPwSTM CR3BP EOMs [1:6] with STM Propagation [6x6]-->[7:42]
%             Integrate with ODE113 for best performance.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1  Int.]  dt time step from ODE113
%       2. Y          [42x1 Dbl.]  State to pass into integ. instant
%                                  [x;y;z;vx;vy;vz;STM(1:36,1)]
%       3. mu         [1x1  Dbl.]  Secondary Body Gravitational Param.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [42x1 Dbl.] State Space and STM from Integration
%

% State Unpacking
x = Y(1);
y = Y(2);
z = Y(3);

dx= Y(4);
dy= Y(5);
dz= Y(6);

r1 = sqrt(  ((x+mu))^2 + (y)^2 + (z)^2);
r2 = sqrt((x-(1-mu))^2 + (y)^2 + (z)^2);

stm = reshape(Y(7:end), 6,6);

% State Derivaive Computation 
d2x = x + 2*dy -(1-mu)*((x+mu)/(r1^3)) - mu*((x-(1-mu))/(r2^3));
d2y = y - 2*dx -(1-mu)*((y)/(r1^3)) - mu*(y/(r2^3));
d2z = -(1-mu)*((z)/(r1^3)) - mu*(z/(r2^3));

dYdt_ = [dx; dy; dz; d2x; d2y; d2z];

% STM Derivative Computation 
Uxx = 1 - ((1-mu)/(r1^3)) - (mu/(r2^3)) + ...
    ((3*(1-mu)*((x+mu)^2))/(r1^5)) + ((3*(mu)*((x-1+mu)^2))/(r2^5));

Uxy = ((3*(1-mu)*(x+mu)*y)/(r1^5)) + ((3*(mu)*(x-1+mu)*y)/(r2^5));
Uyx = Uxy;

Uxz = ((3*(1-mu)*(x+mu)*z)/(r1^5)) + ((3*(mu)*(x-1+mu)*z)/(r2^5));
Uzx = Uxz;

Uyy = 1 + ((3*(1-mu)*(y^2))/(r1^5)) - ((1-mu)/(r1^3)) - ((mu)/(r2^3)) ...
    + ((3*mu*(y^2))/(r2^5));

Uyz = ((3*(1-mu)*y*z)/(r1^5)) + ((3*mu*y*z)/(r2^5));
Uzy = Uyz;

Uzz = ((3*(1-mu)*(z^2))/(r1^5)) - ((1-mu)/(r1^3)) + ((3*mu*(z^2))/(r2^5)) ...
    - (mu/(r2^3));

jUjU = [Uxx Uxy Uxz; Uyx Uyy Uyz; Uzx Uzy Uzz];

At = [zeros(3), eye(3); jUjU, [0 2 0; -2 0 0; 0 0 0]];

dstm = At*stm;
    
% Repack State 
dstm_= reshape(dstm, 36,1);
dYdt = [dYdt_; dstm_];

end