function dYdt = eomCR3BPwSTM(t, Y, mu)
%EOMCR3BPwSTM CR3BP EOMs [1:6] with STM Propagation [6x6]-->[7:42]
%             Integrate with ODE113 for best performance.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1  Int.]  dt time step from ODE113
%       2. Y          [42x1 Dbl.]  State to pass into integ. instant
%       3. mu         [1x1  Dbl.]  Secondary Body Gravitational Param.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [42x1 Dbl.] State Space and STM from Integration
%
x = Y(1);
y = Y(2);
z = Y(3);

dx= Y(4);
dy= Y(5);
dz= Y(6);

r1 = (((x+mu)^2)     + (y^2) + (z^2))^(3/2);
r2 = (((x-(1-mu))^2) + (y^2) + (z^2))^(3/2);

d2x = x + 2*dy -(1-mu)*((x+mu)/r1) - mu*((x-(1-mu))/r2);
d2y = y - 2*dx -(1-mu)*((y)/r1) - mu*(y/r2);
d2z = -(1-mu)*((z)/r1) - mu*(z/r2);

dYdt_ = [dx; dy; dz; d2x; d2y; d2z];

stm = reshape(Y(7:end), 6,6);
dstm= cr3bp_dstm([x;y;z],mu)*stm;

dstm_= reshape(dstm, 36,1);

dYdt = [dYdt_; dstm_];

end