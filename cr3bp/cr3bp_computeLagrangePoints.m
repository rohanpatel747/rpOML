function c3sys = cr3bp_computeLagrangePoints(c3sys, opts)
%CR3BP_COMPUTELAGRANGEPOINTS Finds Lagrange points given system.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure containing at least:  c3sys.mu
%                                where mu is the characteristic mass ratio
%
%       2. 'tolVal'   [1x1]      [Optional] Specified tolerance value.
%                                           Default: 1e-10
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (figure)
%       1. c3sys       [struct]  Input structure but with additional fields
%           .Ln           [3x1]  [x;y;z] ND Rot. Location of Ln for n=1:5
%           .jUjL123check [3x1]  jU/jx evaluated at the L1,2,3 x-coords. to
%                                verify that the point has a y=0 value.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
        c3sys;
        opts.tolVal = 1e-10;
    end

    mu = c3sys.mu;
    tolVal = opts.tolVal;

    
    % Compute L1, L2, and L3
    xi1 = (1-mu) - 0.2;
    xi2 = (1-mu) + 0.2;
    xi3 = ( -mu) - 0.2;
    [x1, fx1] = jUjxNewtonRoot(xi1, xi1+0.005, mu, tolVal);
    [x2, fx2] = jUjxNewtonRoot(xi2, xi2+0.005, mu, tolVal);
    [x3, fx3] = jUjxNewtonRoot(xi3, xi3+0.005, mu, tolVal);


    % Lagrange Points
    c3sys.L1 = [      x1;          0; 0];
    c3sys.L2 = [      x2;          0; 0];
    c3sys.L3 = [      x3;          0; 0];
    c3sys.L4 = [(1/2)-mu;  sqrt(3)/2; 0];
    c3sys.L5 = [(1/2)-mu; -sqrt(3)/2; 0];


    % Check if L1, L2, and L3, points' have y=~0 
    c3sys.jUjL123check = [fx1; fx2; fx3];


    % Newton Root Finding Function for JU/Jx
    function [x, fx] = jUjxNewtonRoot(x,xi,mu,tolVal)
    
        while abs(x-xi) > tolVal
            [fx, dfx] = jUjx(x,mu);
            xi= x;
            x = x - fx/dfx;
        end
    
        % Function computing JU/Jx and d(JU/Jx) for Newton Root Finding
        function [fx, dfx] = jUjx(x,mu)
            fx  = x - (((1-mu)*(x+mu))/(abs(x+mu)^3)) - ...
                ((mu*(x-1+mu))/(abs(x-1+mu)^3));
            dfx = 1 + ((2*(1-mu))/(((x+mu)^2)*abs(x+mu))) + ...
                ((2*mu)/(((x+mu-1)^2)*abs(x+mu-1)));
        end
    end

end