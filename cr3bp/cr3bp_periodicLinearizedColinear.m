function out = cr3bp_periodicLinearizedColinear(c3sys, xyz, e0, n0, opts)
%CR3BP_PERIODICCOLINEAROSCILLATORY Finds Linearized Periodic solution for 
%                                  colinear Lagrange points (L1, L2, L3)
%                                  Solutions entirely in-plane.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sys      [struct]    Structure Containing System L,V,T,mu vals.
%       2. xyz         [3x1]    Lagrange point [X;Y;Z] N.D. Position
%       3. e0          [1x1]    X-comp. Periodic Orbit Initial State
%       4. n0          [1x1]    Y-comp. Periodic Orbit Initial State
%       5. x0          [4x1]    Given full var. [xi0; dxi0; eta0; deta0]
%                               override velocity computation and use this.
%       6. optional
%               tf     [1x1]    Specify final time (default: 2*pi)
%               pltFig [T/F]    Plot Figure
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: Structure 'out' Containing Fields:
%     1&2. .eigV1/3    [1x1]    In-Plane Eigen Values 1 and 3
%     3&4. .alpha1/3   [1x1]    alpha value = (l^2 - Uxx)/(2l)
%     5&6. .Uxx/Uyy    [1x1]    Second partial pseudo-potential XX and YY
%     7&8. .A3/A4      [1x1]    A=inv(B)*x0 values 3 and 4 only
%       9. .xi         [nx1]    Xi  (X) values in time
%      10. .eta        [nx1]    Eta (Y) values in time
%      11. .t          [nx1]    Time values from 0-->tf
%      12. .x0         [4x1]    Initial State from Lagrange Point
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
        c3sys;
        xyz;
        e0;
        n0;
        opts.x0 {mustBeVector} = [NaN; NaN; NaN; NaN];
        opts.tf {mustBeNumeric} = NaN;
        opts.pltFig {mustBeNumericOrLogical} = false;
    end

    % Find Second Partials and Eigen Values at XYZ Point
    pp  = cr3bp_eigValsLinSys(c3sys, xyz);
    l1  = pp.eigVals(1);
    l3  = pp.eigVals(3);
    Uxx = pp.jUjU(1,1);
    Uyy = pp.jUjU(2,2);
    Uxy = pp.jUjU(1,2);
    a1  = ((l1^2)-Uxx)/(2*l1 + Uxy);
    a3  = ((l3^2)-Uxx)/(2*l3 + Uxy);

    % Find Linearized Initial State (Full)
    de0 = (l3*n0)/a3;
    dn0 = a3*l3*e0;
    
    x0  = opts.x0;
    if isnan(x0(1))
        x0_ = [e0; de0; n0; dn0]; 
    else
        x0_ = x0;
    end

    % Find A3 and A4 Values for Linear Prop.
    B    = [   1     1     1     1;
              l1   -l1    l3   -l3;
              a1   -a1    a3   -a3;
           a1*l1 a1*l1 a3*l3 a3*l3];
    A    = inv(B)*x0_;
    A3   = A(3);
    A4   = A(4);

    % Linear Propagation
    if isnan(opts.tf)
        tf = (2*pi)/(abs(sqrt(-1)*l3));
    else
        tf = opts.tf;
    end
    
    T = [0:0.005:tf];
    for i=1:length(T)
        t    = T(i);
        E(i,1) =    A3*exp(l3*t) +    A4*exp(-l3*t);
        N(i,1) = a3*A3*exp(l3*t) - a3*A4*exp(-l3*t);

        E(i,1) = real(E(i,1)) + xyz(1);
        N(i,1) = real(N(i,1)) + xyz(2);
    end

    % Outputs
    out.eigV1  = l1;        out.eigV3  = l3;
    out.alpha1 = a1;        out.alpha3 = a3;
    out.Uxx    = Uxx;       out.Uyy    = Uyy;
    out.A3     = A3;        out.A4     = A4;
    out.xi     = E;
    out.eta    = N;
    out.t      = T;
    out.x0     = x0_;

    % Plot Figure
    if opts.pltFig
        figure()
        ft = tiledlayout(1,1,'Tilespacing','compact','Padding','compact');
        title(ft, ...
        ['Colinear Lagrange Point Periodic Orbit - Linearized System'], ...
            'fontsize',14,'interpreter','latex');
        hold on
        scatter(xyz(1),xyz(2),'r','filled','d');
        plot(E,N,'b','linewidth',1.5);
        hold off
        grid on; box on; axis equal;
        set(gcf,'color','w');
        legend({'Linerization Point', 'Periodic Orbit'}, ...
            'fontsize',12,'interpreter','latex');
        xlabel('X [ND]','fontsize',12,'interpreter','Latex');
        ylabel('Y [ND]','fontsize',12,'interpreter','Latex');
    end

end
