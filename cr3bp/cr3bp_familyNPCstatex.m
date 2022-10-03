function out = cr3bp_familyNPCstatex(c3sys,Vi,dp)
%CR3BP_FAMILYNPCSTATEX Construct Periodic Orbits with y0=0 constraint and
%                      varying the x0 component by (dp).
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct  ]  Structure Containing at least Sys. mu val
%       2. Vi         [7x1 Dbl.]  [Initial State ; Initial Prop.Time]
%       3. dp         [1x1 Dbl.]  Spacing btwn. computed PO family members
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. V          [7x1 Dbl.]  [Corrected State ; Corrected Prop.Time]
%       2. Vi         [7x1 Dbl.]  [Initial State ; Initial Prop.Time]
%       3. Vip        [7x1 Dbl.]  Adjusted (by dp) initial state & time
%       4. fval       [itrx1   ]  Array of constraint norms at each itr.
%       5. itr        [1x1 Dbl.]  Number of Iterations Req. to Compute
%       6. jc         [1x1 Dbl.]  Jacobi Constant of Corrected State
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Notes/References:
%       1. "Computing Periodic Orbit Families, Part 1 & 2" 
%          ASEN6060 Advanced Astrodynamics
%          |--> Natural Parameter Continuation (NPC) Method
%

    mu  = c3sys.mu;
    x0s = Vi(1:6);
    t0  = Vi(end);
    
    ps = x0s(1);
    x0 = [x0s(1)+dp; x0s(2:end)];
    V  = [x0; t0];

    % Single Shooting Loop
    options = odeset('reltol',1e-13,'abstol',1e-13);
    tolval  = 1e-10;
    fval    = 1;
    itr     = 0;
    while fval>tolval

        % State and Time Integration
        tp    = V(end);
        xp    = [V(1:end-1); reshape(eye(6),36,1)];
        [~,x] = ode113(@(t,Y) eomCR3BPwSTM6060(t,Y,mu), [0 tp], xp, options);
        xs    = x(end,1:6).';
        stm   = reshape(x(end,7:end).',6,6);

        % Update Equation
        df = DF(stm,xs,mu);
        f  = F(xs, xp(1:6), ps, dp);
        V  = V - df\f;

        % Loop Termination Condition
        fval(itr+1) = norm(f);
        itr=itr+1;

    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Outputs
    out.V    = V;
    out.Vi   = Vi;
    out.Vip  = [x0; t0];
    out.fval = fval;
    out.itr  = itr;
    out.jc   = cr3bp_computeJacobiConstant(c3sys, V(1:6));

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    % Constraint Function (state-x var. condition)
    function f  = F(xf,xi, ps, dp)
        f   = [xf(1)-xi(1);
               xf(2)-xi(2);
               xf(3)-xi(3);
               xf(4)-xi(4);
               xf(6)-xi(6);
               xi(2);
               xi(1) - (ps+dp)];
    end

    % Jacobian of State/Constraints
    function df = DF(stm,xs,mu)
        eom  = eomCR3BP(0,xs,mu);
        %df   = [stm-eye(6), eom; stm(1,:), eom(1)];
        stm = stm - eye(6);
        stm = [stm(1:4,1:6); stm(6,1:6); [0,1,0,0,0,0]];
        df  = [stm, [eom(1:4); eom(6); 0]; [1,0,0,0,0,0,0]];
    end

 end
