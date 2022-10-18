function out = cr3bp_familyPSAL(c3sys,V,ds,nsp)
%CR3BP_FAMILYPSAL Construct Periodic Orbits with Pseudo-Arc Lenght Method 
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct  ]  Structure Containing at least Sys. mu val
%       2. V          [7x1 Dbl.]  Initial Periodic Orbit (corrected)
%       3. ds         [1x1 Dbl.]  Step along Tangent (Pseudo-arc len. mag.)
%       4. ns         [7x1 Dbl.]  Previous Null Space Vector
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. V          [7x1 Dbl.]  [Corrected State ; Corrected Prop.Time]
%       4. fval       [itrx1   ]  Array of constraint norms at each itr.
%       5. itr        [1x1 Dbl.]  Number of Iterations Req. to Compute
%       6. jc         [1x1 Dbl.]  Jacobi Constant of Corrected State
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Notes/References:
%       1. "Computing Periodic Orbit Families, Part 1 & 2" 
%          ASEN6060 Advanced Astrodynamics
%          |--> Pseudo Arc-Length (PSAL) Method
%

    mu = c3sys.mu;


    % Initial Corrected Periodic Orbit
    poi = cr3bp_correctPO_SS_GVTMCF(c3sys, V);  % should be 1 itr. max.
    Vd  = poi.V;
    ns  = null(poi.df);
    ns  = ns(:,1);

    if ns(1)-nsp(1) < 0
        ns = -ns;
    end


    
    % Step for State
    V = Vd + ds*ns;


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
        dh = DH(stm,xs,ns,mu);
        h  = H([xs;tp],xp,Vd,ns,ds);
        V  = V - dh\h;

        % Loop Termination Condition
        fval(itr+1) = norm(h);
        itr=itr+1;
        %if itr>20; break; end
    
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Outputs
    out.V    = V;
    out.fval = fval;
    out.itr  = itr;
    out.jc   = cr3bp_computeJacobiConstant(c3sys, V(1:6));
    out.ns   = ns;


    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Constraint Vector (H)
    function h  = H(V,xi,Vd,ns,ds)
        h  = [V(1)-xi(1);
              V(2)-xi(2);
              V(3)-xi(3);
              V(4)-xi(4);
              V(6)-xi(6);
                   xi(2);
              (V - Vd).'*ns - ds];
    end
    
    % Jacobian of State/Constraints (DH)
    function dh = DH(stm,xs,ns,mu)
        eom = eomCR3BP(0,xs,mu);
        stm = stm - eye(6);
        stm = [stm(1:4,1:6); stm(6,1:6); [0,1,0,0,0,0]];
        dh  = [stm, [eom(1:4); eom(6); 0]; ns.'];
    end

end