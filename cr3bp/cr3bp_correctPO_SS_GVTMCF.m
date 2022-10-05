function out = cr3bp_correctPO_SS_GVTMCF(c3sys, V)
%CR3BP_CORRECTPO_SS_GVTMCF Correct Periodic Orbit with Single Shooting (SS)
%                          GVT : General Variable Time
%                          MCF : MODIFIED Constraint Formulation
%   Modified Formulation
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure Containing System L,V,T,mu vals.
%       2. V          [7x1]      Initial Free Vector V=[x;y;z;vx;vy;vz;T]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: structure 'out' with fields:
%       1. .V         [7x1]      Corrected Free Vector V=[x;y;z;vx;vy;vz;T]
%       2. .itr       [1x1]      Number of iterations used
%       3. .fval      [itrx1]    Constraint Norm. at each itr.
%       4. .df        [6x7]      Jacobian of State/Constraint at final V
%       5. .stm       [6x6]      State Transition Matrix at final V
%       6. .jc        [1x1]      Jacobi Constant of final V
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6060 Lecture 6 'Computing Periodic Orbit Families, Part I'
%

    mu = c3sys.mu;
    
    % Single Shooting Loop
    options = odeset('reltol',1e-13,'abstol',1e-13);
    tolval  = 1e-10;
    fval    = 1;
    itr     = 0;
    while fval>tolval
        
        % State and Time Integration
        tp    = V(end);
        xp    = [V(1:end-1); reshape(eye(6),36,1)];
        [~,x] = ode113(@(t,Y) eomCR3BPwSTM(t,Y,mu), [0 tp], xp, options);
        xs    = x(end,1:6).';
        stm   = reshape(x(end,7:end).',6,6);
        
        % Update Equation
        df = DF(stm,xs,mu);
        f  = F(xs, xp(1:6));
        V  = V - df.'*((df*df.')\f);
    
        % Loop Termination Condition
        fval(itr+1) = norm(f);
        itr=itr+1;
    
    end
    
    % Outputs
    out.V   = V;
    out.itr = itr;
    out.fval= fval;
    out.df  = df;
    out.stm = stm;
    out.jc  = cr3bp_computeJacobiConstant(c3sys,V(1:6));

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Constraint and Jacobian of State/Constraint Functions
    function f  = F(xf,xi)
        f   = [xf(1)-xi(1);
               xf(2)-xi(2);
               xf(3)-xi(3);
               xf(4)-xi(4);
               xf(6)-xi(6);
               xi(2)];
    end

    function df = DF(stm,xs,mu)
        eom  = eomCR3BP(0,xs,mu);
        stm = stm - eye(6);
        stm = [stm(1:4,1:6); stm(6,1:6); [0,1,0,0,0,0]];
        df  = [stm, [eom(1:4); eom(6); 0]];
    end

end