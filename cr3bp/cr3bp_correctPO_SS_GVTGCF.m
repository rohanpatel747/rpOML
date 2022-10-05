function out = cr3bp_correctPO_SS_GVTGCF(c3sys, V)
%CR3BP_CORRECTPO_SS_GVTGCF Correct Periodic Orbit with Single Shooting (SS)
%                          GVT : General Variable Time
%                          GCF : General Constraint Formulation
% 
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
    
    F       = @(xf,xi) xf-xi;
    DF      = @(stm,xs,mu) [(stm-eye(6)), eomCR3BP(0,xs,mu)]; 
    
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

end