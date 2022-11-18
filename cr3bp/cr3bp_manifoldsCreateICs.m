function out = cr3bp_manifoldsCreateICs(c3sys, npts, xpoi, tpofull, opts)
%CR3BP_MANIFOLDSCREATEICS Creates Initial Conditions along PO for Manifolds
%   This function creates the states and associated monodromy matrices for
%   PO manifold creation. Also, it creates the specific Plus/Minus ICs for
%   stable, unstable, center, and spiral in/out manifolds.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct]   Structure containing system mu value
%       2. npts       [1x1]      Number of states to compute along the PO
%       3. xpoi       [6x1]      N.D. Rot Initial State on PO
%       4. tpofull    [1x1]      N.D. Orbit Period of the PO.
%   Optional Inputs:
%       5. 'computeManifoldIC' [T/F] 
%       6. 'dd'       [1x1]      Scaling Term for Eigen-Vector
%                                [Default: 0.75e-4]
%       7. 'roundingValue [1x1]  Round to This Decimal Place
%                                [Default:5]
%       8. 'idx'      [T/F]      Specify Indices for Stable/US Modes
%                                [Default: False]
%       9. 'usidx'    [1x1]      If idx=true Unstable Index [Default: 1]
%      10. 'sidx'     [1x1]      If idx=true Stable   Index [Default: 2]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: structure 'out' with fields: 
%       1. .xFpo      [nx42]     'n' is the number of points (npts).
%                                xFpo(i,1:6) = State Along PO
%                                xFpo(i,7:42)= Associated Monodromy Matrix
%       2. .xFpoRef   [nnx42]    Given PO's Integrated state for 1 Period
%       3. .xUSp      [nx6]      Unstable Manifold (+) Initial Conditions
%       4. .xUSm      [nx6]      Unstable Manifold (-) Initial Conditions
%       5. .xSp       [nx6]      Stable   Manifold (+) Initial Conditions
%       6. .xSm       [nx6]      Stable   Manifold (-) Initial Conditions
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
        c3sys;
        npts;
        xpoi;
        tpofull;
        opts.computeManifoldIC {mustBeNumericOrLogical} = true;
        opts.dd                {mustBeNumeric}          = 0.75e-4;
        opts.roundingValue     {mustBeNumeric}          = 5;
        opts.idx   = false;
        opts.usidx = 1;
        opts.sidx  = 2;
    end


    % Compute States and STMs along PO
    mu      = c3sys.mu;
    options = odeset('reltol',1e-13,'abstol',1e-13);

    [~,xFpo]    = ode113(@(t,Y) eomCR3BPwSTM6060(t,Y,mu), ...
        [0:(tpofull/npts):tpofull], [xpoi;reshape(eye(6),36,1)], options);
    monodromy   = reshape(xFpo(end,7:42),6,6);

    
    % Find Initial Eigen-Vector and Propagate it using STM
    [eV_,~] = eig(monodromy);
    for i=2:npts
        stm = reshape(xFpo(i,7:42),6,6);
        for j=1:6
            eV_(:,j,i) = stm*eV_(:,j,1);
        end
    end

    
    % Check for Unstable(US)/Stable(S) Manifolds
    rd = opts.roundingValue;
    d = eig(monodromy);
    if opts.idx == false
        for i=1:length(d)
            d_ = d(i);
            d_r= abs(round(real(d_),rd));
            d_i= abs(round(imag(d_),rd));
    
            if     (d_r>1) && (d_i==0)  % US
                usidx = i;
            elseif (d_r<1) && (d_i==0)  % S
                sidx  = i;
            end
        end
    else
        usidx = opts.usidx;
        sidx  = opts.sidx;
    end


    % Compute Initial Conditions of Manifolds
    dd = opts.dd;
    for i=1:npts
        
        xx0 = xFpo(i,1:6).';
        v = eV_(:,:,i);

        % Unstable Manifold State and Integ.
        us_ = v(:,usidx);
        us_ = us_./norm(us_(1:3));
        if us_(1)<0; us_ = -us_; end
        xUSp(i,1:6) = (xx0 + dd.*us_).';
        xUSm(i,1:6) = (xx0 - dd.*us_).';
        
        % Stable Manifold State and Integ.
        s_ = v(:,sidx);
        s_ = s_./norm(s_(1:3));
        if s_(1)<0; s_ = -s_; end
        xSp(i,1:6) = (xx0 + dd.*s_).';
        xSm(i,1:6) = (xx0 - dd.*s_).';

        % Force Real Solutions Only
        xUSp(i,1:6) = real(xUSp(i,1:6));
        xUSm(i,1:6) = real(xUSm(i,1:6));
        xSp(i,1:6)  = real(xSp(i,1:6));
        xSm(i,1:6)  = real(xSm(i,1:6));

    end


    % Outputs
    out.d       = d;
    out.usidx   = usidx;
    out.sidx    = sidx;
    out.eV_     = eV_;

    out.xFpo    = xFpo;

    out.xUSp    = xUSp;
    out.xUSm    = xUSm;
    out.xSp     = xSp;
    out.xSm     = xSm;

end



