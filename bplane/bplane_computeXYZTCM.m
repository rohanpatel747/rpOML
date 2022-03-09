function out = bplane_computeXYZTCM(mu,x_,BTarg,opts)
%BPLANE_COMPUTEXYZTCM Returns DV and state for TCM given desired BP Targ.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]      Central Body Gravitational Parameter (km3/s2)
%       2. x_     [6x1]      State Vector w.r.t. Flyby Body (km and km/s)
%                            [x y; z; vx; vy; vz]
%       3. BTarg  [3x1]      Target B-Plane Vector Including Time to C.A.:
%                            [B*R; B*T; TCA]
%   [OPT]. dv     [1x1]      Perturbation Mag. (km/s) (default: 0.001)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: 'out' structure with fields:
%       1. DV_    [3x1]      Required Delta-V Vector (km/s)
%       2. DV     [1x1]      Norm. of Delta-V Vector (km/s)
%       3. x__    [6x1]      New Full State Vector w.r.t. Flyby Body
%                            [x y; z; vx; vy; vz] (km) and (km/s)
%       4. itt    [1x1]      Number of Iterations to Converge
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. bplaneBRBTfromRV()
%       2. bplane_getJacobian()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
       mu
       x_;
       BTarg;
       opts.dv {mustBeNumeric} = 0.001;
    end

    xi_= x_;                % Initial State
    BRt = BTarg(1);         % Target B*R
    BTt = BTarg(2);         % Target B*T
    TCAt= BTarg(3);         % Target Time to Close Approach
    
    dv = opts.dv;           % Perturbation Magnitude
    i=0;                    % Iteration Count
    tolValBRBT = 1e-6;      % Tolerance Value of BRBT
    tolValTCA  = 1;         % Tolerance Value of TCA   
    DB     = [1e6;1e6;1e6];
    while (abs(DB(1))>tolValBRBT) && (abs(DB(2))>tolValBRBT) %&& (abs(DB(3))>tolValTCA)

        % Nominal B-Plane Values
        [BRn, BTn] = bplaneBRBTfromRV(mu, x_);
        TCAn       = bplaneTCA(mu,x_);

        % Compute Nominal Partials Matrix
        jBjV = bplane_getJacobian(mu, x_, dv, 'npS', 3);

        % BR BT Diff.
        DB = [BTt-BTn; BRt-BRn; TCAt-TCAn];
        DVV= jBjV\DB;

        % Compute New State
        V_ = x_(4:6) + DVV;
        x_ = [x_(1:3); V_];

        % Number of Iterations
        i=i+1;

    end

    % Find DV Vector for TCM
    Dx_ = x_ - xi_;
    DV_ = Dx_(4:6);
    DV  = norm(DV_);

    % Outputs
    out     = struct;
    out.DV_ = DV_;
    out.DV  = DV;
    out.xf_ = x_;
    out.itt = i;

end
