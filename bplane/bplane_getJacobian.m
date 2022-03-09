function jBjV = bplane_getJacobian(mu, x_, dv, opts)
%BPLANEBRBTFROMRV Returns the B-Plane BR and BT Values Only Given and RV
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. mu     [1x1]      Central Body Gravitational Parameter (km3/s2)
%       2. x_     [6x1]      State Vector w.r.t. Flyby Body (km and km/s)
%                            [x y; z; vx; vy; vz]
%       3. dv     [1x1]      Perturbation Magnitude (km/s)
%   [OPT]. npS    [1x1]      Number of States To Perturbate (default: 2)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output:
%       1. jBjV   [npSxnpS]  [nxn] Jacobian Matrix of Partials:
%                            B-Plane States to DV Direction
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    arguments
       mu
       x_;
       dv;
       opts.npS {mustBeInteger} = 2;
    end

    % Nominal B-Plane Values
    [BRn, BTn] = bplaneBRBTfromRV(mu,x_);
    TCAn       = bplaneTCA(mu,x_);

    % X Perturbation Vector Matrix
    dvM = [dv,  0,  0;
            0, dv,  0;
            0,  0, dv];

    % Construct Full BR,BT,TCA Jacobian Matrix (jBjBV) [3x3]
    jBjV = zeros(3,3);
    for i=1:3

        % Perturbation Vector
        dv_ = dvM(i,:).';

        % Perturbed Velocity Vector
        Vp_ = x_(4:6) + dv_;

        % Perturbed B-Plane Parameters
        [BRp,BTp] = bplaneBRBTfromRV( mu, [x_(1:3);Vp_]);
        TCAp      = bplaneTCA(        mu, [x_(1:3);Vp_]);

        % Update Jacobian Matrix
        jBjV(1,i) = (BTp - BTn)/dv;
        jBjV(2,i) = (BRp - BRn)/dv;
        jBjV(3,i) = (TCAp-TCAn)/dv;
        
    end
    
    if opts.npS == 2
        jBjV = jBjV(1:2,1:2);
    end
    
end

















