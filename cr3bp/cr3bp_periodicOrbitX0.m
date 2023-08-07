function out = cr3bp_periodicOrbitX0(X0_, correct_z0, mu)
%CR3BP_PERIODICORBITX0 Adjusts Initial State to Achieve a Periodic Orbit
%
%   Assumptions/Warnings:
%   	1. Method uses Single Shooting (linearization).
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. X0_       [42x1]  Initial State Vector and STM Identity Reshaped
%       2. correct_z0 [T/F]  Apply Corrections to z0 and dy0:
%                               True : corrects initial z and dy0
%                               False: corrects initial x and dy0 
%       3. mu         [1x1]  CR3BP System Gravitational Parameter
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. itr        [1x1]  Number of Iterations Required
%       2. X0_       [42x1]  Initial State Vector (w/ STM)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Single Shooting Handout Interplanetary Mission Design"
%          Professor Kate Davis "SingleShootingHandout.pdf"
%       2. ASEN6008 "Flowchart for constructing simple, Periodic Orbits"
%          Professor Kate Davis "SingleShootingFlowchart.pdf"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % Conditions Setup
    maxittr= 100;
    tolVal = 1e-10;
    options = odeset('Events', @cr3bp_yaxstopcond);

    % Corrections Loop
    ittr = 0;
    uwm(1:2,1) = 1;
    while (norm(uwm(1))>tolVal) && (norm(uwm(2))>tolVal)
        if ittr<maxittr
            
            % Integrate w/ ODE113
            [~,  x,~,~,~] = ode113(@(t,Y) eomCR3BPwSTM(t,Y,mu), [0 8*pi], X0_,options);

            % State and STM at 1st Y-Axis Crossing
            YT2_ = x(end,:);
            stm  = cr3bp_getSTMfromY(YT2_, 1);  % STM(T/2, 0)

            % Get Accels at Time T/2
            dYdT2 = eomCR3BP(0, YT2_, mu);

            % Unwanted Motion
            ddx = -x(end,4);
            ddz = -x(end,6);
            uwm = [ddx; ddz];

            % Correct Unwanted Motion
            if correct_z0
                % Hold x0 fixed, vary z0 and dy0
                Y0dy0_1 = [stm(4,3) stm(4,5);
                           stm(6,3) stm(6,5)];
                Y0dy0_2 = [stm(2,3) stm(2,5)];
            else
                % Hold z0 fixed, vary x0 and dy0
                Y0dy0_1 = [stm(4,1) stm(4,5);
                           stm(6,1) stm(6,5)];
                Y0dy0_2 = [stm(2,1) stm(2,5)];
            end
            uwmc = Y0dy0_1 - (1/dYdT2(2))*[dYdT2(4); dYdT2(6)]*Y0dy0_2;

            % Compute Correction
            uwmCor = inv(uwmc)*uwm;
            if correct_z0
                X0_ = [X0_(1);           X0_(2); X0_(3)+uwmCor(1);
                       X0_(4); X0_(5)+uwmCor(2);           X0_(6)];
            else
                X0_ = [X0_(1)+uwmCor(1);           X0_(2); X0_(3);
                                 X0_(4); X0_(5)+uwmCor(2); X0_(6)];
            end

            % Update I.C.
            X0_(7:42) = reshape(eye(6), 36,1);
            ittr = ittr+1
            
        else
            
            disp(['MAX ITERATIONS (', num2str(maxittr), ...
                ') REACHED. UNABLE TO SOLVE X0_']);
            disp(' ');
            break;
            
        end

    end

    % Outputs
    out = struct;
    out.itr = ittr;
    out.X0_ = X0_;

    
    % Y-Axis Crossing Termination Function
    function [value,isterminal,direction] = stopCond(~,x)
        value = x(2);
        isterminal = 1;
        direction = 0;   
    end

end

