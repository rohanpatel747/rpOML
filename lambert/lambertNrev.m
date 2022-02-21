function out = lambertNrev(sv1, sv2, dt, mu, type)
%LAMBERTNREV Compute N Revolution Lambert Transfer Given Traj. Type
%
%   Assumptions/Warnings:
%       1. Solution Tollerance Value:   [ 1 second ]
%   	2. Method breaks down for very short or long dt values.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sv1  [6x1]or[3x1] Initial State or Position Vector (km)
%       2. sv2  [6x1]or[3x1] Final   State or Position Vector (km)
%       3. dt   [1x1]        Transfer TOF (s)
%       4. mu   [1x1]        Central Body Gravitational Parameter (km3/s2)
%       5. type [1x1]        Trajectory Type (i.e. 3 or 4) for a 1 rev.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. vi     [3x1]      Departure Velocity Vector (km/s)
%       2. vf     [3x1]      Arrival   Velocity Vector (km/s)
%       3. psiV   [mx1]      Universal Variable (Psi) Value @ Ittr. (rad.)
%       4. psimin [1x1]      Psi Value at Branch Break (Min. btwn Types)
%       5. tofV   [mx1]      Computed TOF @ Ittr. (s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Lambert's Problem" Professor Kate Davis
%       "Lambert-Handout-3.pdf"
%       2. Bate, Muller, and White "Fundamentals of Astrodynamics"
%       3. Vallado "Fundamentals of Astrodynamics and Applications"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   TODO:
%       1. Iteration counter to prematurely terminate the solver if
%          solution cannot be found
%


    tolVal = 1;     % Solution Tollerance Value (seconds)


    % Find Traj. Type Psi Limits
    rev = lambert_getRevfromType(type);
    [psiLL,psiUU] = lambert_getPsiNRev(rev);
    psiVals = psiLL+0.001:0.05:psiUU;
    for i=1:length(psiVals)
       tofVals(i) =  lambert_getTOFfromPsi(psiVals(i), sv1, sv2, mu);
    end

    % Find min for Type break
    [~,idxMin] = min(tofVals);
    psiMin = psiVals(idxMin);

    % Set Initial Bounds of Psi Based on Type
    if rem(type,2) == 1
        psiL = psiMin;
        psiU = psiLL;
        tofVals = tofVals(1:idxMin);
        psiVals = psiVals(1:idxMin);
    else
        psiL = psiMin;
        psiU = psiUU;
        tofVals = tofVals(idxMin:end);
        psiVals = psiVals(idxMin:end);
    end
    
    % A really good guess for Psi:
    for i=1:length(tofVals)
       dtofV(i) = abs(tofVals(i)-dt); 
    end
    [~,idxMin2] = min(dtofV);
    psi = psiVals(idxMin2);
    tofV= tofVals(idxMin2);
    
    psiV = psi;
    psiMinCalc = psiMin;
    
    
    
    
    % 0-Rev Lambert Algorithm
    r0_ = sv1(1:3); r0 = norm(r0_);
    rf_ = sv2(1:3); rf = norm(rf_);
    
    nu0 = atan2(r0_(2), r0_(1));
    nuf = atan2(rf_(2), rf_(1));
    dnu = nuf-nu0;
    if dnu<0
       dnu = (2*pi) + dnu; 
    end

    DM=1;
    if dnu > pi; DM=-1; end
   
    cdnu = dot(r0_,rf_)/(r0*rf);

    A = DM*sqrt(r0*rf*(1+cdnu));

    if (A==0) || (dnu==0)
       disp('Trajectory Invalid and Cannot be Computed'); 
       v0 = [NaN; NaN; NaN];
       vf = v0;
    else

        % Psi will always be >1 so default c2/c3 conditions
        c2 = (1.0-cos(sqrt(psi)))/psi;
        c3 = (sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);

        y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));
        xi = sqrt(y/c2);
        dti = ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));

        N = 0.8;

        l=2;
        while abs(dti-dt) > tolVal

            y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));

            if (A>0.00) && (y<0.0)
                psi = N*(1/c3)*(1 - ((sqrt(c2)/A)*(r0+rf)));
                [c2,c3] = lambert_getc2c3(psi);
                y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));
            end

            xi = sqrt(y/c2);

            dti= ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));

            % bisection method assumes positive slope (odd Types have neg
            % slope) We dont change this logic, but we flip our bounds of
            % psiL and psiU (lines 82-83 to turn a neg. slope to positive)
            if dti <= dt
                psiL = psi;
            else
                psiU = psi;
            end


            psi = (psiU + psiL)/2;
            [c2,c3] = lambert_getc2c3(psi);
            
            
            psiV(l) = psi;
            tofV(l) = dti;
            l=l+1;
        end


        f = 1-(y/r0);
        dg= 1-(y/rf);
        g = A*sqrt(y/mu);

        v0 = (rf_  - f*r0_)/g;
        vf = (dg*rf_ - r0_)/g;

    end
     
    % Outputs
    out        = struct;
    out.vi     = v0;
    out.vf     = vf;
    out.psiV   = psiV;
    out.psimin = psiMinCalc;
    out.tofV   = tofV;
    
end

