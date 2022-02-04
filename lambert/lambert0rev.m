function out = lambert0rev(sv1, sv2, dt, mu)
%LAMBERT0REV Compute 0 Revolution Lambert Transfer
%
%   Assumptions/Warnings:
%   	1. Method breaks down for very short or long dt values.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sv1 [6x1]or[3x1] Initial State or Position Vector (km)
%       2. sv2 [6x1]or[3x1] Final   State or Position Vector (km)
%       3. dt  [1x1]        Transfer TOF (s)
%       4. mu  [1x1]        Central Body Gravitational Parameter (km3/s2)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. vi  [3x1]        Departure Velocity Vector (km/s)
%       2. vf  [3x1]        Arrival   Velocity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Lambert's Problem" Professor Kate Davis
%       "Lambert-Handout-3.pdf"
%       2. Bate, Muller, and White "Fundamentals of Astrodynamics"
%       3. Vallado "Fundamentals of Astrodynamics and Applications"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    r0_ = sv1(1:3); r0 = norm(r0_);
    rf_ = sv2(1:3); rf = norm(rf_);
    
    nu0 = atan2(r0_(2), r0_(1));
    nuf = atan2(rf_(2), rf_(1));
    dnu = nuf-nu0;
    %dnu = atan2(norm(cross(r0_,rf_)), dot(r0_,rf_));
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

        c2   = 1/2;
        c3   = 1/6;
        phi  = 0;
        phiU = 4*pi^2;
        phiL = -4*pi;

        y = r0 + rf + ((A*(phi*c3 - 1))/sqrt(c2));
        xi = sqrt(y/c2);
        dti = ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));

        N = 0.8;

        while abs(dti-dt) > 1e-5

            y = r0 + rf + ((A*(phi*c3 - 1))/sqrt(c2));

            if (A>0.00) && (y<0.0)
                phi = N*(1/c3)*(1 - ((sqrt(c2)/A)*(r0+rf)));
                [c2, c3] = getc2c3(phi);
                y = r0 + rf + ((A*(phi*c3 - 1))/sqrt(c2));
            end

            xi = sqrt(y/c2);

            dti= ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));


            if dti <= dt
                phiL = phi;
            else
                phiU = phi;
            end


            phi = (phiU + phiL)/2;
            [c2, c3] = getc2c3(phi);

        end


        f = 1-(y/r0);
        dg= 1-(y/rf);
        g = A*sqrt(y/mu);


        v0 = (rf_  - f*r0_)/g;
        vf = (dg*rf_ - r0_)/g;

    end
    
    
    % Outputs
    out    = struct;
    out.vi = v0;
    out.vf = vf;
    

    function [c2, c3] = getc2c3(phi)
        if phi > 1e-6
            c2 = (1.0-cos(sqrt(phi)))/phi;
            c3 = (sqrt(phi)-sin(sqrt(phi)))/sqrt(phi^3);
        elseif phi < -1e-6
            c2 = (1.0-cosh(sqrt(-phi)))/phi;
            c3 = (sinh(sqrt(-phi)) - sqrt(-phi))/(sqrt((-phi)^3));
        else
            c2 = 1/2;
            c3 = 1/6;
        end
    end

end