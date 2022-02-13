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
%   TODO:
%       1. Itteration counter to prematurely terminate the solver if
%          solution cannot be found
%
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
        psi  = 0;
        psiU = 4*pi^2;
        psiL = -4*pi;

        y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));
        xi = sqrt(y/c2);
        dti = ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));

        N = 0.8;

        while abs(dti-dt) > 1e-5

            y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));

            if (A>0.00) && (y<0.0)
                psi = N*(1/c3)*(1 - ((sqrt(c2)/A)*(r0+rf)));
                [c2,c3] = lambert_getc2c3(psi);
                y = r0 + rf + ((A*(psi*c3 - 1))/sqrt(c2));
            end

            xi = sqrt(y/c2);

            dti= ((xi^3)*c3 + A*sqrt(y))/(sqrt(mu));


            if dti <= dt
                psiL = psi;
            else
                psiU = psi;
            end


            psi = (psiU + psiL)/2;
            [c2,c3] = lambert_getc2c3(psi);

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
    

%     function [c2, c3] = getc2c3(psi)
%         if psi > 1e-6
%             c2 = (1.0-cos(sqrt(psi)))/psi;
%             c3 = (sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
%         elseif psi < -1e-6
%             c2 = (1.0-cosh(sqrt(-psi)))/psi;
%             c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/(sqrt((-psi)^3));
%         else
%             c2 = 1/2;
%             c3 = 1/6;
%         end
%     end

end