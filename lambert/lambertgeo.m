function out = lambertgeo(mu, x1, x2, dt, N, Long, type, fullState)
%LAMBERTGEO Solve Lambert's Problem and output Transfer DVs using Geometric
%Method. (CURRENTLY ONLY FOR ELLIPTICAL ORBITS)
%   
%   Inputs:
%       mu        = [1x1] Central Body Gravitational Parameter (km3/s2)
%       x1        = [6x1] Initial State Vector [r1; v1]  
%       x2        = [6x1] Final State Vector [r2; v2]
%       dt        = [1x1] Transfer Time of Flight (seconds)
%       N         = [1x1] Number of Revolutions (0 revs. only supported right now)
%       Long      = [T/F] True if Transfer ta>pi, False if Transfer ta<pi
%       type      = [str] Conic Type ('elliptical','hyperbolic','parabolic')
%       fullstate = [T/F] Output Boolean
%                           True - Returns Full State Structure of Transfer
%                                  Orbit Includes: 
%                                       fXt.lambertTOF = Transfer TOF (s)
%                                       fXt.lambertDTA = Change in True Anomaly (deg.)
%                                       fXt.lambertN   = Number of Revolutions
%                           False- Returns array with [DV1; DV2] where both
%                                  DV1and2 are [1x3] inertial vectors in km/s

    disp('*** LAMBERTGEOMETRIC(): Currently for Elliptical Orbits Only! ***'); disp(' ');
    
    
    % Process State Vectors
    r1_ = x1(1:3); r1 = norm(r1_);
    v1_ = x1(4:6); v1 = norm(v1_);
    r2_ = x2(1:3); r2 = norm(r2_);
    v2_ = x2(4:6); v2 = norm(v2_);

    
    % Space Triangle Properties
    dta = acos(dot(r1_,r2_)/(r1*r2));   % Change in True Anomaly from r1_ to r2_
    if Long; dta = 2*pi - dta; end      % Looking for trajectory true anomaly > pi? 
    c = norm(r2_ - r1_);                % Chord sqrt((r1^2) + (r2^2) - (2*r1*r2*cosdta))
    s = (1/2)*(r1+r2+c);                % Semi-Perimeter

    
    % Minimum Energy Transfer Arc
    a_min     = s/2;
    alpha_min = pi;
    beta_min  = 2*asin(sqrt((s-c)/(2*a_min)));
                if dta > pi; beta_min = -beta_min; end
    n_min     = sqrt(mu/(a_min^3));
    TOF_min   = (1/n_min)*(2*N*pi + (alpha_min - beta_min) - (sin(alpha_min) - sin(beta_min)));

    
    % Parabolic TOFp check
    if dta < pi; TOFp = (1/3)*(sqrt(2/mu))*((s^(3/2)) - ((s-c)^(3/2)));
           else; TOFp = (1/3)*(sqrt(2/mu))*((s^(3/2)) + ((s-c)^(3/2))); end
    if dt>TOFp;  disp('elliptical');
           else; disp('hyp'); end
           
    
    
    % Iteratively Solve Lambert's Eqn. for SMA (a)
    options = optimoptions('fsolve','Display','none','TolFun', 1e-6);
    a = fsolve(@lambertsEqn,a_min, options);
    [alpha, beta] = getAlphaBeta(a);

    
    % Find Orbit Properties
    p = (((4*a)*(s-r1)*(s-r2))/(c^2))*((sin((alpha + beta)/2))^2);
    e = sqrt(1-(p/a));
    E = -mu/(2*a);
    h = sqrt(mu*p);
    
    
    % Find Transfer DVs using f&g functions
    f = 1 - (r2/p)*(1-cos(dta));
    g = ((r2*r1)/sqrt(mu*p))*sin(dta);
    df= sqrt(mu/p)*tan(dta/2)*(((1-cos(dta))/p)-(1/r2)-(1/r1));
    dg= 1-(r1/p)*(1-cos(dta));

    v1_t = (1/g)*(r2_ - f*r1_);
    v2_t = (1/g)*(dg*r2_ - r1_);
    dv1_ = v1_t - v1_; dv1 = norm(dv1_);
    dv2_ = v2_t - v2_; dv2 = norm(dv2_);

    
    % Outputs
    if fullState
        out = create_state([r1_;v1_t],  mu, 'rv', false, false, true);
        out.lambertTOF = dt;
        out.lambertDTA = dta*(180/pi);
        out.lambertN   = N;
    else
        out = [dv1_';dv2_'];
    end
    

    % Itteratively solve this function using fsolve() to get the SMA (a)
    function dTOF = lambertsEqn(a)
        [alpha, beta] = getAlphaBeta(a);
        n   = sqrt(mu/(a^3));
        TOF = (1/n)*((2*N*pi) + alpha - beta - (sin(alpha) - sin(beta)));
        dTOF= abs(TOF-dt);
    end


    % Find Associated Alpha and Beta Values
    function [alpha, beta] = getAlphaBeta(a)

        
        % Principle Alpha and Beta Values
        if strcmp(type, 'elliptical')
            [alpha_p, beta_p] = apbpElip(a);
        elseif strcmp(type, 'hyperbolic')
            [alpha_p, beta_p] = apbpHyp(a);
        end % parabolic arcs don't depend on alpha and beta

        
        % Conditions for Alpha and Beta
        if dt < TOF_min
            % Lower e (Closer F')
            if dta < pi
                alpha = alpha_p;
                beta  = beta_p;
            else
                alpha = alpha_p;
                beta  = -beta_p;
            end
        else
            % Higher e (Further F'')
            if dta < pi
                alpha = (2*pi) - alpha_p;
                beta  = beta_p;
            else
                alpha = (2*pi) - alpha_p;
                beta  = -beta_p;
            end
        end

    end


    % Principal Alpha and Beta - Elliptical Orbits
    function [alpha_p, beta_p] = apbpElip(a)
        alpha_p = 2*asin(sqrt(s/(2*a)));
        beta_p  = 2*asin(sqrt((s-c)/(2*a)));
    end


    % Principal Alpha and Beta - Hyperbolic Orbits
    function [alpha_p, beta_p] = apbpHyp(a)
       alpha_p = 2*sinh(sqrt((s)/(2*abs(a)))); 
       beta_p  = 2*sinh(sqrt((s-c)/(2*abs(a))));
    end


end

