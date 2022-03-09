function out = getAnomalyandDt(k,displayOut)
%GETANOMALYANDDT Gets Eccentric, Mean Anomaly, and DT from Periapsis
%
%   Assumptions/Warnings:
%   	1. k.ta (Theta) must be from [0 360] DEGREES
%       2. Parabolic Orbits not Supported
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. k                Input Structure Containing Fields:
%          k.T                  Period (s)
%          k.e                  Eccentricity
%          k.ta                 True Anomaly (degrees)
%          k.mu                 CB Grav Param (km3/s2)
%          k.h                  Angular Momentum (km2/s)
%          k.a                  SMA (km)
%          k.r_                 Position State Vector (km) [3x1]
%       2. displayOut       Boolean to show debug data (true/false)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: 'out' structure with fields:
%       1. ea     [1x1]      Eccentric Anomaly (rad.)
%       2. ma     [1x1]      Mean      Anomaly (rad.)
%       3. dtperi [1x1]      Time to/since Periapsis (s)
%                               Elliptical: always positive
%                               Hyperbolic: (-) Time to   Periapsis
%                                           (+) Time from Periapsis
%       4. taInf  [1x1]      True Anomaly at Infinity (rad.) (hyp. only)
%       5. unitNames         String Array of Output Quantity Units
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    eccTol = 5;                                             % Eccentricity Tolerance
    n = 2*pi/k.T;                                           % Mean Motion (2pi/period)

    if round(k.e,eccTol) < 1
        % Ellipitcal
        EcAnom = acos((1-(norm(k.r_)/k.a))/k.e);            % Eccentric Anomaly (rad)
        if k.ta > 180; EcAnom = -EcAnom; end

        MnAnom = EcAnom - k.e*sin(EcAnom);                  % Mean Anomaly (rad)
        dtPeri = MnAnom/n;                                  % Time since periapsis (s)
        if dtPeri < 0; dtPeri = k.T + dtPeri; end
        
    elseif round(k.e,eccTol) > 1
        % Hyperbolic
        e = k.e;
        ta= k.ta * (pi/180);
        mu= k.mu;
        h = k.h;
        
        tainf = acos(-1/e);
        if ta>pi; tainf = 2*pi - tainf; end
        
        F  = log((sqrt(e+1) + sqrt(e-1) * tan(ta/2))/(sqrt(e+1) - sqrt(e-1) * tan(ta/2)));
        Mh = e*sinh(F) - F;
        dtPeri = (h^3/(mu^2))*(1/((e^2 - 1)^(3/2)))*Mh;
        EcAnom = F;
        MnAnom = Mh;
        hypOrb = true;

    elseif round(k.e,eccTol) == 1
        % Parabolic
        disp('***ERROR***')
        disp('Parabolic Orbit Eccentric, Mean Anomalies, and dt not coded yet');
        disp(' ');
    end

    % Display Outputs
    if displayOut
        disp('-------------------------------------------')
        disp('Input Elements')
        disp(' ')
        disp(['Radius:              ',num2str(norm(k.r_)), ' km'])
        disp(['Semimajor Axis:      ',num2str(k.a), ' km'])
        disp(['Eccentricity:        ',num2str(k.e)])
        disp(['True Anomaly:        ',num2str(k.ta), ' deg'])
        disp(['Orbit Period:        ',num2str(k.T), ' s'])
        disp(' ')
        disp(' ')
        disp('Output Anomaly Properties')
        disp(' ')
        disp(['Eccentric Anomaly:   ',num2str(EcAnom), ' rad'])
        disp(['Mean Anomaly:        ',num2str(MnAnom), ' rad'])
        disp(['dt from Periapsis:   ',num2str(dtPeri), ' s'])
        disp(' ')
        disp('-------------------------------------------')
    end
    
    % Outputs
    out = struct;
    out.ea = EcAnom;
    out.ma = MnAnom;
    out.dt = dtPeri;
    if hypOrb; out.taInf = tainf; end
    
    out.unitNames = [
        "Ec    : Eccentric Anomaly (rad)";
        "Mn    : Mean Anomaly (rad)";
        "dt    : Time since Periapsis (s) (always positive for elliptical orbits)";
        "tainf : True Anomaly at Infinity (rad) for hyp. orbits only";
    ];

end


