function out = getAnomalyandDt(k,displayOut)
%GETANOMALYANDDT Gets Eccentric, Mean Anomaly, and DT from Periapsis
%   Detailed explanation goes here

    eccTol = 5;

    %disp(' ');
    %disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
    %disp('getAnomalyandDT : TA (Theta) must be from [0 360] degrees');
    %disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
    %disp(' ');

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


    elseif round(k.e,eccTol) == 1
        % Parabolic

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
    out.ea = EcAnom*(180/pi);
    out.ma = MnAnom*(180/pi);
    out.dt = dtPeri;
    
    out.unitNames = [
        "Ec : Eccentric Anomaly (deg)";
        "Mn : Mean Anomaly (deg)";
        "dt : Time since Periapsis (s) (always positive)";
    ];

end


