function [out] = conv_state2ele(xi,mu,displayOut)
%CONV_STATE2ELE Summary of this function goes here
%   Detailed explanation goes here

    r_ = xi(1:3);
    v_ = xi(4:6);
    vr = dot(v_,r_)/norm(r_);
    h_ = cross(r_,v_);
    h = norm(h_);
    i = acos(h_(3)/h)*(180/pi);
    if (90<i) && (i<270)
        i_ = 'retrograde';
    else
        i_ = 'prograde';
    end
    N = cross([0;0;1],h_);
    omega = acos(N(1)/norm(N))*(180/pi);
    if N(2) < 0
        omega = 360-omega;
    end
    evec = (1/mu)*((norm(v_)^2 - mu/norm(r_))*r_ - norm(r_)*vr*v_);
    emag = sqrt(dot(evec,evec));
    w = acos(dot(N,evec)/(norm(N)*emag))*(180/pi);
    if evec(3) < 0
        w = 360-w;
    end
    theta = acos(dot(evec,r_)/(emag*norm(r_)))*(180/pi);
    if dot(r_,v_) < 0
        theta = 360-theta;
    end
    rp = (h^2/mu)*(1/(1+emag));
    ra = (h^2/mu)*(1/(1-emag));
    a = 0.5*(rp+ra);
    T = ((2*pi)/sqrt(mu))*a^(3/2);
    T_hrs = T/3600;
    E = ((norm(v_)^2)/2) - mu/norm(r_);
    p = a*(1-emag^2);
    
    
    if isnan(omega)
        disp(' ');
        disp('* * * C A U T I O N * * *')
        disp('    RAAN is NaN. Value set to 0.00 instead.')
        omega = 0.000;
    end
    
    if isnan(w)
        disp('    Argument of Periapsis is NaN. Value set to 0.00 instead.')
        disp(' ');
        disp(' ');
        w = 0.000;
    end
    
    
    
    if displayOut
        disp('-------------------------------------------')
        disp('Input State Vector: [x;y;z;vx;vy;vz]')
        disp(' ')
        disp(xi)
        disp(['Grav. Param. CB:     ',num2str(mu), ' (km3/s2)'])
        disp('-------------------------------------------')
        disp(['Energy:              ',num2str(E), ' km^2/s^2'])
        disp(['Angular Momentum:    ',num2str(h), ' km^2/s'])
        disp(' ')
        disp(['Inclination:         ',num2str(i), ' deg'])
        disp(['Inclination Type:    ',i_])
        disp(' ')
        disp(['RAAN:                ',num2str(omega), ' deg'])
        disp(['Eccentricity:        ',num2str(emag)])
        disp(['Arg. of Periapsis:   ',num2str(w), ' deg'])
        disp(['True Anomaly:        ',num2str(theta), ' deg'])
        disp(' ')
        disp(['Periapsis Radius:    ',num2str(rp), ' km'])
        disp(['Apoapsis Radius:     ',num2str(ra), ' km'])
        disp(['Semimajor Axis:      ',num2str(a), ' km'])
        disp(['Semilatus Rectum:    ',num2str(p), ' km'])
        disp(['Orbit Period:        ',num2str(T_hrs), ' hrs'])
        disp('-------------------------------------------')
    end
    
    out = struct;
    out.xi = xi;
    out.mu = mu;
    out.r_ = r_;
    out.v_ = v_;
    out.h_ = h_;
    out.h  = h;
    out.a  = a;
    out.e_ = evec;
    out.e  = emag;
    out.i  = i;
    out.o = omega;
    out.w  = w;
    out.ta = theta;
    out.rp = rp;
    out.ra = ra;
    out.T = T;
    out.orbitType = i_;
    out.nodesVec = N;
    out.E = E;
    out.p = p;

    out.units = [
        "length : km";
        "velocity : km/s";
        "angular momentum : km2/s";
        "angle : deg";
        "time : seconds";
        "energy : km2/s2";
    ];
    out.unitNames = [
        "a : semi-major axis";
        "e : eccentricity";
        "i : inclination";
        "omega : RAAN";
        "w : argument of periapsis";
        "theta : True Anomaly";
        "p : Semi-latus Rectum";
        "rp : Radius at Periapsis";
        "ra : Radius at Apoapsis";
        "T : Period";
        "nodesVec : Node Vector";
    ];
    
end