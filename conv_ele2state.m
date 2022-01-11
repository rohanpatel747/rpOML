function out = conv_ele2state(aeiowta, mu, isDeg, displayOut)
%CONV_ELE2STATE Summary of this function goes here
%   Detailed explanation goes here

    if isDeg
       aeiowta(3:6) = aeiowta(3:6)*(pi/180);
    end
    
    a = aeiowta(1); 
    e = aeiowta(2);
    i = aeiowta(3); 
    o = aeiowta(4); 
    w = aeiowta(5); 
    ta= aeiowta(6);

    p = a*(1 - (e^2));
    r = p / (1 + (e*cos(ta)));

    % Perifocal State Vector
    r_ = [r*cos(ta); r*sin(ta); 0.000];
    v_ = sqrt(mu/p)*[-sin(ta); (e+cos(ta)); 0.00];
    xp = [r_;v_];
    

    % Inertial State Vector
    % 3-1-3 Rotation Sequence using raan, inclination, and arg.peri.
    R3w = [cos(w), sin(w), 0;
          -sin(w), cos(w), 0;
                0,      0, 1];

    R1i = [1,      0,      0;
           0, cos(i), sin(i);
           0,-sin(i), cos(i)];

    R3o = [cos(o), sin(o), 0;
          -sin(o), cos(o), 0;
                0,      0, 1];

    R313 = (R3w*R1i*R3o)';
    xi = [R313*r_; R313*v_];


    if displayOut
        disp(' ')
        disp('-------------------------------------------')
        disp('Input Orbital Elements: [a;e;i;omega;w;TA]')
        disp(' ')
        disp(aeiowta)
        disp(['Grav. Param. CB:     ',num2str(mu), ' (km3/s2)'])
        disp('-------------------------------------------')
        disp('Perifocal State Vector')
        disp(' ')
        disp(xp)
        disp(' ')
        disp('Inertial State Vector')
        disp(' ')
        disp(xi)
        disp('-------------------------------------------')
        disp(' ')
    end
    
    out = struct;
    out.xp = xp;
    out.xi = xi;
    out.R313 = R313;
    out.unitNames = [
        "xp : Perifocal State Vector";
        "xi : Inertial State Vector";
        "R313 : Rotation Matrix";
    ];

end

