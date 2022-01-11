function l1result = l1(ctrbdy,trgbdy1,trgbdy2,bdy1et,bdy2et,tof,flag)
%{
    Single Pair Lambert Calculation Call
        Input:
            ctrbdy   = Center Body (1x1) Char
            trgbdy1  = Departure Body (1x1) Char
            trgbdy2  = Arrival Body (1x1) Char
            bdy1et   = Ephemeris Time at Departure (1x1) Double
            bdy2et   = Ephemeris Time at Arrival (1x1) Double
            tof      = Time of Flight (1x1) Double
            flag     = Condition for output data type
        Output:
            l1result = Results from Lambert Calculation (1x6) Cell Array
                
                Cell 1 = Orbital Elements of Departure Body (1x6) Double
                Cell 2 = Orbital Elements of Arrival Body (1x6) Double
                Cell 3 = Delta V Vector
%}
    mu = 1.32712*10^11;     % km^3/s^2

    pb1 = mice_spkezr(trgbdy1, bdy1et, 'J2000', 'NONE', ctrbdy);
    pb2 = mice_spkezr(trgbdy2, bdy2et, 'J2000', 'NONE', ctrbdy);
    p1 = [pb1.state];
    p2 = [pb2.state];

    r1 = [p1(1),p1(2),p1(3)];
    v1 = [p1(4),p1(5),p1(6)]; 
    oev1 = eci2orb1(mu, r1, v1);

    r2 = [p2(1),p2(2),p2(3)];
    v2 = [p2(4),p2(5),p2(6)];
    oev2 = eci2orb1(mu, r2, v2);

    for i = 1:1:3  
        sv1(i) = r1(i);
        sv1(i + 3) = v1(i);
        sv2(i) = r2(i);
        sv2(i + 3) = v2(i);   
    end

    [vito, vfto] = glambert(mu, sv1, sv2, tof, 0);
    oevtoi = eci2orb1(mu, r1, vito');
    oevtof = eci2orb1(mu, r2, vfto');

    % dV vector calcualation
    dvi(1) = vito(1) - v1(1);
    dvi(2) = vito(2) - v1(2);
    dvi(3) = vito(3) - v1(3);
    dvf(1) = vfto(1) - v2(1);
    dvf(2) = vfto(2) - v2(2);
    dvf(3) = vfto(3) - v2(3);
  
    %l1result = {oevtoi,oevtof,dvi,dvf,r1,r2,v1,v2};
    l1results = {dvi,dvf,r1,r2,v1,v2};  %reduced output for saving comp. time
    
end