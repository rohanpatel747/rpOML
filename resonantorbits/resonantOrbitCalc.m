function out = resonantOrbitCalc(in)
%RESONANTORBITCHECK Plots 1 Resonant Orbit (2 Flybys) of a Body 
%
%   Assumptions/Warnings:
%   	1. This script only works for res. orbits w.r.t planets (not moons)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. resbdy [1x1]     Resonant Orbit w.r.t this Planet (ID Number)
%       2. prev   [1x1]     Number of Revolutions by the Planet
%       3. srev   [1x1]     Number of Revolutions by the Spacecraft
%       4. mucb   [1x1]     Central Body Gravitational Parameter (km3/s2)
%       5. tga1   [1x1]     Julian Epoch of Incoming Flyby
%       6. tga2   [1x1]     Julian Epoch of Outgoing Flyby
%       7. vi1_   [3x1]     Pre-Resonant Orbit Incoming Vinfinity Vector (km/s)
%       8. vi2_   [3x1]     Post-Resonant Orbit Outgoing Vinfinity Vector (km/s)
%       9. rminfb [1x1]     Minimum Flyby Radius (km)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: 'out' structure with following fields:
%       1. T      [1x1]     Resonant Orbit Period (s)
%       2. a      [1x1]     Resonant Orbit SMA (km)
%       3. vi1p_  [3x1]     Post 1st Flyby V-Infinity Vector (km/s)
%       4. vi2m_  [3x1]     Pre  2nd Flyby V-Infinity Vector (km/s)
%       5. rca1   [1x1]     1st Flyby Close Approach Radius (km)
%       6. rca2   [1x1]     2nd Flyby Close Approach Radius (km)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. getStatePlanet()
%       2. getPlanetName()
%       3. bplanefromVi1Vi2()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % Inputs
    resbdy = in.resbdy;
    prev   = in.prev;
    srev   = in.srev;
    mucb   = in.mucb;
    tga1   = in.tga1;
    tga2   = in.tga2;
    vinf1_ = in.vi1_;
    vinf2_ = in.vi2_;
    rminfb = in.rminfb;


    % Body States
    xga1 = getStatePlanet(resbdy,tga1).x;
    xga2 = getStatePlanet(resbdy,tga2).x;

    rp1_ = xga1(1:3).';  rp1 = norm(rp1_);
    vp1_ = xga1(4:6).';  vp1 = norm(vp1_);
    rp2_ = xga2(1:3).';  rp2 = norm(rp2_);
    vp2_ = xga2(4:6).';  vp2 = norm(vp2_);


    % V-Infinity Conditioning
    vinf1_= vinf1_(:);
    vinf2_= vinf2_(:);
    vinf1 = norm(vinf1_);
    vinf2 = norm(vinf2_);
    vinf = (vinf1+vinf2)/2;


    % Resonant Orbit Properties and Angle Between Vp and Vinf.
    murb = constants().(getPlanetName(resbdy)).mu;
    T    = constants().(getPlanetName(resbdy)).t;
    T    = T*(prev/srev);                                 % res. orbit period
    a    = (mucb*((T/(2*pi))^2))^(1/3);                   % res. orbit SMA
    vsc  = sqrt(mucb*( (2/rp1) - (1/a)));                 % res. orbit entry speed
    ta   = acos((-vsc^2 + vinf^2 + vp1^2)/(2*vinf*vp1));  % Ang. btwn Vp and Vinf


    % Transform from VNC Frame to Inertial for Incoming Encounter
    V1_ = vp1_/vp1;
    N1_ = cross(rp1_,vp1_)/norm(cross(rp1_,vp1_));
    C1_ = cross(V1_,N1_);
    VNC1= [V1_, N1_, C1_];


    % Use Given Locus Angle to Find B-Plane and Rp Parameters
    phi = in.phi*(pi/180); j=1;
    for i=1:length(phi)

        % Outgoing Vinf From First Flyby
        vinf1pVNC_ = vinf*[cos(pi-ta);
                           sin(pi-ta)*cos(phi(i));
                          -sin(pi-ta)*sin(phi(i))];
        vinf1p_ = VNC1*vinf1pVNC_;
        rca1(i) = bplanefromVi1Vi2(murb, vinf1_, vinf1p_).rp;

        % Incoming Vinf From Second Flyby
        vinf2m_ = vinf1p_ + vp1_ - vp2_;
        rca2(i) = bplanefromVi1Vi2(murb, vinf2m_, vinf2_).rp;

        if (rca1(i)>=rminfb) && (rca2(i)>=rminfb)
           ruseable(j,1:3) = [phi(i)*(180/pi),rca1(i),rca2(i)];
           j=j+1;
        end
    end

    % Outputs
    out = struct;
    out.T = T;
    out.a = a;
    out.entryspeed = vsc;
    out.angVpVinf1 = ta;
    out.VNC = VNC1;
    out.phi = phi;
    out.vi1p_ = vinf1p_;
    out.rca1  = rca1;
    out.vi2m_ = vinf2m_;
    out.rca2  = rca2;
    
end