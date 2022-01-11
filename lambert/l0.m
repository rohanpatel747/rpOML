function l1result = l0(dataout,p1,p2,tof,mu)
%l0 Single Rev Lambert Fit
%       Calculates a single lambert fit
%
%        Input:
%           dataflag = Flag For Requested Data Output
%           p1       = Departure State [6x1] array [x,y,z,vx,vy,vz]
%           p2       = Arrival State [6x1] array [x,y,z,vx,vy,vz]
%           tof      = Time of Flight [1x1] double in seconds
%           mu       = Gravitational Parameter of Central body (km^3/s^2)
%
%        Output Conditions:
%           dataflag = 1
%               l1result(1) = DeltaV_dep mag   
%               l1result(2) = DeltaV_arr mag
%
%           dataflag = 2
%               l1result(1,:) = [deltaVx, deltaVy, deltaVz] @ Departure
%               l1result(2,:) = [deltaVx, deltaVy, deltaVz] @ Arrival
%               l1result(3,:) = [x, y, z] @ Departure
%               l1result(4,:) = [x, y, z] @ Arrival
%               l1result(5,:) = [Vx, Vy, Vz] @ Departure
%               l1result(6,:) = [Vx, Vy, Vz] @ Arrival
%

    r1 = [p1(1),p1(2),p1(3)];
    v1 = [p1(4),p1(5),p1(6)]; 
    %oev1 = eci2orb1(mu, r1, v1);

    r2 = [p2(1),p2(2),p2(3)];
    v2 = [p2(4),p2(5),p2(6)];
    %oev2 = eci2orb1(mu, r2, v2);

    for i = 1:1:3  
        sv1(i) = r1(i);
        sv1(i + 3) = v1(i);
        sv2(i) = r2(i);
        sv2(i + 3) = v2(i);   
    end
    
    nrev = 0;
    [vito, vfto] = glambert(mu, sv1, sv2, tof, nrev);

    % dV vector calcualation
    dvi(1) = vito(1) - v1(1);
    dvi(2) = vito(2) - v1(2);
    dvi(3) = vito(3) - v1(3);
    dvf(1) = vfto(1) - v2(1);
    dvf(2) = vfto(2) - v2(2);
    dvf(3) = vfto(3) - v2(3);
    
    if dataout == 1
        l1result = [norm(dvi),norm(dvf)];
    elseif dataout == 2
        l1result = [dvi;dvf;r1;r2;vito;vfto];
    end
    
end