function dcm = dcm_iot2rth(i,o,ta, w)
%DCM_IOT2RTH Given [radians] inclination, RAAN (omega), TA (theta), Arg.ofPeri (w) find DCM
%   Input angles must be in RADIANS
%   Recall that for the r,t,h frame:
%       dcm(1:3,1) = r/norm(r)  = r_hat
%       dcm(1:3,2) = cross(h,r) = theta_hat
%       dcm(1:3,3) = h/norm(h)  = h_hat

    t = ta+w;

    dcm = [cos(o)*cos(t)-sin(o)*cos(i)*sin(t),	-cos(o)*sin(t)-sin(o)*cos(i)*cos(t),	sin(o)*sin(i);
           sin(o)*cos(t)+cos(o)*cos(i)*sin(t),	-sin(o)*sin(t)+cos(o)*cos(i)*cos(t),   -cos(o)*sin(i);
           sin(i)*sin(t)                     ,   sin(i)*cos(t)                     ,    cos(i)];

end

