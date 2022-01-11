function l2(r1,r2,v1,v2)

mu = 1.32712*10^11;     % km^3/s^2
aunit = 149597870.691;  % km
ticstep = 10;           % 10 Day step
tbcoef = 1;

% Compute Orbital Elements of Departure and Arrival Planets
oevi = eci2orb1(mu,r1,v1);
oevf = eci2orb1(mu,r2,v2);




% Find Departure Body Periods
oev1 = eci2orb1(mu, r1, v1);
period1 = 2 * pi * oev1(1) * sqrt(oev1(1) / mu) / 86400;
xve = oev1(1) / aunit;

oev2 = eci2orb1(mu, r2, v2);
period2 = 2 * pi * oev2(1) * sqrt(oev2(1) / mu) / 86400;
if (oev2(1) > oev1(1))
   xve = oev2(1) / aunit;
end

% Set Up Data Tick Marks
npts1 = fix(period1 / deltat);
npts2 = fix(period2 / deltat);
npts3 = fix((jdate2 - jdate1) / deltat);
[rti, vti] = orb2eci(mu, oevtoi);

end