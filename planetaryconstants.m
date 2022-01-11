function out = planetaryconstants()
%CONSTANTS Planetary Gravitational Coefficients and Radii


% Conversions
out.aukm    = 1.49597870700e8;
out.dayyear = 365.242189;
out.G       = 6.67259e-11;
out.c       = 299792458;


% Gravitational Constants (mu) (km3/s2)
out.muSun     = 1.32712440018e11;
out.muMercury = 22030;
out.muVenus   = 3.24858599e5;
out.muEarth   = 3.98600433e5;
out.muLuna    = 4903;
out.muMars    = 4.28283100e4;
out.muJupiter = 1.266865361e8;
out.muSaturn  = 3.7931208e7;
out.muUranus  = 5.7939513e6;
out.muNeptune = 6.835100e6;
out.muPluto = 8.71e2;


% Planetary Radii (km)
out.rMercury = 2440;
out.rVenus   = 6051.8;
out.rEarth   = 6378.14;
out.rLuna    = 1737;
out.rMars    = 3396.19;
out.rJupiter = 71492;
out.rSaturn  = 60268;
out.rUranus  = 25559;
out.rNeptune = 24764;
out.rPluto   = 1188.3;


% Planetary Orbit Semi-Major Axis (km)
out.smaMercury = 57.91e6;
out.smaVenus   = 108.2e6;
out.smaEarth   = 149.6e6;
out.smaLuna    = 384.4e3;
out.smaMars    = 227.9e6;
out.smaJupiter = 778.6e6;
out.smaSaturn  = 1.433e9;
out.smaUranus  = 2.872e9;
out.smaNeptune = 4.495e9;
out.smaPluto   = 5.87e9;


% Planetary Orbital Period (s)
d = 86400;
y = d*out.dayyear;
T = [87.97*d, 224.7*d, 365.256*d, ...
     1.881*y, 11.86*y,  29.460*y, 84.01*y, 164.8*y, 247.7*y];

out.tMercury = T(1);
out.tVenus   = T(2);
out.tEarth   = T(3);
out.tLuna    = 27.322*d;
out.tMars    = T(4);
out.tJupiter = T(5);
out.tSaturn  = T(6);
out.tUranus  = T(7);
out.tNeptune = T(8);
out.tPluto   = T(9);


% Units
out.units.r   = "km";
out.units.sma = "km";
out.units.t   = "s";
out.units.mu  = "km3/s2";
out.units.aukm= "km";
out.units.dyr = "days";
out.units.G   = "(N m2)/kg2";
out.units.c   = "km/s";


% Warnings
disp('- - - Caution - - -')
disp('SMA and T taken from Curtis');
disp('Mercury and Moon Info Not Given by Professor');
disp('- - - - - - - - - -')
disp(' ');
end

