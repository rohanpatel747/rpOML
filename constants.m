function out = constants()
%CONSTANTS Astrodynamics Constants


% Conversions
out.aukm    = 1.49597870700e8;
out.dayyear = 365.242189;
out.G       = 6.67259e-11;
out.c       = 299792458;


% Gravitational Constants (mu) (km3/s2)
out.Sun.mu     = 1.32712440018e11;
out.Mercury.mu = 22030;
out.Venus.mu   = 3.24858599e5;
out.Earth.mu   = 3.98600433e5;
out.Luna.mu    = 4903;
out.Mars.mu    = 4.28283100e4;
out.Jupiter.mu = 1.266865361e8;
out.Saturn.mu  = 3.7931208e7;
out.Uranus.mu  = 5.7939513e6;
out.Neptune.mu = 6.835100e6;
out.Pluto.mu = 8.71e2;


% Planetary Radii (km)
out.Mercury.r = 2440;
out.Venus.r   = 6051.8;
out.Earth.r   = 6378.1363;
out.Luna.r    = 1737;
out.Mars.r    = 3396.19;
out.Jupiter.r = 71492;
out.Saturn.r  = 60268;
out.Uranus.r  = 25559;
out.Neptune.r = 24764;
out.Pluto.r   = 1188.3;


% Planetary Orbit Semi-Major Axis (km)
out.Mercury.sma = 57.91e6;
out.Venus.sma   = 108.2e6;
out.Earth.sma   = 1.00*out.aukm;
out.Luna.sma    = 384.4e3;
out.Mars.sma    = 1.52368*out.aukm;
out.Jupiter.sma = 778.6e6;
out.Saturn.sma  = 1.433e9;
out.Uranus.sma  = 2.872e9;
out.Neptune.sma = 4.495e9;
out.Pluto.sma   = 5.87e9;


% Planetary Orbital Period (s)
d = 86400;
y = d*out.dayyear;
T = [87.97*d, 224.7*d, 365.256*d, ...
     1.881*y, 11.86*y,  29.460*y, 84.01*y, 164.8*y, 247.7*y];

out.Mercury.t = T(1);
out.Venus.t   = T(2);
out.Earth.t   = T(3);
out.Luna.t    = 27.322*d;
out.Mars.t    = T(4);
out.Jupiter.t = T(5);
out.Saturn.t  = T(6);
out.Uranus.t  = T(7);
out.Neptune.t = T(8);
out.Pluto.t   = T(9);


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

