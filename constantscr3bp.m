function sys = constantscr3bp()
%CONSTANTSCR3BP Returns Set of Constants Used for the CR3BP
%               Including mu and non-dimm to dim. parameters
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs: (none)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: Structure 'out' Containing Fields:
%       1. Fieldname is the CR3BP sys. name with subfields:
%           a) mu  - Secondary Mass Parameter (use this)
%           b) mu1 - Primary   Mass Parameter
%           c) L   - Length    Conversion (km)
%           d) V   - Velocity  Conversion (km/s)
%           e) T   - Time      Conversion (s)
%           f) b1  - State of Primary   Body (km and km/s) [NON-DIM.]
%           g) b2  - State of Secondary Body (km and km/s) [NON-DIM.]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Description:
%       Conversion Rules from Non-Dim. to Dim. System:
%           Distance: d' = L*d
%           Velocity: s' = V*s
%           Time: t' = t*(T/(2*pi))
%           Mass Parameter: mu = m2/(m1+m2) assuming m1>m2
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. System Parameters (Lo Section 2.2, p.25_bk)
%

    c = constants();



    sys.SunJupiter.mu      = 9.537e-4;
    sys.SunEarthMoon.mu    = 3.04042340206596e-06;
    sys.EarthMoon.mu       = 0.012150585609624;
    sys.MarsPhobos.mu      = 1.667e-08;
    sys.JupiterIo.mu       = 4.704e-05;
    sys.JupiterEuropa.mu   = 2.528e-05;
    sys.JupiterGanymede.mu = 7.804e-05;
    sys.JupiterCallisto.mu = 5.667e-05;
    sys.SaturnMimas.mu     = 6.723e-08;
    sys.SaturnTitan.mu     = 0.0002366;
    sys.NeptuneTriton.mu   = 0.0002089;
    sys.PlutoCharon.mu     = 0.1097;

    sys.SunJupiter.L       = 7.784e8;
    sys.SunEarthMoon.L     = 149597870.691; 
    sys.EarthMoon.L        = 384747.962856037;
    sys.MarsPhobos.L       = 9380;
    sys.JupiterIo.L        = 421800;
    sys.JupiterEuropa.L    = 671100;
    sys.JupiterGanymede.L  = 1070000;
    sys.JupiterCallisto.L  = 1883000;
    sys.SaturnMimas.L      = 185600;
    sys.SaturnTitan.L      = 1222000;
    sys.NeptuneTriton.L    = 354800;
    sys.PlutoCharon.L      = 19410;

    sys.SunJupiter.V       = 13.102;
    sys.SunEarthMoon.V     = 29.784;
    sys.EarthMoon.V        = 1.025;
    sys.MarsPhobos.V       = 2.144;
    sys.JupiterIo.V        = 17.39;
    sys.JupiterEuropa.V    = 13.78;
    sys.JupiterGanymede.V  = 10.909;
    sys.JupiterCallisto.V  = 8.226;
    sys.SaturnMimas.V      = 14.367;
    sys.SaturnTitan.V      = 5.588;
    sys.NeptuneTriton.V    = 4.402;
    sys.PlutoCharon.V      = 0.222;

    sys.SunJupiter.T       = 3.733e8;
    sys.SunEarthMoon.T     = 3.147e7;
    sys.EarthMoon.T        = 2361000;
    sys.MarsPhobos.T       = 27490;
    sys.JupiterIo.T        = 152400;
    sys.JupiterEuropa.T    = 306000;
    sys.JupiterGanymede.T  = 616500;
    sys.JupiterCallisto.T  = 1438000;
    sys.SaturnMimas.T      = 81170;
    sys.SaturnTitan.T      = 1374000;
    sys.NeptuneTriton.T    = 506400;
    sys.PlutoCharon.T      = 550300;
    
    sys.EarthMoon.rb1      = c.Earth.r   / sys.EarthMoon.L;
    sys.EarthMoon.rb2      = c.Luna.r    / sys.EarthMoon.L;
    sys.SunEarthMoon.rb1   = c.Sun.r     / sys.SunEarthMoon.L;
    sys.SunEarthMoon.rb2   = c.Earth.r   / sys.SunEarthMoon.L;
    sys.SunJupiter.rb1     = c.Sun.r     / sys.SunJupiter.L;
    sys.SunJupiter.rb2     = c.Jupiter.r / sys.SunJupiter.L;
    

    fns = fieldnames(sys);
    for i=1:length(fns)
        sys.(fns{i}).mu1 = 1 - sys.(fns{i}).mu;
        sys.(fns{i}).b1  = [ - sys.(fns{i}).mu; 0.00; 0.00 ;0.00 ;0.00; 0.00];
        sys.(fns{i}).b2  = [1- sys.(fns{i}).mu; 0.00; 0.00 ;0.00 ;0.00; 0.00];
        sys.(fns{i}).name= fns{i};
    end

end