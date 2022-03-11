function out = getInterplanetaryLambertInfo(depBody, arrBody, depDate, arrDate, opts)
%GETINTERPLANETARYLAMBERTINFO Returns Transfer Parameters btwn. Planets
%
%   Assumptions/Warnings:
%   	1. 0 Rev (Type I/II Transfers) Only
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. depBody    [1x1]             Departure Planet Number (1-9) 
%       2. arrBody    [1x1]             Arrival   Planet Number (1-9)
%       3. depDate    [1x1]             Departure Julian Date OR
%                     ['dd-mmm-yyyy']   Character Format
%       4. arrDate    [1x1]             Arrival   Julian Date OR
%                     ['dd-mmm-yyyy']   Character Format
%   [OPT]. disp       [T/F]             Print Data to Command Window
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. x1     [6x1]        Departure Body State Vector (km, km/s)
%       2. x2     [6x1]        Arrival   Body State Vector (km, km/s)
%       3. vi1_   [3x1]        Departure V-Infinity Vector (km/s)
%       4. vi2_   [3x1]        Arrival   V-Infinity Vector (km/s)
%       5. vi1    [1x1]        Departure V-Infinity Magnitude (km/s)
%       6. vi2    [1x1]        Arrival   V-Infinity Magnitude (km/s)
%       7. tof    [1x1]        Time of Flight (s)
%       8. depc3  [1x1]        Departure C3 Energy (km2/s2)
%       9. deprla [1x1]        Departure R.A. of Launch Asy. (deg.)
%      10. depdla [1x1]        Departure D.A. of Launch Asy. (deg.)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. MATLAB: juliandate()
%       2. getStatePlanet()
%       3. lambert0rev()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % Arguments
    arguments
       depBody;
       arrBody;
       depDate;
       arrDate;
       opts.disp {mustBeNumericOrLogical} = true;
    end

    
    % Constants
    mu = 132712440018;  % from constants()

    
    % Accept Either Julian Date or Calendar Date
    if ischar(depDate)
        jddep = juliandate(depDate,'dd-mmm-yyyy');
    elseif isnumeric(depDate)
        jddep = depDate;
    end
    if ischar(arrDate)
        jdarr = juliandate(arrDate,'dd-mmm-yyyy');
    elseif isnumeric(arrDate)
        jdarr = arrDate;
    end  

    
    % Lambert's Problem
    xd = getStatePlanet(depBody, jddep, 'meeus').x;
    xa = getStatePlanet(arrBody, jdarr, 'meeus').x;

    rd  = xd(1:3);
    vd  = xd(4:6);
    ra  = xa(1:3);
    va  = xa(4:6);
    tof = (jdarr - jddep)*86400;

    l = lambert0rev(rd, ra, tof, mu);
    
    
    % Outgoing/Incomig Vinfinity 
    vi1_ = l.vi - vd;
    vi2_ = l.vf - va;
    
    viox = vi1_(1);
    vioy = vi1_(2);
    vioz = vi1_(3);

    
    % Output Data
    out         = struct;
    out.x1      = xd;
    out.x2      = xa;
    out.vi1_    = vi1_;    out.vi1 = norm(vi1_);
    out.vi2_    = vi2_;    out.vi2 = norm(vi2_);
    out.tof     = tof;
    out.depc3   = (out.vi1)^2;
    out.deprla  = atan2(vioy,viox) * (180/pi);
    out.depdla  = asin(vioz/out.vi1) * (180/pi);

    
    % Print to Command Window
    if opts.disp
        disp(' ');
        disp('_____________________________________________________________');
        disp(['Transfer Between Bodies: ',num2str(depBody), ' - ', num2str(arrBody)])
        
        disp(' ');
        disp(['    Dep Date   : ', datestr(datetime(jddep,'convertfrom','juliandate')), ...
            '      ', num2str(jddep,12)]);
        disp(['    Arr Date   : ', datestr(datetime(jdarr,'convertfrom','juliandate')), ...
            '      ', num2str(jdarr,12)]);
        disp(['    TOF (days) : ', num2str(out.tof/86400,12)]);
        
        disp(' ');
        disp('    V-Infinity:')
        disp(['        Outgoing (km/s) : ', num2str(out.vi1,12)]);
        disp(out.vi1_);
        disp(['        Incoming (km/s) : ', num2str(out.vi2,12)]);
        disp(out.vi2_);
        
        if depBody == 3
            disp(' ');
            disp(['    Dep. C3  (km2/s2) :   ', num2str(out.depc3,12)]);
            disp(['    Dep. RLA (deg)    :   ', num2str(out.deprla,12)]);
            disp(['    Dep. DLA (deg)    :   ', num2str(out.depdla,12)]);
        end
        disp('_____________________________________________________________');
        disp(' ');
    end

end