function out = cspice_queryStateSingle(bdyID, ti, opt)
%CSPICE_QUERYSTATE Computes States & Times (JD, Cal, ET) of .bsp Traj. File
%
%   Assumptions/Warnings:
%   	1. Does not consider light time corrections.
%       2. Default Assume Values (can be overriden)
%       2.      'frame' : MJ2000 Ecliptic
%       3.     'center' : "Solar System Barycenter"
%       4.     'dtDays' : 10 (10 day spacing) 
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. bdyID [1x1 Char] NAIF Body ID (ex. '-660' <-- common MALTO out) 
%       2. ti               Start Time (Julian Date # or Cal. 'DD-MMM-YYYY')
%       3. tf               End   Time (Julian Date # or Cal. 'DD-MMM-YYYY')
%       4. Optional Inputs
%          'frame'  [char]   Reference Frame (default "ECLIPJ2000")
%                                   Recall "J2000" = "ICRF" in SPICE
%       5. 'ctrbdy' [char]   Center Body of Frame (default S.S.Barycenter)
%       6. 'dtDays' [1x1]    Timestep (days)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%         1. bdyID         [1x1]   Queried Body ID Number
%       2/3. eti/etf       [1x1]   Start/Stop Ephemeris Time (Since J2000)
%       4/5. tijd/tfjd     [1x1]   Start/Stop Julian Date
%       6/7. tical/tfcal   [chr]   Start/Stop Cal. Date (CSPICE format)
%         8. ctrbdy        [chr]   Name of Center Body
%         9. frame         [chr]   Name of Frame
%        10. x             [nx6]   State Vectors Matrix (x,y,z,vx,vy,vz)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/et2utc_c.html
%       2. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html#Julian%20Dates
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

arguments
    bdyID
    ti
    opt.dtDays {mustBeNumeric} = 10;
    opt.ctrBdy {mustBeTextScalar} = 'SOLAR SYSTEM BARYCENTER';
    opt.frame  {mustBeTextScalar} = 'ECLIPJ2000';
end

bdyID = num2str(bdyID);
ctrbdy = opt.ctrBdy;
frame  = opt.frame;
dtDays = opt.dtDays;


% Accept Either Julian Date or Calendar Date
if ischar(ti)
    tijd = juliandate(depDate,'dd-mmm-yyyy');
elseif isnumeric(ti)
    tijd = ti;
end


% Convert to Ephemeris Time
eti = cspice_str2et(['JD',num2str(tijd,10)]);
et(1,:) = eti;


% Query States
sc1 = mice_spkezr(bdyID, et.', frame, 'NONE', ctrbdy);
sc  = [sc1.state];


% Outputs
out        = struct;
out.bdyID  = bdyID;
out.tiet   = eti;
out.tijd   = tijd;
out.tical  = cspice_et2utc( eti, 'C', 3);
out.ctrBdy = ctrbdy;
out.frame  = frame;
out.x(:,1) = sc(1,:);
out.x(:,2) = sc(2,:);
out.x(:,3) = sc(3,:);
out.x(:,4) = sc(4,:);
out.x(:,5) = sc(5,:);
out.x(:,6) = sc(6,:);

end

