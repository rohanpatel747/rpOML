function ssdout = cr3bp_querySSDPOs(sys, family, libPt, branch, opts)
%CR3BP_QUERYSSDPOS Query JPL's Solar System Dynamics Website for PO ICs.
%
% Query Database Link: https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
%__________________________________________________________________________
% Family Options (family) (Type: char):
%   Libration Point Oriented (Specify Lagrange Pt.):
%       1. 'lyapunov'  - Lyapunov
%       2. 'halo'      - Halo
%       3. 'vertical'  - Vertical
%       4. 'axial'     - Axial
%       5. 'longp'     - Long  Period
%       6. 'short'     - Short Period (w/o 'p')
%       7. 'butterfly' - Butterfly
%       8. 'dragonfly' - Dragonfly
%   Secondary Oriented:
%       1. 'dro'       - Distant Retrograde Orbit
%       2. 'dpo'       - Distant Prograde   Orbit
%       3. 'lpo'       - Low     Prograde   Orbit
%   Resonant Orbits:
%       1. 'resonant'  - Resonant
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Libration Point Options (libPt) (Type: char: '1','2','3','4',or'5'):
%       1. Halo              Orbits - 1,2,or3
%       2. Long/Short Period Orbits - 4,or5
%       3. Axial/Vertical    Orbits - 1,2,3,4,or5
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Branch for Orbit Family (branch) (Type: char):
%       1. 'N' or 'S' - Halo, Dragonfly, Butterfly Orbits
%       2. 'E' or 'W' - Low Prograde Orbits
%       3. 'nm'       - 'pq' Integers (ex: 1:2-->'12') Resonant Orbits
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Optional Inputs for Bounds:
%       1. 'minJ' or 'maxJ' - Min/Max Jacobi Constant of Search
%       2. 'minT' or 'maxT' - Min/Max Orbit Period (must specify units)
%       3. 'Tunits'         - Orbit Period Units (s, h, d, TU)
%       4. 'minS' or 'maxS' - Min/Max Orbit Stability.
%__________________________________________________________________________
% Examples
%
%   1) All Earth-Moon Northern Halo Orbits about L2:
%       ssdout = cr3bp_querySSDPOs('earth-moon','halo','2','N')
%                                   
%   2) Earth-Moon DPOs with Period Greater than 10 Hours:
%       ssdout =
%       cr3bp_querySSDPOs('earth-moon','dpo','','','minT',10,'Tunits','h');
%
%

    arguments
        sys;
        family;
        libPt;
        branch;
        opts.minJ    {mustBeNumeric} = NaN;
        opts.maxJ    {mustBeNumeric} = NaN;
        opts.minT    {mustBeNumeric} = NaN;
        opts.maxT    {mustBeNumeric} = NaN;
        opts.Tunits  {mustBeText}    = 'TU';
        opts.minS    {mustBeNumeric} = NaN;
        opts.maxS    {mustBeNumeric} = NaN;
    end

    runssd = true;
    minJ = opts.minJ; maxJ = opts.maxJ;
    minT = opts.minT; maxT = opts.maxT;
    minS = opts.minS; maxS = opts.maxS;


    % Libration Point Conditions
    if contains(family,'halo')
        if isempty(libPt) || contains(libPt,'4') || contains(libPt,'5')
            runssd = false; disp('Aborting Query:')
            disp('    Halo Orbits Require Libration Point Input: 1, 2, or 3');
        end
    end
    if contains(family,'short') || contains(family,'longp')
        if isempty(libPt) || contains(libPt,'3') || contains(libPt,'2') || contains(libPt,'1')
            runssd = false; disp('Aborting Query:')
            disp('    Short/Long Period Orbits Require Libration Point Input: 4 or 5');
        end
    end
    if contains(family,'axial') || contains(family,'vertical')
        if isempty(libPt)
            runssd = false; disp('Aborting Query:')
            disp('    Axial/Vertical Orbits Require Libration Point Input: 1, 2, 3, 4, or 5');
        end
    end
    if contains(family,'dro') || contains(family,'dpo') || contains(family,'lpo')
        libPt = '';
    else
        libPt = ['&libr=',libPt];
    end
    

    % Branch Conditions
    if contains(family,'halo') || contains(family,'dragonfly') || contains(family,'butterfly') 
        if contains(branch,'N') || contains(branch,'S'); else
            runssd = false; disp('Aborting Query:')
            disp('    Halo/Butterfly/Dragonfly Orbits Require Branch Input: N or S');
        end
    end
    if contains(family,'lpo') 
        if contains(branch,'E') || contains(branch,'W'); else
            runssd = false; disp('Aborting Query:')
            disp('    Low Prograde Orbits Require Branch Input: E or W');
        end
    end
    if contains(family,'resonant')
        if contains(branch,{'11';'12';'13';'14';'21';'23';'31';'32';'34';'41';'43'}); else
            runssd = false; disp('Aborting Query:')
            disp('    Resonant Orbits Require Branch Input:');
            disp(['       pq = ','11 ','12 ','13 ','14 ','21 ','23 ','31 ','32 ','34 ','41 ','43'])
        end
    end
    if isempty(branch)
        branch = '';
    else
        branch = ['&branch=',branch];
    end


    % Min/Max Conditions
    if isnan(minJ); jmin=''; else; jmin=['&jacobimin=',num2str(minJ)]; end
    if isnan(maxJ); jmax=''; else; jmax=['&jacobimax=',num2str(maxJ)]; end
    if isnan(minT); tmin=''; else; tmin=['&periodmin=',num2str(minT)]; end
    if isnan(maxT); tmax=''; else; tmax=['&periodmax=',num2str(maxT)]; end
    if isnan(minS); smin=''; else; smin=['&stabmin='  ,num2str(minS)]; end
    if isnan(maxS); smax=''; else; smax=['&stabmax='  ,num2str(maxS)]; end

    if ~isnan(minT) || ~isnan(maxT)
        if ~contains(opts.Tunits,{'s','h','d','TU'})
            runssd = false; disp('Aborting Query:')
            disp('    Invalid Units for Time Bounds. Valid Units are:');
            disp('        s, h, d, TU');
        end
        tunits = ['&periodunits=',opts.Tunits];
    else
        tunits = '';
    end
    
    
    % Create Web Request and Process Data
    url = 'https://ssd-api.jpl.nasa.gov/periodic_orbits.api?';
    url = [url,'sys=',sys,'&family=',family,libPt,branch,jmin,jmax,tmin,tmax,smin,smax,tunits];

    if runssd
        ssd = webread(url);

        if contains('data',fields(ssd))
            for i=1:str2num(ssd.count)
                dta = ssd.data{i};
                data(i,1) = str2num(ssd.data{i}{1});    % State x,y,z,vx,vy,vz
                data(i,2) = str2num(ssd.data{i}{2});
                data(i,3) = str2num(ssd.data{i}{3});
                data(i,4) = str2num(ssd.data{i}{4});
                data(i,5) = str2num(ssd.data{i}{5});
                data(i,6) = str2num(ssd.data{i}{6});
                data(i,7) = str2num(ssd.data{i}{8});    % Period (8 in ssd) 
                data(i,8) = str2num(ssd.data{i}{7});    % Jacobi Constant
                data(i,9) = str2num(ssd.data{i}{9});    % Stability Index
            end

            if class(ssd.limits.period) == 'double'
                limits.period    = [(ssd.limits.period(1));
                                    (ssd.limits.period(2));];
                limits.jacobi    = [(ssd.limits.jacobi(1));
                                    (ssd.limits.jacobi(2));];
                limits.stability = [(ssd.limits.stability(1));
                                    (ssd.limits.stability(2));];
            else
                limits.period    = [str2num(ssd.limits.period{1});
                                    str2num(ssd.limits.period{2});];
                limits.jacobi    = [str2num(ssd.limits.jacobi{1});
                                    str2num(ssd.limits.jacobi{2});];
                limits.stability = [str2num(ssd.limits.stability{1});
                                    str2num(ssd.limits.stability{2});];
            end

            ssdout = struct;
            ssdout.signature      = ssd.signature;
            ssdout.ssdsystem      = ssd.system;
            ssdout.ssdmu          = str2double(ssd.system.mass_ratio);
            ssdout.family         = ssd.family;
            if contains('libration_point',fields(ssd))
            ssdout.librationpoint = ssd.libration_point;
            end
            if contains('branch',fields(ssd))
            ssdout.branch         = ssd.branch;
            end
            if contains('filter',fields(ssd))
            ssdout.filter         = ssd.filter;
            end
            ssdout.limits         = limits;
            ssdout.data           = data;
            ssdout.fields         = 'X, Y, Z, dX, dY, dZ, Per., Jacobi, Stability';
            ssdout.url            = url;

        elseif contains('warning',fields(ssd))
            disp('Query Failed')
            ssdout.url     = url;
        end
    else
        disp('Query Not Configured Correctly.')
        ssdout.url     = url;
    end

end
