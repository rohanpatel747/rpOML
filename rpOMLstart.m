function rpOMLstart(optional)
%RPOMLSTART Initialize rpOMLibrary
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs: 
%       1. 'cpsice',true/false [bool] - Should rpOML load the CPSICE MICE
%                                       toolkit?
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Syntax:
%   (add either line to the top of your script)
%
%       Load only rpOML:            rpOMLstart();
%       Load rpOML and CSPICE:      rpOMLstart('cspice',true);
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   [OPTIONAL JPL CSPICE MICE LIBRARY]
%       Currently, rpOML is not dependent on the JPL CPSICE MICE Toolkit.
%       However, it is a powerful astrodynamics library that can be used
%       for orbital mechanics and space mission design purposes. It can be
%       loaded through this script by setting up the optional path inputs.
%
%       Paths/Variable Inputs:
%           1. deVersion  (./jpl/spk/) Planetary Ephemeris BSP
%               Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
%               Update Char Variable with latest version number.
%           2. pckVersion (./jpl/spk/) Planetary Constants TPC
%               Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
%               Update Char Variable with latest version number.
%           3. tlsVersion (./jpl/spk/) Frames File TLS
%               Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
%               Update Char Variable with latest version number.
%           4. libdir  (./dev/git/rpOML) Library Directory Path
%           5. miceDir (./jpl/mice/)     Path to CSPICE Mice Toolkit
%           6. spkDir  (./jpl/spk/)      Path to MICE Data Files (.bsp, .tpc, .tls)

%
%% INPUTS
arguments
    optional.cspice {mustBeNumericOrLogical} = false;
end

% _________________________________________________________________________
% [REQUIRED] INPUT : RPOML LIBRARY PATH:
if ismac
    libDir = '/Users/rohan/dev/git/rpOML/';
elseif ispc
    libDir = 'C:\Users\rohan\dev\git\rpOML\';
end
% _________________________________________________________________________

% [OPTIONAL] INPUT : JPL CSPICE MICE LIBRARY PATHS/VERSIONS:
if optional.cspice
    % Planetary SPK and Constants Versions
    deVersion = '440'; % cooresponds to the downloaded file: (de440.bsp)
    pckVersion= '10';  % cooresponds to the downloaded file: (naif0010.tls)
    tlsVersion= '12';  % cooresponds to the downloaded file: (pck00012.tpc)

    % Additional Required Paths if CSPICE is to be used.
    if ismac
        miceDir= '/Users/rohan/dev/jpl/mice/';
        spkDir = '/Users/rohan/dev/jpl/spk/';
    elseif ispc
        miceDir= 'C:\Users\rohan\dev\jpl\mice\';
        spkDir = 'C:\Users\rohan\dev\jpl\spk\';
    end
end








%% Initialize Library w/ or w/o CPSICE MICE Toolkit
addpath(genpath(libDir));
disp('rpOML        Loaded');

if optional.cspice
    addpath(genpath(miceDir));
    addpath(genpath(spkDir));
    cspice_kclear;
    deDir  = [spkDir,     'de',  deVersion, '.bsp'];
    tlsDir = [spkDir, 'naif00', tlsVersion, '.tls'];
    pckDir = [spkDir, 'pck000', pckVersion, '.tpc'];
    cspice_furnsh({deDir,tlsDir,pckDir});
    disp('CSPICE MICE  Loaded');
end

disp('- - - - - - - - - - - -');
disp(' ');
end

