function rpOMLstart()
%RPOMLSTART Initialize rpOMLibrary
%   Inputs (none)
%   Paths/Variable Inputs:
%       1. deVersion  (./jpl/spk/) Planetary Ephemeris BSP
%          Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
%          Update Char Variable with latest version number.
%       2. pckVersion (./jpl/spk/) Planetary Constants TPC
%          Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/
%          Update Char Variable with latest version number.
%       3. tlsVersion (./jpl/spk/) Frames File TLS
%          Download Latest Here: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/
%          Update Char Variable with latest version number.
%
%       4. libdir  (./dev/git/rpOML) Library Directory Path
%       5. miceDir (./jpl/mice/)     Path to CSPICE Mice Toolkit
%       6. spkDir  (./jpl/spk/)      Path to MICE Data Files (.bsp, .tpc, .tls)

%% INPUTS
% Planetary SPK and Constants Versions
deVersion = '440';  % cooresponds to (de440.bsp)
pckVersion= '10';   % cooresponds to (naif0010.tls)
tlsVersion= '12';   % cooresponds to (pck00012.tpc)

% Required Paths
if ismac
    libDir = '/Users/rohan/dev/git/mice/';
    miceDir= '/Users/rohan/dev/jpl/mice/';
    spkDir = '/Users/rohan/dev/jpl/rpOML/';
elseif ispc
    libDir = 'C:\Users\rohan\dev\git\rpOML\';
    miceDir= 'C:\Users\rohan\dev\jpl\mice\';
    spkDir = 'C:\Users\rohan\dev\jpl\spk\';
end


%% Initialize Library and MICE Toolkit

% Add Directories to Path
addpath(genpath(libDir));
addpath(genpath(miceDir));
addpath(genpath(spkDir));

% Clear/Load CSPICE MICE Library
cspice_kclear;

deDir  = [spkDir,     'de',  deVersion, '.bsp'];
tlsDir = [spkDir, 'naif00', tlsVersion, '.tls'];
pckDir = [spkDir, 'pck000', pckVersion, '.tpc'];
cspice_furnsh({deDir,tlsDir,pckDir});

% Print Statements
disp('rpOML        Loaded');
disp('CSPICE MICE  Loaded');
disp('- - - - - - - - - - - -');
disp(' ');
end

