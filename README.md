# rpOML
## MATLAB Orbital Mechanics Library
Created: 11JAN22
<br /> **NOTE: This library is a work in progress and will be frequently updated** <br />
This library is compiled from work done for several astrodynamics and mission design CU Boulder graduate courses.

## Overview
This library of MATLAB functions and scripts is designed around CU Boulder’s ASEN5050, ASEN6008, and ASEN6060 courses. The intent is to catalog common astrodynamics functions and tools that are used in trajectory design. Optionally, the JPL CSPICE MICE [Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) can be loaded through this library. The toolkit is also useful to process JPL published constants, frames, and ephemerides. rpOML primarily consists of 2-Body orbital mechanics functions and interplanetary mission design related algorithms.

## Library Contents
Please see the rpOML_readme.pdf file for a complete list of the scripts, functions, and examples included in this repository.

* Planetary Constants and Astrodynamics Quantities
* CR3BP Non-Dim. (ND) Characteristic Quantities
* Query Planetary States (from Ephemeris or Meeus Alg.)
* Numerical Propagation EOMs: 2BP, CR3BP, and NBP
* Kepler's Equations for 2BP Propagation
* Lambert's Problem (0 rev and N rev solutions)
* 2BP State Conversions
* CR3BP conversions, Jacobi Constant, STM, and Periodic Orbits
* Launch Vehicle Payload to C3
* Porkchop Plots
* Flybys in 2D, 3D, and Orbit Resonance
* B-Plane Targeting
* Interplanetary Broad Trajectory Search

## Installation
1. Clone the following repository
2. In the download there is a file called rpOMLstart.m . Open this file in MATLAB and follow the commented block of code at the top to set the correct path(s).
3. Extract the "/mice/" folder from Step (2) anywhere you'd like. This folder contains the entire toolkit.
5. [Optional] The JPL CSPICE MICE Toolkit can be loaded in with this library. To do so, follow the instructions in the file from step (2).

## Using the Library
1. Copy/Paste  rpOMLstart.m to your project’s working directory.
2. At the beginning of your script add the line: rpOMLstart();
3. If CSPICE is to be initialized add this line instead: rpOMLstart('cspice',true);

Note: Step (3) will clear the CSPICE kernal and reinitialze it everytime it is called. Consider this when attempting to use "rpOMLstart('cspice',true);" in functions or subroutines. 

Examples of commonly used functions are shown in the "./examples/" directory.
