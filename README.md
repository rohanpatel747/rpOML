# rpOML
## MATLAB Orbital Mechanics Library
C: 11JAN22  <br />
**NOTE: This library is a work in progress and will be frequently updated** <br />
Library designed around CU Boulder's ASEN5050 and ASEN6008 courses.

## Overview
This library of MATLAB functions and scripts is designed around CU Boulder’s ASEN5050 and ASEN6008 courses. The intent is to catalog common astrodynamics functions and tools that are used in trajectory design and orbital mechanics. Optionally, the JPL CSPICE MICE [Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html) can be loaded through this library. Currently, rpOML is not dependent on CSPICE, but in future development, it is likely that astrodynamics functions from JPL’s library will be utilized here. The toolkit is also useful to process JPL published constants, frames, and ephemerides. rpOML primarily consists of 2-Body orbital mechanics functions and interplanetary mission design related algorithms.

## Installation
1. Clone the following repository
2. In the download there is a file called rpOMLstart.m . Open this file in MATLAB and follow the commented block of code at the top to set the correct path(s).
3. Extract the "/mice/" folder from Step (2) anywhere you'd like. This folder contains the entire toolkit.
5. [Optional] The JPL CSPICE MICE Toolkit can be loaded in with this library. To do so, follow the instructions in the file from step (2).

## Using the Library
1. Copy/Paste  rpOMLstart.m to your project’s working directory.
2. At the beginning of your script add the line: rpOMLstart();

Note: Step (2) will clear the CSPICE kernal and reinitialze it everytime it is called. Consider this when attempting to use "rpOMLstart();" in functions or subroutines. 

Examples of commonly used functions are shown in the "_examples.m" script file.


## Library Contents
Please see the rpOML_readme.pdf file for a list of the scripts and functions included in this repository.
