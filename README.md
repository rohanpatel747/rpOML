# rpOML
## Orbital Mechanics Library



## Installation
1. Clone the following repository
2. Download JPL's NAIF CSPICE MICE Toolkit found [here](https://naif.jpl.nasa.gov/naif/toolkit.html).
3. Extract the "/mice/" folder from Step (2) anywhere you'd like. This folder contains the entire toolkit.
5. Open the "rpOMLstart_copy.m" file and follow the instructions on downloading NAIF data files, setting data file versions, and paths for the library.

## Using the Library
1. Copy/Paste "rpOMLstart_copy.m" to your project's working directory and rename it to "rpOMLstart.m"
2. At the beginning of your script add the line: rpOMLstart();

Note: Step (2) will clear the CSPICE kernal and reinitialze it everytime it is called. Consider this when attempting to use "rpOMLstart();" in functions or subroutines. 

Examples of commonly used functions are shown in the "examples.m" script file.


## Library Contents
### Scripts
| Name | Description |
| ----------- | ----------- |
| ./LVPerf.m  | Various Launch Vehicle Curves (C3 versus Payload Mass)

### Functions
| Name | Description |
| ----------- | ----------- |
| ./create_state.m | Given a state vector OR Kepler's elements, create a "full state" structure (auto-calculates the other set of coordinates) |
| ./dcm_iot2rth.m | Given orbital element angles, find the Directional Cosine Matrix |
| ./conv_ele2state.m | Converts Keplerian elements to inertial state vector | 
| ./conv_state2ele.m | Converts inertial state vector to Keplerian elements |
| ./twobp/prop2bp.m | Integrate 2BP EOMs (optional plotting available) |
| ./twobp/propKepElip.m | (Elliptical) Find new eccentric and true anamonlies given full state IC and DT |
| ./twobp/propKepHyp.m | (Hyperbolic)  Find new inertial state given full state IC and DT |

(more sub-functions are present that can be called on independently)
