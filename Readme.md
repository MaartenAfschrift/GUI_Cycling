# GUI Spring Cycling

This GUI serves as an introductory example of how musculoskeletal modeling and simulation can help understand cycling movement and support bike design optimization.

The cycling model is based on [[Limiting radial pedal forces greatly reduces maximal power output and efficiency in sprint cycling; an optimal control study | Journal of Applied Physiology | American Physiological Society](https://doi.org/10.1152/japplphysiol.00733.2021)] and runs using either the VU muscle model or the KUL muscle model (https://doi.org/10.1007/s10439-016-1591-9).

## Installation

1. Download or clone this repository.

2. *(Optional)* Install CasADi in MATLAB: [https://web.casadi.org/](https://web.casadi.org/)  
   If you skip this step, the script will attempt to install the correct CasADi version for you.

3. Run `CyclingApp_exported.m`.

## Notes

This project is primarily intended as a **tutorial or showcase**. It is not suitable for research purposes in its current form. To reduce computation time, several shortcuts have been implemented—such as a reduced number of mesh points and relaxed solver tolerances—which come at the cost of simulation accuracy.

Enjoy exploring!
