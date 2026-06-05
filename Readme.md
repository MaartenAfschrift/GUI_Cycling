# GUI Spring Cycling

This GUI serves as an introductory example of how musculoskeletal modeling and simulation can help understand cycling movement and support bike design optimization.

The cycling model is based on [[Limiting radial pedal forces greatly reduces maximal power output and efficiency in sprint cycling; an optimal control study | Journal of Applied Physiology | American Physiological Society](https://doi.org/10.1152/japplphysiol.00733.2021)] (but runs here using the VU the KUL muscle model) (https://doi.org/10.1007/s10439-016-1591-9).

You can run this in matlab (origin implementation) or in python. Note that the matlab2python conversion was done by github copilot and not extensively checked:

## Matlab Installation

1. Download or clone this repository.

2. *(Optional)* Install CasADi in MATLAB: [https://web.casadi.org/](https://web.casadi.org/)  
   If you skip this step, the script will attempt to install the correct CasADi version for you.

3. Run `CyclingApp_exported.m`.

## Python installation

create conda environment:

```bash
conda env create -f environment.yml
conda activate cycling_gui
```

go to the instllation folder. In my case:

```bash
cd C:\Users\mat950\sim\cycling_gui
```

run the python script

```bash
python -m cycling_gui 
```



## Notes

This project is primarily intended as a **tutorial or showcase**. It is not suitable for research purposes in its current form. To reduce computation time, several shortcuts have been implemented—such as a reduced number of mesh points and relaxed solver tolerances—which come at the cost of simulation accuracy.

Enjoy exploring!
