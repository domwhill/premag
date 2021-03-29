# Premagnetised paper plots 
## Setup instructions
Dependencies:
- [anaconda](https://www.anaconda.com/products/individual#Downloads)
- [gnu make] - this is supplied with most operating systems

1. Copy (or sym link) the input data to a folder in the root directory entitled `data`.

2. In the commandline of the root directory enter:
 ```make all``` 
 This should build the conda based python 2.7 environment and generate the plots required. 
 Figures will be saved to a folder entitled `output_images/`.

# Notes on simulation data structure
Scripts expect simulation data to be stored in folder `data/` of root of the repo.
Lookup data for simulations of diffenet magnetic field strengths/scale lengths is contained in modules/paths.json file.
Folders are organised first by 1D and 2D simulations (1D_RUNS, 2D_RUNS); then by scale length (LT1, LT2, etc.)
 LT1 = 20um conduction zone. LT2 = 40um conduction zone; finally by magnetic field strength.

# Generating figures
Type the following command in the commandline in root directory:

```make build_plots```

