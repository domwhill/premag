# Premagnetised paper plots 
## Setup instructions
Dependencies:
- [anaconda](https://www.anaconda.com/products/individual#Downloads)
- [gnu make] - this is supplied with most operating systems

1. Add input data to folder in root entitled `data/`
2. Build python environment.  Navigate to the root directory of this repo in command line and execute the following command:
```
make install_environment
```
this will build the conda based python2.7 environemnt required for generating the figures.

# Generating figures
Type the following command in the commandline in root directory:

```make build_plots```

this should run python scripts to build plots which are saved to folder `output_images/`.