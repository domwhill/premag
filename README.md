# ANALYSIS MODULE

### Scripts for plotting premagnetised IMPACT data
Most of the plotting scripts may be found in PLOTTERS_PAPER folder

### Setting internal paths to simulation data
Scripts currently assume test data is located in a folder one level lower than src 
directory entitled '../TF'. 


This may be modified to path on users computer by modifying variable 'data_path' to point to simulation data location.
'data_path' is set in class 'directory_paths' of file MODULES/house_keeping.py.

Names of simulation runs (for different scale lengths/applied field strengths)
 are set in json file MODULES/paths.json.


#Environment
Environment: python2.7


# todo:
check development branch for missing scripts.
- gen_dyqyRL looks like it is more upt otd ate tehere.
- plot_f_log.py - looks like it will be needed.