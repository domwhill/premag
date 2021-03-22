#!/bin/bash
source activate pop_env

# Fig. 1 lineouts from IMPACT profiles
python PLOTTERS_PAPER/plot_onedim_2col_POP.py

# Fig 2 f0 

# Fig 3 a-d
# 3a
python PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py 'SH x' 'Te' $save_dir
# 3b
python PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py 'vN x' 'Te' $save_dir
# 3c
python PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py 'RL y' 'Bz' $save_dir
# 3d
python PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py 'E y' 'Bz' $save_dir

# Fig. 4 - dT along x and y at different point in conduction zone
python PLOTTERS_PAPER/plot_dT_4plot.py
# Fig. 5 - Hall-parameter contour plot + Righi-Leduc heating 

# Fig. 6 - Hall-parameter contour plot + self-generated magnetic field generation

# Fig 7 - dqy lineouts