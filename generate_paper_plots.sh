#!/bin/bash
source activate pop_env
save_dir=output
mkdir -p $save_dir

# Fig. 1 lineouts from IMPACT profiles
python PLOTTERS_PAPER/plot_onedim_2col_POP.py

# Fig 2 f0 
python PLOTTERS_PAPER/plot_f_log.py
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
python PLOTTERS_PAPER/plot_dqy_contour.py

# Fig. 6 - Hall-parameter contour plot + self-generated magnetic field generation
python PLOTTERS_PAPER/plot_dqy_contour.py -v "bier"
# Fig 7 - dqy lineouts
python PLOTTERS_PAPER/plot_dqy_3plots.py -v "RL"