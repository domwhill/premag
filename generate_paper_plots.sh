#!/bin/bash -e
# activate conda environment
environment_name=pop_env
if [ "${CONDA_DEFAULT_ENV}" != "${environment_name}" ]; then
  echo "activating conda environment ${environment_name}"
  source activate ${environment_name}
fi

save_dir="./output_images/"
echo "saving images in ${save_dir}"
echo mkdir -p ${save_dir}
mkdir -p ${save_dir}

# Fig. 1 lineouts from IMPACT profiles
python plotters/plot_initial_profile_lineouts.py $save_dir

# Fig 2 f0 
python plotters/plot_f_log.py $save_dir
# Fig 3 a-d
# 3a
python plotters/plot_transport_lineouts_vs_x.py 'SH x' 'Te' $save_dir
# 3b
python plotters/plot_transport_lineouts_vs_x.py 'vN x' 'Te' $save_dir
# 3c
python plotters/plot_transport_lineouts_vs_x.py 'RL y' 'Bz' $save_dir
# 3d
python plotters/plot_transport_lineouts_vs_x.py 'E y' 'Bz' $save_dir

# Fig. 4 - dT along x and y at different point in conduction zone
python plotters/plot_dT_4plot.py $save_dir

# Fig. 5 - Hall-parameter contour plot + Righi-Leduc heating 
python plotters/plot_dyqrl_bier_contour.py -v "RL" -o $save_dir

# Fig. 6 - Hall-parameter contour plot + self-generated magnetic field generation
python plotters/plot_dyqrl_bier_contour.py -v "bier" -o $save_dir

# Fig 7 - dqy lineouts
python plotters/plot_dq_vs_y_3plots.py $save_dir
