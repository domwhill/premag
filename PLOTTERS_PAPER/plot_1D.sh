save_dir='LT1_lineouts'
mkdir ../OUTPUT/$(save_dir)

python2 plot_qRL_abs_1D_POP.py 'SH x' 'Te' $save_dir
python2 plot_qRL_abs_1D_POP.py 'RL y' 'Te' $save_dir
python2 plot_qRL_abs_1D_POP.py 'vN x' 'Te' $save_dir
python2 plot_qRL_abs_1D_POP.py 'E x' 'Bz' $save_dir
python2 plot_qRL_abs_1D_POP.py 'E y' 'Bz' $save_dir


#python2 plot_qRL_abs_1D_POP.py 'SH x' 'Bz' $save_dir
python2 plot_qRL_abs_1D_POP.py 'RL y' 'Bz' $save_dir
#python2 plot_qRL_abs_1D_POP.py 'vN x' 'Bz' $save_dir
#python2 plot_qRL_abs_1D_POP.py 'SH x' 'Te' '$save_dir'
#python2 plot_qRL_abs_1D_POP.py 'SH x' 'Te' '$save_dir'

