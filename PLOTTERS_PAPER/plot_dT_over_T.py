'''
27/08/2019
- Rehash of plot_q_lineout.py


'''
import numpy as np, sys, os, getpass, site, re
sys.path.extend(["./"])
import matplotlib.pyplot as plt
import MODULES.figure_prl_twocol as fprl
import matplotlib.gridspec as GS
import matplotlib.ticker as ticker
from pylab import *
import MODULES.kinetic_ohmslaw_module_varZ as q_mod
import PLOTTERS.gen_Te_lineout as TEL
import gen_dyqy_RL_varZ as qrl
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
import MODULES.tsi_module as tsi

#---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
#-----
# functions
fpre = lambda path_in: path_in.split('/')[-1]


def b_lab(bz_in):
    if int(bz_in) == -1:
        lab = 'no B'
    else:
        lab = r'$%i\,\si{T}$' % (bz_in)
    return lab


#---> file inputs
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir
norm_dir = paths.norm_dir
log_file = norm_dir + 'norm.log'
[T_ref, n_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir, Z_ref, Ar=6.51)

Te_factor = 2.0 * T_ref * 1e-3

path_list = []

#--- generating  path_list
scale_len = 1
lambda_p = 5
#bz_list = [-1,50.0,100.0]
lab_list = []
bz_in = 400.0
pert_amp = '100p'
path = paths.get_path(scale_len, bz_in, lambda_p, pert_amp=pert_amp)

#--->  getting path ->
t_list, tc_list = cf.get_t_list(path, var='Te')
time = '06'    #t_list[-1]
save_tag = 'LT%i_%i_%s' % (scale_len, bz_in, time)

save_name = '%sdU_im_%s.png' % (save_path, save_tag)

#aesthetics
cmap = cm.viridis
style_list = [':', '--', '-', '-', '--', ':']
mstyle_list = [None, None, None, 'x', '^', 'o']

#---> functions
dict_T = cf.load_data_all(path, fpre(path), time)

x_grid = dict_T['x_grid'][1:-1] * cd5.xstep_factor
y_grid = dict_T['y_grid'] * cd5.xstep_factor
if scale_len == 2:
    xmin, xmax = 0.0, 25.0
else:
    xmin, xmax = 0.0, 20.0
ymin, ymax = y_grid[0], y_grid[-1]



def custom_im(ax, xgrid, data, lab):
    cmap = kwargs.pop('cmap', cmap)
    bool_t = (xgrid >= xmin) * (xgrid < xmax)
    #treat_data = lambda a: np.sign(data[bool_t,:])*np.log10(np.abs(a[bool_t,:]))
    treat_data = lambda a: a[bool_t]
    im = ax.imshow(treat_data(data), aspect='auto', extent=[ymin, ymax, xmax, xmin], cmap=cmap)
    plt.colorbar(im, ax=ax, label=lab)


def custom_contour(ax, xgrid, data, lab, **kwargs):
    colors = kwargs.pop('colors', 'k')
    bool_t = (xgrid >= xmin) * (xgrid < xmax)
    #treat_data = lambda a: np.sign(data[bool_t,:])*np.log10(np.abs(a[bool_t,:]))
    treat_data = lambda a: a[bool_t]
    std_data = np.std(treat_data(data))

    im = ax.contour(
        treat_data(data),
        levels=[-std_data, 0, std_data],
        colors=(colors),
        extent=[ymin, ymax, xmin, xmax],
    )
    #ax.set_ylim(xmax,xmin)




#<<<--- transport inputs

# data in
T_data = np.transpose(dict_T['Te'])
Bz_data = np.transpose(dict_T['Bz'])
data = cf.get_U_dev_frac(T_data)

dyqy, qlab = qrl.get_divqRL(path, time)
data_k = np.transpose(dyqy * 1.0)

dyqy_c, qlab_c = qrl.get_divqRL(path, time, switch_kinetic_on=False)
data_c = np.transpose(dyqy_c * 1.0)

fig = fprl.newfig_generic_twinx(1.0, scale_width=1.0, scale_ratio=1.0)
ax = fig.add_subplot(111)
lab = '$\\frac{\delta T_e}{T_e}$'
custom_im(ax, x_grid, data, lab)
custom_contour(ax, x_grid, data_c, lab, colors='b')
custom_contour(ax, x_grid, data_k, lab)

#------------------->>>
#plt.show()

plt.savefig(save_name, dpi=600)
print('saving as: ', save_name)
print(' copy and paste:\n open -a preview ' + save_name)
os.system('open -a preview ' + save_name)
#plt.show()
#plt.close()
