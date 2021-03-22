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
import pdb
from pylab import *
import MODULES.kinetic_ohmslaw_module_varZ as q_mod
import PLOTTERS.gen_Te_lineout as TEL
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
import MODULES.plot_utils as utils
import MODULES.tsi_module as tsi
'''

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
[T_ref,n_ref,Z_ref,Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir,Z_ref,Ar=6.51)


Te_factor = 2.0*T_ref*1e-3
#--->  getting path ->
lambda_p = 5

path_list = []
#--- generating  path_list
scale_len = 1
save_tag = 'LT%i_' % (scale_len)
bz_list = [-1,50.0,100.0]
lab_list = []
for ip in range(len(bz_list)):    
    path_loc = paths.get_path(scale_len,bz_list[ip],lambda_p)
    path_list.append(path_loc)
    lab_list.append(b_lab(bz_list[ip]))
'''
#aesthetics
cmap = cm.viridis

#<<<--- transport inputs
time = '07'
iy = 20
c_list = ['r', 'g', 'b']

#-------------------------------------->

fpre = lambda path_in: path_in.split('/')[-1]


def get_amp(data):
    return data - np.average(data)


def set_yax(ax, xlab, ylab, xlim):
    ax.grid(c='gray')
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_xlim(xlim[0], xlim[-1])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.tick_params(which='both', direction='in')


class plot_dT_vs_x:

    def __init__(self, run_obj_list, time_list):
        self.style_list = [':', '--', '-', '-', '--', ':']
        self.mstyle_list = [None, None, None, 'x', '^', 'o']
        self.run_obj_list = run_obj_list
        self.path_list = self.run_obj_list.paths
        self.norm = self.run_obj_list.run_objs[0].norm
        self.time_list = time_list
        self.xmin, self.xmax = -2.0, 20.0

    def plot(self, ax_t):
        ax_Te = ax_t.twinx()
        p1dt = TEL.plot_custom_xlineout_amp_tevol(
            fig,
            ax_t,
            self.path_list,
            var_amp='Te',
            time_list=self.time_list,
            style_list=self.style_list,
            mstyle_list=self.mstyle_list,
            cmap=['gray', 'b', 'k'],
            axleg=ax_t,
            leg_dict=[])    #{'title': self.norm.tlab, 'some_other thing': 1.0})
        # leg = ax_t.legend(p1dt, leg_list,loc = 'bottom right')

        var = 'Te'
        p_list = []
        time = self.time_list[-1]
        bz_list = []
        for pp, run_obj in enumerate(self.run_obj_list.run_objs):
            path = run_obj.path
            dict_T = cf.load_dict(path, fpre(path), var, time)
            data_grid = dict_T['x_grid'] * self.norm.xstep_factor
            data_plot = dict_T['mat'][:, iy] * (2.0 * run_obj.norm.T0 / 1e3)
            p2, = ax_Te.plot(data_grid, data_plot, c='k', linestyle=self.style_list[pp])
            p_list.append(p2)
            bz_list.append(run_obj.bz)

            # plot crosses
            self.plot_crosses(ax_Te, run_obj, data_grid, data_plot)

        ax_Te.set_ylabel(r'$T_e\,[\si{keV}]$')
        # legend
        id_sort = np.argsort(bz_list)
        leg2 = ax_Te.legend(list(np.array(p_list)[id_sort]),
                            list(np.array(self.run_obj_list.tags)[id_sort]))
        leg2.get_frame().set_linewidth(0.0)
        leg2.get_frame().set_facecolor('none')

        ax_t.set_xlim(self.xmin, self.xmax)
        ax_Te.tick_params(axis='both', which='both', direction='in')
        ax_t.tick_params(axis='both', which='both', direction='in')

    def plot_crosses(self, ax, run_obj, x_grid, data):
        idx_list = run_obj.get_idx_list()
        for ix, idx in enumerate(idx_list):
            ax.scatter(x_grid[idx], data[idx], c=c_list[ix], marker='x')

        pass


class plot_dT_vs_y:

    def __init__(self, run_obj, **kwargs):
        # assuming norms are the same for all run objs which is the case for all premag runs.
        self.run_obj = run_obj
        self.norm = self.run_obj.norm
        self.time = kwargs.get('time', '06')

    def plot(self, ax, time):
        # --->
        self.time = time
        # plot one simulation per plot
        p_list = []
        lab_list = []

        p_list = []
        lab_list = []

        # --> load data
        path = self.run_obj.path
        leg_tag = self.run_obj.tag
        idx_list = self.run_obj.get_idx_list()
        #-->--->
        dict_all = cf.load_data_all(path, cf.fpre(path), self.time)
        data = dict_all['Te'].T
        y_grid = dict_all['y_grid'][1:-1]

        for idx, ix in enumerate(idx_list):
            color = c_list[idx]
            p1, = ax.plot(y_grid * self.norm.xstep_factor,
                          get_amp(data[ix, :]),
                          c=color,
                          linestyle='-')

            p_list.append(p1)
            lab_list.append(leg_tag)

            xlab = r'y [$\si{\micro\meter}$]'
            ylab = r'$\delta T_e$ [\si{keV}]'
            ymin, ymax = y_grid[1] * self.norm.xstep_factor, y_grid[-2] * self.norm.xstep_factor

            set_yax(ax, xlab, ylab, [ymin, ymax])


#-------------------------------------->

fig = fprl.newfig_generic_twinx(1.0, scale_width=1.0, scale_ratio=1.4)
gs = GS.GridSpec(2, 3, height_ratios=[1.0, 0.4])
gs.update(wspace=0.5, hspace=0.5, left=0.15, right=0.88)
ax_t = plt.subplot(gs[0, :])
ax1 = plt.subplot(gs[1, 0])
ax2 = plt.subplot(gs[1, 1])
ax3 = plt.subplot(gs[1, 2])

time_list = ['06']
bz_list = [0.0, 50.0, 400.0]
run_obj_list = utils.run_obj_list('bz', scale_length=1, bz_in=bz_list, pert_amp='1p')
time_max = run_obj_list.tmax_index
# sort by bz

plot_dT_vs_x(run_obj_list, time_list).plot(ax_t)
plot_dT_vs_y(run_obj_list.run_obj_dict['0']).plot(ax1, time_max)
plot_dT_vs_y(run_obj_list.run_obj_dict['50']).plot(ax2, time_max)
plot_dT_vs_y(run_obj_list.run_obj_dict['400']).plot(ax3, time_max)

ax2.set_ylabel('')
ax3.set_ylabel('')
save_name = '%s_%s.png' % (run_obj_list.save_tag, '4plotdT')
#------------------->>>

#====================================================

#=======
#plt.show()

plt.savefig(save_name, dpi=600)
print('saving as: ', save_name)
print(' copy and paste: open -a preview ' + save_name)

#plt.show()
#plt.close()
