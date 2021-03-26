'''
27/08/2019

plot lineouts of dTe (the perturbed electron temperature).



'''
import sys
sys.path.extend(["./"])
import MODULES.figure_latex as fprl
import matplotlib.gridspec as GS
from pylab import *
import MODULES.lineout_plotting as TEL
import MODULES.chfoil_module as cf
import MODULES.plot_utils as utils

#aesthetics
cmap = cm.viridis
# impact y index to plot (for x lineout iy = 0..40, pick iy=20 for as central coord,)
iy = 20
c_list = ['r', 'g', 'b']

#<<<--- transport inputs

fpre = lambda path_in: path_in.split('/')[-1]


def get_save_folder():
    if len(sys.argv) > 1:
        save_folder = sys.argv[1]
        if save_folder[-1] != '/':
            save_folder = save_folder + '/'
    else:
        save_folder = ''
    return save_folder


def get_amp(data):
    return data - np.average(data)


def set_yax(ax, xlab, ylab, xlim):
    ax.grid(c='gray')
    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_xlim(xlim[0], xlim[-1])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax.tick_params(which='both', direction='in')


class PlotdTx:
    """Class to plot a line out dTe vs spatial coordinate along the x direction.

    """

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
        p1dt = TEL.plot_custom_xlineout_amp_tevol(fig,
                                                  ax_t,
                                                  self.path_list,
                                                  var_amp='Te',
                                                  time_list=self.time_list,
                                                  style_list=self.style_list,
                                                  mstyle_list=self.mstyle_list,
                                                  cmap=['gray', 'b', 'k'],
                                                  axleg=ax_t,
                                                  leg_dict=[])

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


class PlotdTy:
    """Plot a lineout of dT (perturbed electron temperature) along the y axis

    """

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


if __name__ == "__main__":
    save_path = get_save_folder()
    # times to plot - 06 is IMPACT time dump number 6
    time_to_plot = '06'
    # magnetic field strengths (T) being plot
    bz_list = [0.0, 50.0, 400.0]


    fig = fprl.newfig_generic_twinx(1.0, scale_width=1.0, scale_ratio=1.4)
    gs = GS.GridSpec(2, 3, height_ratios=[1.0, 0.4])
    gs.update(wspace=0.5, hspace=0.5, left=0.15, right=0.88)
    ax_t = plt.subplot(gs[0, :])
    ax1 = plt.subplot(gs[1, 0])
    ax2 = plt.subplot(gs[1, 1])
    ax3 = plt.subplot(gs[1, 2])

    run_obj_list = utils.RunInfoList('bz',
                                     scale_length=1,
                                     bz_in=bz_list,
                                     pert_amp='1p',
                                     time=time_to_plot)
    # sort by bz - note plot_dT_vs_x expects a list of times)
    PlotdTx(run_obj_list, [time_to_plot]).plot(ax_t)
    PlotdTy(run_obj_list.run_obj_dict['0']).plot(ax1, time_to_plot)
    PlotdTy(run_obj_list.run_obj_dict['50']).plot(ax2, time_to_plot)
    PlotdTy(run_obj_list.run_obj_dict['400']).plot(ax3, time_to_plot)

    ax2.set_ylabel('')
    ax3.set_ylabel('')
    save_name = '%s%s_%s.png' % (save_path, run_obj_list.save_tag, '4plotdT')
    plt.savefig(save_name, dpi=600)
    print('saving as: ', save_name)
    print(' copy and paste: open -a preview ' + save_name)
