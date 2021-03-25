'''
Plot of div.qRL or Biermann contours overlayed on the perturbed Hall parameter image plot.

Usage:
# plto biermann
python PLOTTERS_PAPER/plot_dqy_contour.py -v "bier"
# plot RL
python PLOTTERS_PAPER/plot_dqy_contour.py -v "RL"

'''
import argparse
import sys
sys.path.extend(["./"])
from pylab import *
import matplotlib.gridspec as GS
import MODULES.kinetic_ohmslaw_module_varZ as kohb
import MODULES.figure_prl_twocol as fprl
import MODULES.chfoil_module as cf
from MODULES.plot_utils import RunInfoList
import MODULES.gen_dyqy_RL_varZ as qrl
from matplotlib import ticker

#---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
#-----
# functions
fpre = lambda path_in: path_in.split('/')[-1]

def get_command_line_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-v", "--variable",
        default="RL",
        help="Transport term to plot RL x = x component of Righi-Leduc, SH y = y component of Spitzer-HArm/diffusive heat flow,"
        "vN y = y-component of Nernst")
    args = parser.parse_args()
    if args.variable not in ("RL", "bier"):
        raise ValueError("Input arg -v must be either RL or bier")
    return args.variable


def clear_yax(ax):
    fprl.clear_yax(ax)
    ax.set_ylabel('')
    pass


def format_yax(ax):
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%i'))
    pass


def b_lab(bz_in):
    if int(bz_in) == -1:
        lab = 'no B'
    else:
        lab = r'$%i\,\si{T}$' % (bz_in)
    return lab


def custom_im(ax, xgrid, data, lab, lim, **kwargs):

    xmin, xmax, ymin, ymax = lim
    cmap = kwargs.pop('cmap', 'RdBu_r')

    bool_t = (xgrid >= xmin) * (xgrid < xmax)
    #treat_data = lambda a: np.sign(data[bool_t,:])*np.log10(np.abs(a[bool_t,:]))
    treat_data = lambda a: a[bool_t]
    im = ax.imshow(treat_data(data), aspect='auto', extent=[ymin, ymax, xmax, xmin], cmap=cmap)
    return im


def colourbar(fig, im, ax, cax, lab):
    cbar = fig.colorbar(im, ax=ax, cax=cax, label=lab)
    #cbar = fprl.sci_fmt_colorbar(fig, im, ax, cax= cax, lab=lab)
    return cbar


def custom_contour(ax, xgrid, data, lab, lim, **kwargs):

    xmin, xmax, ymin, ymax = lim
    colors = kwargs.pop('colors', 'k')

    bool_t = (xgrid >= xmin) * (xgrid < xmax)
    #treat_data = lambda a: np.sign(data[bool_t,:])*np.log10(np.abs(a[bool_t,:]))
    treat_data = lambda a: a[bool_t]
    std_data = np.std(treat_data(data))
    # set levels
    lvls = kwargs.pop('levels', [-std_data, 0, std_data])

    im = ax.contour(
        treat_data(data),
        levels=lvls,
        colors=(colors),
        extent=[ymin, ymax, xmin, xmax],
    )
    # contour labels
    # ref: https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/contour_label_demo.html#sphx-glr-gallery-images-contours-and-fields-contour-label-demo-py
    strs = [r'-$\sigma$', '0', r'+$\sigma$']
    fmt = {}
    for l, s in zip(im.levels, strs):
        fmt[l] = s
    ax.clabel(im, im.levels, inline=True, fmt=fmt)

    return lvls


def ensure_list(s):
    # Ref: https://stackoverflow.com/a/56641168/
    return s if isinstance(
        s, list) else list(s) if isinstance(s,
                                            (tuple, set, np.ndarray)) else [] if s is None else [s]


def convert_lists_to_set(a_list, b_list, c_list, d_list):
    var_list = set((a, b, c, d) for a in ensure_list(a_list) for b in ensure_list(b_list)
                   for c in ensure_list(c_list) for d in ensure_list(d_list))
    return var_list


class PlotContour():

    def __init__(self, run_obj, **kwargs):
        # --- generating  path_list
        self.run_obj = run_obj
        self.path = run_obj.path
        self.norm = self.run_obj.norm
        self.transport_opt = kwargs.pop('transport_opt', 'bier')
        self.time = kwargs.pop('time', run_obj.time_max)
        self.var = kwargs.pop('var', 'Te')
        #--> lab
        varlab = self.convert_var_to_var_lab(self.var)
        self.lab = '$\\frac{\delta %s}{%s}$' % (varlab, varlab)    # aesthetics
        self.lab = '$\delta %s /\\langle %s \\rangle$' % (varlab, varlab)    # aesthetics
        self.kwargs = {}

    def load(self):
        self.dict_all = cf.load_data_all(self.path, fpre(self.path), self.time)
        self.init_grids(self.dict_all)
        data_var = np.transpose(self.dict_all[self.var])
        if self.var == 'wt':
            self.data = cf.get_U_dev(data_var) * 1e2
            self.lab = r'$10^{-2}\delta \chi$'
            self.kwargs['cmap'] = cm.Greys
        else:
            self.data = cf.get_U_dev_frac(data_var)

        if self.transport_opt == 'RL':
            self.load_RL()
        else:
            self.load_bier()

    def load_bier(self):
        biermann_c, biermann_k = kohb.get_kinetic_and_classicB_c(self.path, fpre(self.path),
                                                                 self.time)
        self.data_c = biermann_c
        self.data_k = biermann_k

    def load_RL(self):

        dyqy, qlab = qrl.get_divqRL(self.path, self.time, switch_kinetic_on=True)
        self.data_k = np.transpose(dyqy * 1.0)

        dyqy_c, qlab_c = qrl.get_divqRL(self.path, self.time, switch_kinetic_on=False)
        self.data_c = np.transpose(dyqy_c * 1.0)

    def run(self, fig, ax):
        self.load()
        self.plot(fig, ax)

    def save(self, save_name):
        plt.savefig(save_name, dpi=600)
        print('saving as: ', save_name)
        print(' copy and paste:\n open -a preview ' + save_name)

    def convert_var_to_var_lab(self, var):
        if len(var) > 1:
            varlab = '%s_%s' % (var[0], var[1])
        else:
            varlab = var
        return varlab

    def plot(self, fig, ax):
        lim = self.get_lim()

        try:
            assert len(self.x_grid) == np.shape(self.data_k)[0]
        except AssertionError:
            x_grid_b = self.x_grid[1:-1]

        self.im = custom_im(ax, self.x_grid, self.data, self.lab, lim, **self.kwargs)
        # contour plots - set contour levels using classical data
        lvls = custom_contour(ax, x_grid_b, self.data_c, self.lab, lim=lim, colors='xkcd:green')
        lvls = custom_contour(ax,
                              x_grid_b,
                              self.data_k,
                              self.lab,
                              lim=lim,
                              colors='xkcd:red',
                              levels=lvls)

        ax.set_xlim(self.y_grid[1], self.y_grid[-2])
        ax.tick_params(which='both', direction='in')
        ax.set_ylabel(self.norm.xlab_rel)
        ax.set_xlabel(r'y [$\si{\micro\meter}$]')

    def init_grids(self, dict_T):
        # ---> functions

        self.x_grid = dict_T['x_grid'][1:-1] * self.norm.xstep_factor
        self.y_grid = dict_T['y_grid'][1:-1] * self.norm.xstep_factor
        if self.run_obj.scale_length == 2:
            self.xmin, self.xmax = 0.0, 30.0
        else:
            self.xmin, self.xmax = 0.0, 20.0
        self.ymin, self.ymax = self.y_grid[0], self.y_grid[-1]

    def get_lim(self):
        lim = [self.xmin, self.xmax, self.ymin, self.ymax]
        return lim


def main():
    #  transport_opt can be 'RL' = Righi-Leduc heat flow divergence or 'bier' = biermann
    transport_opt = get_command_line_args()

    #<<<--- transport inputs
    time = "06"
    #------------------->>>
    bz_list = [50.0, 400.0]

    # Generate run obj list
    obj_list = RunInfoList(var_type='bz', scale_length=1, bz_in=bz_list, lambda_p=5, pert_amp='1p')

    fig = fprl.newfig_generic(1.4, 0.9)
    plt.subplots_adjust(left=0.10, right=0.88, bottom=0.18, top=0.9, wspace=0.1)
    gs1 = GS.GridSpec(1, 3, width_ratios=[1.0, 1.0, 0.08])
    ax1 = fig.add_subplot(gs1[0])
    ax2 = fig.add_subplot(gs1[1])
    ax3 = fig.add_subplot(gs1[2])

    # for obj in run_obj list
    # plot contour
    cont1 = PlotContour(obj_list.run_objs[0],
                        time=time,
                        var='wt',
                        transport_opt=transport_opt)
    cont2 = PlotContour(obj_list.run_objs[1],
                        time=time,
                        var='wt',
                        transport_opt=transport_opt)
    cont1.run(fig, ax1)
    cont2.run(fig, ax2)

    clim1 = cont1.im.get_clim()
    cont1.im.set_clim(clim1[0],
                      clim1[-1] * 1.1)    # make plot slightly more grey by extending limits
    cont2.im.set_clim(cont1.im.get_clim())
    colourbar(fig, cont2.im, ax3, cax=ax3, lab=cont1.lab)
    clear_yax(ax2)
    format_yax(ax1)

    save_name = '%s_%s_bier.png' % (obj_list.save_tag, obj_list.run_objs[0].save_tag)
    fig.savefig(save_name)
    print(' copy and paste:\n open -a preview ' + save_name)

if __name__=="__main__":
    main()
