'''Figure 1 of paper.

    This file plots the inital profiles then the evolved profiles for a 50T
    sim. in format for 2column paper

    Usage instructions

        python PLOTTERS_PAPER/plot_initial_profile_lineouts.py output_folder_to_save_plot/
'''
import os
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl

sys.path.extend(["./"])

import MODULES.house_keeping as hk
import MODULES.chfoil_module as cf
import MODULES.figure_latex as fprl

paths = hk.directory_paths()
norm_path = paths.norm_dir
cfg = cf.ConversionFactors(norm_path)

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
#---------------
# user parameters
mult = 6
cfg.yi = 1

dim = '1D'
bz_in = 50
# scale length (1,2,3) - (1 = 20um conduction zone sim,1.5, 2 = 40um conduction zone sim etc. )
s_len = 1
path1 = paths.get_path(s_len, 50, lambda_p=5, dim=dim)
path2 = paths.get_path(s_len, 400, lambda_p=5, dim=dim)
# List of paths to simulations to plot
path_list = [path1]
# plot limits
xmin, xmax = -10.0, 25.0


def get_save_name():
    save_suffix = 'one_dim_2col' + time
    if len(sys.argv) > 1:
        save_folder = sys.argv[1]
        if save_folder[-1] != '/':
            save_folder = save_folder + '/'
    else:
        save_folder = ''
    savename = save_folder + 'figure1_initial_profiles_' + save_suffix + '.png'
    return savename


def plot_ax(ax1, x_grid, data, norm_const, c='b', tlab='00', cfg=cfg):
    if len(np.shape(data)) > 1:
        if len(data[0, :]) < 3:
            x_gr = x_grid[:len(data[:, 0])]
            ax1.plot(x_gr, data[:, 0] * norm_const, c=c, label=r'\textit{$' + str(tlab) + '$}')
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        else:
            x_gr = x_grid[:len(data[:, 0])]
            ax1.plot(x_gr, data[:, cfg.yi] * norm_const, c=c, label=r'\textit{$' + str(tlab) + '$}')
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
    else:
        x_gr = x_grid[:len(data)]
        ax1.plot(x_gr, data[:, cfg.yi] * norm_const, c=c, label=r'\textit{$' + str(tlab) + '$}')

    return


if __name__ == "__main__":
    # variables to plot lineouts of.
    var_list = ['Cx', 'n', 'Te', 'Bz', 'wt']

    fig = fprl.newfig_generic(1.4, 0.9)
    plt.subplots_adjust(left=0.1, right=0.92, bottom=0.18, top=0.9, wspace=0.1)

    axn0 = fig.add_subplot(121)
    ax20 = axn0.twinx()
    axn1 = fig.add_subplot(122)
    ax21 = axn1.twinx()
    ax_list = [(axn0, ax20), (axn1, ax21)]

    lstyle_list = ['-', '--', ':']
    var_extra = ''
    var_extra_lab = r''
    var_extra_ylab = r''
    # simulation times to plot
    tt_list = [0, 15]
    for path in path_list:
        for itt, tt in enumerate(tt_list):
            axn = ax_list[itt][0]
            ax2 = ax_list[itt][1]
            lstyle = lstyle_list[path_list.index(path)]

            #--- load data
            time = '%02i' % tt
            fprefix = path.split('/')[-1]
            Te_dict = cf.load_dict(path, fprefix, 'Te', time)
            n_dict = cf.load_dict(path, fprefix, 'n', time)
            Cx_dict = cf.load_dict(path, fprefix, 'Cx', time)
            Bz_dict = cf.load_dict(path, fprefix, 'Bz', time)

            T_data = Te_dict['mat']
            n_data = n_dict['mat']
            Cx_data = Cx_dict['mat']
            Bz_data = Bz_dict['mat']
            time_col = Te_dict['time']

            x_grid_ccg = n_dict['x_grid'] * cfg.xstep_factor
            power = 0
            mod = ''
            # rescaling units to make lineouts all fit easily on one plot
            n_mult, T_mult, C_mult, B_mult = 1.0 / 40.0, 2.5, 20.0, 1.2

            ne_lab = r'$n_e$[\SI{%1.1e}{\centi\meter^{-3}}]' % (cfg.n0 * (n_mult**-1))
            Te_lab = r'$T_e$[\SI{%1.1f}{\kilo\electronvolt}]' % ((T_mult**-1) * 2.0 * cfg.T0 * 1e-3)
            Cx_lab = r'$C_x$[\SI{%1.1e}{\kilo\meter\per\second}]' % ((cfg.v_te * 1e-3) *
                                                                     (C_mult**-1))
            Bz_lab = r'$B_z$ [\SI{%1.1e}{T}]' % ((cfg.Bz0**-1) * (B_mult**-1))

            ne_norm = n_data[cfg.cl_index:cfg.c_index, cfg.yi] * n_mult
            Te_norm = T_data[cfg.cl_index:cfg.c_index, cfg.yi] * T_mult
            C_norm = Cx_data[cfg.cl_index:cfg.c_index, cfg.yi] * C_mult
            B_norm = Bz_data[cfg.cl_index:cfg.c_index, cfg.yi] * B_mult

            b1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],
                           B_norm,
                           'k',
                           linestyle=lstyle,
                           label='$B_z$')
            n1, = ax2.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],
                           ne_norm,
                           c='g',
                           linestyle=lstyle,
                           label=r'$n_e$')
            t1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],
                           Te_norm,
                           c='r',
                           linestyle=lstyle,
                           label='$T_e$')
            Cx1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],
                            C_norm,
                            c='b',
                            linestyle=lstyle,
                            label='$C_x$')

            if var_extra != '':
                wt_dict = cf.load_dict(path, fprefix, var_extra, time)
                wt_data = wt_dict['mat']
                wt_lab = var_extra_lab
                wt_norm = wt_data[cfg.cl_index:cfg.c_index, cfg.yi]
                wt1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],
                                wt_norm,
                                c='k',
                                linestyle=lstyle,
                                label='$\omega \tau$')
            ylim = axn.get_ylim()
            axn.set_ylim(ylim[0], ylim[1] * 1.2)
            axn.set_xlim(xmin, xmax)
            time_si = r'%i $\si{ps}$' % (time_col * cfg.tstep_factor)
            if itt == 1:
                dx_mult = 0.38
                dy_mult = 0.28
            else:
                dx_mult = 0.3
                dy_mult = 0.24
            cf.annotate_time(axn,
                             lett=time_si,
                             dx_mult=dx_mult,
                             dy_mult=dy_mult,
                             loc='top',
                             fontsize=8)

        # axn.axhline(n_crit_norm_oneplot,c='k',linestyle='--')
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)

    yn1, yn2 = axn1.get_ylim()
    ya1, ya2 = ax21.get_ylim()

    axn0.set_ylim(yn1, yn2)
    ax20.set_ylim(ya1, ya2)

    ax_all = [axn0, ax20, axn1, ax21]
    map(lambda ax: ax.tick_params(direction='in'), ax_all)
    if var_extra != '':
        p_list = [n1, t1, Cx1, b1, wt1]
        leg_lab_list = [ne_lab, Te_lab, Cx_lab, Bz_lab, wt_lab]
    else:
        p_list = [n1, t1, Cx1, b1]
        leg_lab_list = [ne_lab, Te_lab, Cx_lab, Bz_lab]

    leg = axn0.legend(p_list,
                      leg_lab_list,
                      ncol=1,
                      fontsize=7,
                      loc='upper left',
                      framealpha=1.0,
                      frameon=False)

    axn0.set_xlabel(r'$x-x_{abl}$ [$\si{\micro\meter}$ ]')
    axn1.set_xlabel(r'$x-x_{abl}$ [$\si{\micro\meter}$ ]')
    axn0.set_ylabel(r'$T_e,\, C_x,\, B_z' + var_extra_ylab + '$')
    axn1.set_yticks([])
    ax20.set_yticks([])

    ax21.set_ylabel(r'$n_e$')
    axn0.set_xlim(xmin, xmax)
    axn1.set_xlim(xmin, xmax)

    savename = get_save_name()
    print 'saving fig as ', savename
    print 'copy to terminal: \n open -a preview ' + savename

    plt.savefig(savename, dpi=400)

    axn.cla()
    ax2.cla()
    plt.close()