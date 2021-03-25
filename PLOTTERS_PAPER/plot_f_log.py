'''
    Plots fo relative to a gaussian with same temperature and density
    for x coord in x_list and y_coords in y_list
'''
import os
import sys
from pylab import *
sys.path.extend(["./"])

from MODULES.plot_utils import RunInfo
from MODULES.plot_utils import GetSimData
import MODULES.figure_prl_twocol as fprl

def get_save_folder():
    if len(sys.argv) > 1:
        save_folder = sys.argv[1]
        if save_folder[-1] != '/':
            save_folder = save_folder + '/'
    else:
        save_folder = ''
    return save_folder


def get_fo_M(v_grid, ne, Te):

    fo_un_norm = np.exp(-(v_grid**2) * (0.5 / Te))
    dv = abs(v_grid[1:] - v_grid[:-1])
    vstep = dv[0]
    ne_norm = np.sum(fo_un_norm * (v_grid**2)) * 4.0 * np.pi * vstep
    fo_final = ne * (ne_norm**-1) * fo_un_norm
    return fo_final


def plot_f0(ax, v_grid, fo, ne, Te):

    fo_M = get_fo_M(v_grid, ne, Te)
    ax.semilogy(v_grid, fo, c='b', linestyle='-')
    ax.semilogy(v_grid, fo_M, c='b', linestyle='--')
    ax.set_xlabel(r'v [$\si{v_{th}}$]')
    ax.set_xlim(0.0, 6.0)


if __name__ == "__main__":
    save_path = get_save_folder()

    lambda_p = 5.0
    sim_obj = RunInfo(scale_length=1, bz_in=50.0, lambda_p=lambda_p, pert_amp='1p', dim='1D')
    sim_data_obj = GetSimData(sim_obj, time='15')
    save_name = '%sfigure2_%s_%i_f0.png' % (save_path, sim_obj.save_tag, int(sim_data_obj.time))

    v_grid = sim_data_obj.v_grid
    fo = sim_data_obj.fo
    Te = sim_data_obj.Te
    ne = sim_data_obj.ne
    x_grid_um = sim_data_obj.x_grid * sim_data_obj.cfg.lambda_mfp * 1e6
    Te_keV = Te * sim_data_obj.Te_ref * 2.0 * 1e-3
    xlab = r'$x - x_{abl}\,[\si{\micro\meter}]$'
    # ---> 3 plots setting
    #fig = fprl.newfig_generic_2yscale(1.4,
    #                                    scale_width=1.2,scale_ratio=0.5)#(1.1,scale_width=1.5,scale_ratio=0.5)#plt.figure()
    #fig.subplots_adjust(left=0.1, right=0.9, wspace = 0.6, top=0.9, bottom=0.28)
    # ---> 2 plots setting
    plot_Te_on = True
    if plot_Te_on:
        fig = fprl.newfig_generic_2yscale(
            1.4, scale_width=1.2,
            scale_ratio=0.5)
        fig.subplots_adjust(left=0.1, right=0.9, wspace=0.4, top=0.9, bottom=0.28)
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
    else:
        fig = fprl.newfig_generic(1.4, 0.9)
        plt.subplots_adjust(left=0.15, right=0.9, bottom=0.18, top=0.9, wspace=0.1)
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)

    ix_list = [100, 170]
    iy = 1
    ix = ix_list[0]
    lim = [1e-12, 10.0]
    plot_f0(ax1, v_grid, fo[:, iy, ix], ne[iy, ix], Te[iy, ix])
    ax1.set_ylabel(r'$f_0$')
    ax1.set_ylim(lim)

    ix = ix_list[1]
    plot_f0(ax2, v_grid, fo[:, iy, ix], ne[iy, ix], Te[iy, ix])
    ax2.set_ylim(lim)
    fprl.clear_yax(ax2)

    if plot_Te_on:
        ax3.plot(x_grid_um, Te_keV[iy, :], c='k')
        for ix in ix_list:
            ax3.scatter(x_grid_um[ix], Te_keV[iy, ix], c='r', marker='x')
        ax3.set_xlim(-5.0, 20.0)
        ax3.set_xlabel(xlab)
        ax3.set_ylabel(r'$T_e \,[\si{keV}]$')
        ax3.tick_params(which='both', direction='in')

    for ax in [ax1, ax2]:
        ax.tick_params(which='both', direction='in')

    print('---> saving as: %s' % (save_name))
    fig.savefig(save_name)
