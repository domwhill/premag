'''
    Module to plot lineouts of variables along y/x directions.

'''

import numpy as np
import re
import os
import sys
import matplotlib.pyplot as plt

sys.path.extend(["./"])
import modules.chfoil_module as cfoil
from modules.chfoil_module import cd5_switches
import modules.house_keeping as hk
import modules.figure_latex as fprl
from pylab import *

SI_on = cd5_switches.SI_on
save_on = cd5_switches.save_on
hlines_on = cd5_switches.hlines_on
grid_on = cd5_switches.grid_on
horizontal_on = cd5_switches.horizontal_on
separate_plots_on = cd5_switches.separate_plots_on
marker_on = False    #cd5_switches.marker_on
MSIZE = 4
init_path = os.getcwd()

norm_path = hk.directory_paths().norm_dir
cd5 = cfoil.ConversionFactors(norm_path)

cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)

SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab
ylab = cd5.ylab
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
yi = cd5.yi
divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
lh_style = cd5.lh_style
dashes = cd5.dashes


def fpre(pathy):
    return pathy.split('/')[-1]


def plot_custom_xlineout_amp_tevol(fig,
                                   ax,
                                   path_list,
                                   var_amp='Te',
                                   time_list=('15'),
                                   vert_lines='False',
                                   cmap=None,
                                   style_list=('-', '--', ':', '-', '--', ':'),
                                   mstyle_list=(None, None, None, 'x', '^', 'o'),
                                   leg_dict=(),
                                   axleg=(),
                                   **kwargs):
    '''
        p2 = plot_custom_xlineout(fig,ax,path_list,var_amp='Te')
    '''
    # --- init stuff

    leg_list = []
    lim_data = 0.01

    # ---- loop over paths
    if cmap is None:
        c_list = cm.plasma(np.linspace(0, 1, len(time_list)))
    elif isinstance(cmap, (list, tuple)):
        c_list = cmap
    else:
        c_list = cmap(np.linspace(0, 1, len(time_list)))

    for pp in range(len(path_list)):
        p_list = []
        leg_list = []
        for tt in range(len(time_list)):
            time = time_list[tt]
            var = var_amp
            fname = cfoil.construct_fname(path_list[pp], fpre(path_list[pp]), var, time)
            dict_T = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), var, time)

            if var == 'wt':
                dict_ne = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), 'n', time)
                ne = dict_ne['mat']
                dict_fo = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), 'fo', time)
                fo = dict_fo['mat']
                nv, ny, nx = np.shape(fo)

                dict_Z = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), 'Z', time)
                prof_Z = dict_Z['mat']
                print('-- loading wt =', np.shape(ne), nx, ny)
                Z2ni = cfoil.trim_array(ne, nx, ny) * cfoil.trim_array(prof_Z, nx, ny)
                hmat = cfoil.trim_array(dict_T['mat'], nx, ny)
                T_data = hmat / Z2ni
            else:
                T_data = dict_T['mat']

            time_col = float(dict_T['time']) * tstep_factor
            x_c_grid = dict_T['x_grid'] * xstep_factor
            min_data = np.min(T_data[cl_index:c_index, :])
            dict = cfoil.calc_norms_2(var, min_data, normal_class=cd5, forced_power=[0])
            norm_const = dict['norm_const']
            c_title = r'$\delta$ ' + dict['title']

            final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
            final_lab = r'$' + final_lab + '$'

            dashes = (4, 2)

            lstyle = style_list[pp]
            mstyle = mstyle_list[pp]
            print 'lstyle = ', lstyle, mstyle
            data_amp_2D = cfoil.get_U_dev_abs(T_data[cl_index:c_index, :])
            data_amp = np.max(data_amp_2D, axis=1)

            p1, = ax.semilogy(x_c_grid[cl_index:c_index],
                              np.abs(data_amp * norm_const),
                              c=c_list[tt],
                              linestyle=lstyle,
                              marker=mstyle,
                              markersize=MSIZE,
                              markevery=3,
                              label=final_lab)

            if lstyle == '--':
                p1.set_dashes(dashes)
            tlab = r'%i\,ps' % time_col
            p_list.append(p1)
            leg_list.append(tlab)
    ax.grid(color='0.5', linestyle='-')
    if len(leg_dict) != 0:
        axleg.legend(p_list, leg_list, bbox_to_anchor=[1.0, 0.5])
    ax.set_xlabel(xlab)
    ax.set_ylabel(c_title)
    return p_list
