'''
    plots the divergence of q_RL
'''
import numpy as np, re, os, sys, getpass, matplotlib.pyplot as plt
import getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')

import MODULES.chfoil_module as cf
from chfoil_module import conv_factors_eos
from chfoil_module import cd5_switches
import plot_comparison as pc
import MODULES.figure_prl_twocol as fprl
import MODULES.kinetic_ohmslaw_module_varZ as kohb
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
#import load_names as LN

userid = getpass.getuser()

#path_tag = cfoil.retrieve_path_tag(path_list[0])

SI_on = cd5_switches.SI_on
save_on = cd5_switches.save_on
hlines_on = cd5_switches.hlines_on
grid_on = cd5_switches.grid_on
horizontal_on = cd5_switches.horizontal_on
separate_plots_on = cd5_switches.separate_plots_on

cwdpath = os.getcwd()
cd5 = cf.conv_factors_custom(cwdpath)
Z0 = cd5.Z
T0 = cd5.T0
n0 = cd5.n0
Bz0 = cd5.Bz_ref
c_index = cd5.c_index
cl_index = cd5.cl_index
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab
xlab_rel = cd5.xlab_rel
ylab = cd5.ylab
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
yi = cd5.yi

divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
lh_style = cd5.lh_style
dashes = cd5.dashes
BZ_FACTOR = Bz0    #dict_B['norm_const'] #Tesla
TIME_FACTOR = tau_ei    #s
BIER_FACTR = (BZ_FACTOR / TIME_FACTOR) * (1e-12)
dtBz_unit = r'[$\si{Tps^{-1}}$]'
MSIZE = 4
#---------
'''
loc_nspace = LN.load_names()
path1 = loc_nspace.single_B
path1_noB = loc_nspace.single_noB
path1_speckle_tc5 = loc_nspace.speckle_B
path1_speckle_tc5_noB = loc_nspace.speckle_noB
save_path = loc_nspace.save_path
'''
save_path = './'


def load_qdata(path, time='10'):
    #-------- LOAD q data ----------------------------------------------
    kohnew = {}

    tt = int(time)
    string_tt_glb = '%2.2i' % int(time)
    #string_tt_glb = '%02i' %  int(time)
    time = string_tt_glb
    print ' ---------- tt ==== ', tt
    print '--------'
    print '\n\nkinetic model time = ', time
    dict = kohb.get_kinetic_heatflow_b(path, str(time))
    kohnew[path] = dict

    return kohnew


def extract_ohms(path, time, x_limit=73):
    '''
    dict = extract_ohms(path,time)
    '''
    #data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time)
    #data_x2d_vN_k,data_y2d_vN_k = data_x2d_vN['kinetic'],data_y2d_vN['kinetic']
    #data_x2d_vN_c,data_y2d_vN_c = data_x2d_vN['classical'],data_y2d_vN['classical']
    ##classic_biermann,kinetic_biermann = kbier.get_kinetic_and_classicB(path,path,time)
    #E_P_classical,E_P_kinetic =kbier.get_kinetic_and_classic_E_pressure(path,fpre(path),time,xlim=73)
    #E_x,E_y = kbier.get_hallfield(grid,rho_vyx,Bz_vyx,jx_vyx,jy_vyx,fo)
    #E_TEwedgex,E_TEwedgey = kbier.get_thermoelectricEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)

    #E_P_c_x, E_P_c_y = E_P_classical['E_P_x'],E_P_classical['E_P_y']
    #E_P_k_x, E_P_k_y = E_P_kinetic['E_P_x'],E_P_kinetic['E_P_y']
    fprefix = fpre(path)
    dict_kinetic = kbier.get_kinetic_E(path, fprefix, time, xlim=73)

    ##data_x2d_vN,data_y2d_vN =
    ##alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time):

    dict = {}
    dict['E_betawedge_x'] = {}
    dict['E_betawedge_x']['data'] = dict_kinetic['E_TEwedge_x']
    dict['E_betawedge_x']['label'] = r'$E_{\beta_{\wedge},x}$'

    dict['E_betawedge_y'] = {}
    dict['E_betawedge_y']['data'] = dict_kinetic['E_TEwedge_y']
    dict['E_betawedge_y']['label'] = r'$E_{\beta_{\wedge},y}$'

    dict['E_hall_x'] = {}
    dict['E_hall_x']['data'] = dict_kinetic['E_hall_x']
    dict['E_hall_x']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},x}$'

    dict['E_hall_y'] = {}
    dict['E_hall_y']['data'] = dict_kinetic['E_hall_y']
    dict['E_hall_y']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},y}$'

    dict['E_P_x'] = {}
    dict['E_P_x']['data'] = dict_kinetic['E_P_x']
    dict['E_P_x']['label'] = r'$E_{\nabla P,x}$'

    dict['E_P_y'] = {}
    dict['E_P_y']['data'] = dict_kinetic['E_P_y']
    dict['E_P_y']['label'] = r'$E_{\nabla P,y}$'

    dict['alpha_jx'] = {}
    dict['alpha_jx']['data'] = dict_kinetic['alpha_jx']
    dict['alpha_jx']['label'] = r'$\alpha_{\perp}\mathbf{j}_x$'

    dict['alpha_jy'] = {}
    dict['alpha_jy']['data'] = dict_kinetic['alpha_jy']
    dict['alpha_jy']['label'] = r'$\alpha_{\perp}\mathbf{j}_y$'

    func_l = lambda data: np.transpose(data)    #[:,:]
    for key in dict.keys():
        dict[key]['data'] = func_l(dict[key]['data'])

    return dict


def fpre(path):
    return path.split('/')[-1]


def plot_dBdt_bier(fig, ax, cax, path, time):

    #kohnew = load_qdata(path,time)
    classic_biermann, kinetic_biermann = kohb.get_kinetic_and_classicB(path, fpre(path), time)
    dict_t = cf.load_dict_1D(path, fpre(path), 'fo', time)
    x_grid_SI = dict_t['x_grid'] * xstep_factor
    y_grid_SI = dict_t['y_grid'] * xstep_factor
    time_col = dict_t['time'] * tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8, 8)    # x,y inches
    lab_dict['lims'] = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[c_index], x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'], lab_dict['ylab'] = xlab, ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'    #'hot'
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #dyqy = np.transpose(dyqy)#q_SH_y[path][tt,:,:]
    data_bier_k = np.transpose(kinetic_biermann)[::-1]
    data_bier_c = np.transpose(classic_biermann)[::-1]

    data_bier_k *= BIER_FACTR
    data_bier_c *= BIER_FACTR
    vmin, vmax = np.min(data_bier_k[:, cl_index:c_index]), np.max(data_bier_k[:, cl_index:c_index])
    print ' np.shape(shape bier k) = ', np.shape(data_bier_k)
    lims_im = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cl_index], x_grid_SI[c_index]]    #
    lims_rev = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[c_index], x_grid_SI[cl_index]]
    lims_im = [x_grid_SI[cl_index], x_grid_SI[c_index], y_grid_SI[0], y_grid_SI[-1]]    #

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(data_bier_k[:, cl_index:c_index],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    norm=norm,
                    aspect=asp,
                    extent=lims_im)
    levels = np.arange(vmin, vmax, (vmax - vmin) / 4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    #for ll in range(len(levels)):
    #    levels[ll] *= divq_factor
    #
    levels = [0.0]

    #ylab_bier = r'$\partial_t \mathbf{B}|_{B}$'
    claby = r'$\partial_t \mathbf{B}|_{B}$ ' + dtBz_unit
    CSy = ax.contour(data_bier_k[::-1, cl_index:c_index],
                     levels=levels,
                     colors=('k'),
                     extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')

    #--- plot hlines
    lim = ax.get_ylim()
    if not lh_style == 'None':
        cf.plot_xvlines(ax, x_grid_SI, lim, lineout_list, linestyle=lh_style, dashes=dashes)

    ax.set_xlabel(xlab_rel)
    ax.set_ylabel(ylab)
    #---
    c2 = fig.colorbar(imy, cax=cax, ax=ax, format='%1i', label=claby)
    c2.add_lines(CSy)

    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()

    return imy


def plot_ylineout_dtBz_custom(fig,
                              axy,
                              path_list,
                              var_amp='Te',
                              time='15',
                              mstyle_list=[None, None, None, 'x', '^', 'o'],
                              style_list=['-', '--', ':', '-', '--', ':'],
                              leg_on=True,
                              dict_list=[]):
    # --- init stuff
    if len(dict_list) != 0:
        color_lineout = dict_list['color_lineout']
        lineout_list = dict_list['lineout_list']

    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data
    take_amp_on = True
    lab_type = 'B'
    n = 4    # number of ks

    #===

    cmap = 'RdBu_r'    #'hot'
    var = 'q_RL'
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #----?

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    print('--- plot_ylineout_qRLy_custom ----')
    print('path_list = ', path_list)
    for pp in range(len(path_list)):
        var = var_amp

        path = path_list[pp]
        fprefix = path_list[pp].split('/')[-1]
        print('--- path = ', path, fprefix)
        ##-------------->>
        #classic_biermann,kinetic_biermann = kohb.get_kinetic_and_classicB(path,fpre(path),time)
        data_bier_c, data_bier_k = kohb.get_kinetic_and_classicB(path, fpre(path), time)

        dict_t = cf.load_dict(path, fpre(path), 'fo', time)
        #x_grid_SI = dict_t['x_grid']*xstep_factor
        #y_grid_SI = dict_t['y_grid']*xstep_factor
        time_col = dict_t['time'] * tstep_factor
        #data_bier_k = np.transpose(kinetic_biermann)[::-1]
        #data_bier_c = np.transpose(classic_biermann)[::-1]

        data_bier_k *= BIER_FACTR
        data_bier_c *= BIER_FACTR
        # Check out the array orientations here..
        T_data = data_bier_k
        T_data_c = data_bier_c

        var = 'dtB'
        tt = int(time)

        #<<<<---------------------
        print '\n --> time = ', time_col, ' ps  = ', dict_t['time'], ' tcol <-----\n'
        x_c_grid = dict_t['x_grid'] * xstep_factor
        y_c_grid = dict_t['y_grid'][1:-1] * xstep_factor

        print('dtb 283 = ', np.shape(data_bier_k), np.shape(data_bier_c))
        print(' shape x y grid = ', np.shape(x_c_grid), np.shape(y_c_grid))
        min_data = np.min(T_data)
        #dict = cf.calc_norms_2(var, min_data)
        norm_const = 1.0
        #c_title = dict['title']
        #c_fmt = dict['c_fmt']
        #var_name = dict['var_name']
        #units = dict['units']

        final_lab = r'$\partial_t \mathbf{B}|_{B}$ ' + dtBz_unit
        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]
        take_amp_on = False
        for xl in range(len(lineout_list)):
            x_lineout = lineout_list[xl]
            ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            if take_amp_on:
                T_dev = T_data[lineout_list[xl], :] * norm_const - np.average(
                    T_data[lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data[lineout_list[xl], :]))

                T_dev_c = T_data_c[lineout_list[xl], :] * norm_const - np.average(
                    T_data_c[lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data_c[lineout_list[xl], :]))

            else:
                T_dev = T_data[lineout_list[xl], :] * norm_const
                T_dev_c = T_data_c[lineout_list[xl], :] * norm_const

            p2, = axy.plot(y_c_grid,
                           T_dev,
                           linestyle=lstyle,
                           c=color_lineout[xl],
                           marker=mstyle,
                           markersize=MSIZE,
                           markevery=3)
            p2c, = axy.plot(y_c_grid,
                            T_dev_c,
                            linestyle='--',
                            c=color_lineout[xl],
                            marker=mstyle,
                            markersize=MSIZE,
                            markevery=3)

            p2.set_linewidth = 2
            leg2_list.append(p2)
            lab_lineout_list.append(ylab_lineout)
        leg_list.append(p2)
        lab_list.append(final_lab)

    if take_amp_on:
        ydifflab = r'$\delta \partial_t B_z$'
    else:
        ydifflab = final_lab
    axy.set_ylabel(ydifflab)
    axy.set_xlabel(ylab)
    if leg_on:
        leg_list = [p2c, p2]
        lab_list = [r'$\delta \partial_t B_c$', r'$\partial_t B_k$']
    axy.grid(color='0.5', linestyle='-')

    return p2


if __name__ == "__main__":
    #-- generate new axis for color bar with gridspec here....
    fig, ax = fprl.newfig(1.0)
    path, time = 'chfoil_default5_heat11_2D1--cx1', '15'
    #matplotlib.rc('font', **{'family':"sans-serif"})

    params = {'text.latex.preamble': [r'\usepackage{siunitx}']}
    plt.rcParams.update(params)

    cax = []
    imy = plot_dyqy_RL(fig, ax, cax, path, time)
    plt.show()
