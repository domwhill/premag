'''
    plots the divergence of q_RL
'''
import numpy as np, re, os, sys, matplotlib.pyplot as plt

sys.path.extend(["./"])
import MODULES.chfoil_module as cf
import MODULES.figure_prl_twocol as fprl
import MODULES.kinetic_ohmslaw_module_varZ as kohb
import MODULES.house_keeping as hk
from matplotlib import ticker

norm_path = hk.directory_paths().norm_dir
cfg = cf.ConversionFactors(norm_path)

#---------
'''
loc_nspace = LN.load_names()
path1 = loc_nspace.single_B
path1_noB = loc_nspace.single_noB
path1_speckle_tc5 = loc_nspace.speckle_B
path1_speckle_tc5_noB = loc_nspace.speckle_noB
save_path = loc_nspace.save_path
'''
save_path = '../PLOTTERS/'


def repack_2D(path, time, cfg=cfg):
    '''
        dict_c,dict_k = repack_2D(path,time)
    '''
    dict_qc, dict_qk = kohb.get_q_abs(path,
                                      time)    #get_q_ratio(path,time)#get_q_individ(path,time)
    v_nx_k, v_ny_k, v_nx_c, v_ny_c = kohb.get_Nernst_abs(path, time)

    dict_out_c = {}
    dict_out_k = {}
    var_list = ['SH x', 'SH y', 'RL x', 'RL y', 'E x', 'E y', 'vN x', 'vN y']

    for var in var_list:
        if var[0] != 'v':
            dict_out_k[var] = dict_qk[var]
        elif var == 'vN x':
            dict_out_k[var] = v_nx_k
        elif var == 'vN y':
            dict_out_k[var] = v_ny_k

    for var in var_list:
        if var[0] != 'v':
            dict_out_c[var] = dict_qc[var]
        elif var == 'vN x':
            dict_out_c[var] = v_nx_c
        elif var == 'vN y':
            dict_out_c[var] = v_ny_c

    dict_out_k['U'] = dict_qk['U']
    dict_out_k['Te'] = dict_qk['Te']
    dict_out_k['Bz'] = dict_qk['Bz']
    dict_out_k['x_grid'] = dict_qk['x_grid']
    return dict_out_c, dict_out_k


def load_qdata(path, time='10', cfg=cfg):
    #-------- LOAD q data ----------------------------------------------
    kohnew = {}

    tt = int(time)
    string_tt_glb = '%2.2i' % int(time)
    time = string_tt_glb
    print '\n\nkinetic model time = ', time
    dict = kohb.get_kinetic_heatflow_b(path, str(time))
    kohnew[path] = dict

    return kohnew


def get_divqRL(path, time, cfg=cfg, **kwargs):
    '''
        dyqy, qlab = get_divqRL(path,time,cfg=cfg,**kwargs)
    :param path:
    :param time:
    :param cfg:
    :param kwargs:
    :return:
    '''
    switch_kinetic_on = kwargs.pop('switch_kinetic_on', True)

    kohnew = load_qdata(path, time)
    dict_c, dict_k = repack_2D(path, time)

    x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
    y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
    # --- do everything in units of 10 eV/ps

    if switch_kinetic_on:
        data_yc = dict_k['RL y']
        kc_str = 'k'
    else:
        data_yc = dict_c['RL y']
        kc_str = 'c'
    data_y = data_yc * 1.0

    dxqy, dyqy = cf.get_grad(x_grid_SI, y_grid_SI, data_y)
    dyqy = np.transpose(dyqy)    # q_SH_y[path][tt,:,:]
    MULT = -1.0

    dyqy *= MULT * cfg.divq_factor * 0.1
    qlab = r'-$\partial_y q_{y,RL,%s}$ %s' % (kc_str, cfg.divq_unit)

    return dyqy, qlab


def extract_ohms(path, time, x_limit=73, cfg=cfg):
    '''
    dict = extract_ohms(path,time)
    '''
    #data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time)
    #data_x2d_vN_k,data_y2d_vN_k = data_x2d_vN['kinetic'],data_y2d_vN['kinetic']
    #data_x2d_vN_c,data_y2d_vN_c = data_x2d_vN['classical'],data_y2d_vN['classical']
    ##classic_biermann,kinetic_biermann = kohb.get_kinetic_and_classicB(path,path,time)
    #E_P_classical,E_P_kinetic =kohb.get_kinetic_and_classic_E_pressure(path,fpre(path),time,xlim=73)
    #E_x,E_y = kohb.get_hallfield(grid,rho_vyx,Bz_vyx,jx_vyx,jy_vyx,fo)
    #E_TEwedgex,E_TEwedgey = kohb.get_thermoelectricEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)

    #E_P_c_x, E_P_c_y = E_P_classical['E_P_x'],E_P_classical['E_P_y']
    #E_P_k_x, E_P_k_y = E_P_kinetic['E_P_x'],E_P_kinetic['E_P_y']
    fprefix = fpre(path)
    dict_kinetic = kbier.get_kinetic_E(path, fprefix, time, xlim=73)

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

    func_l = lambda data: np.transpose(data)
    for key in dict.keys():
        dict[key]['data'] = func_l(dict[key]['data'])

    return dict


def fpre(path):
    return path.split('/')[-1]


def plot_dyqy_RL(fig, ax, cax, path, time, cfg=cfg):

    kohnew = load_qdata(path, time)
    data_y = np.transpose(kohnew[path]["q RL" + '_y']['data'])    #q_SH_y[path][tt,:,:]

    x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
    y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
    time_col = kohnew[path]['time'] * cfg.tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8, 8)    # x,y inches
    lab_dict['lims'] = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.c_index], x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'], lab_dict['ylab'] = cfg.xlab, cfg.ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'    #'hot'
    var = 'q_RL'
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    dxqy, dyqy = cf.get_grad(x_grid_SI, y_grid_SI, data_y)
    dyqy = np.transpose(dyqy)    #q_SH_y[path][tt,:,:]

    dyqy *= cfg.divq_factor * 0.1

    MULT = -1.0
    #vmin,vmax = np.min(MULT*dyqy),np.max(MULT*dyqy)
    #vmin,vmax = -1,2
    multtemp = 0.9
    vmin, vmax = np.min(MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp, np.max(
        MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp
    lims_im = [x_grid_SI[cfg.cl_index], x_grid_SI[cfg.c_index], y_grid_SI[0], y_grid_SI[-1]]    #

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(MULT * dyqy[:, cfg.cl_index:cfg.c_index],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    norm=norm,
                    aspect=asp,
                    extent=lims_im)
    levels = np.arange(vmin, vmax, (vmax - vmin) / 4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        #print 'levels ll = ', ll, levels[ll]
        levels[ll] *= cfg.divq_factor

    #
    levels = [0.0]
    claby = r'-$\partial_y q_{y,RL}$ ' + cfg.divq_unit
    CSy = ax.contour(MULT * dyqy[::-1, cfg.cl_index:cfg.c_index] * cfg.divq_factor,
                     levels=levels,
                     colors=('k'),
                     extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')

    #ax.set_title('dyqy at ' + str(time_col))

    #--- plot hlines
    lim = ax.get_ylim()
    if not cfg.lh_style == 'None':
        cf.plot_xvlines(ax,
                        x_grid_SI,
                        lim,
                        cfg.lineout_list,
                        linestyle=cfg.lh_style,
                        dashes=cfg.dashes)

    ax.set_xlabel(cfg.xlab_rel)
    ax.set_ylabel(cfg.ylab)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #---
    #c2 = fig.colorbar(imy,cax=cax,ax=ax,
    #                    format= '%1i',
    #                    label=claby)
    c2 = fprl.custom_colorbar(fig, imy, ax, cax, claby)

    c2.add_lines(CSy)

    save_dict = {}
    save_dict['x_grid_SI'] = x_grid_SI[cfg.cl_index:cfg.c_index]
    save_dict['y_grid_SI'] = y_grid_SI
    save_dict['lab'] = claby
    save_dict['data'] = dyqy[::-1, cfg.cl_index:cfg.c_index] * cfg.divq_factor
    np.save('test_dyqy.npy', save_dict)
    return imy


def plot_dyqy_RL_c(fig, ax, cax, path, time, cfg=cfg, **kwargs):
    switch_kinetic_on = kwargs.pop('switch_kinetic_on', True)

    kohnew = load_qdata(path, time)
    dict_c, dict_k = repack_2D(path, time)

    x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
    y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
    time_col = kohnew[path]['time'] * cfg.tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8, 8)    # x,y inches
    lab_dict['lims'] = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.c_index], x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'], lab_dict['ylab'] = cfg.xlab, cfg.ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'    #'hot'
    var = 'q_RL'
    print ' shape = ', np.shape(kohnew[path][var + '_x']['data'])
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    data_y = np.transpose(kohnew[path][var + '_y']['data'])
    if switch_kinetic_on:
        data_yc = dict_k['RL y']

    else:
        data_yc = dict_c['RL y']
    data_y = data_yc * 1.0
    dxqy, dyqy = cf.get_grad(x_grid_SI, y_grid_SI, data_y)
    dyqy = np.transpose(dyqy)    #q_SH_y[path][tt,:,:]

    dyqy *= cfg.divq_factor * 0.1

    MULT = -1.0
    multtemp = 0.9
    vmin, vmax = np.min(MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp, np.max(
        MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp

    lims_im = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.cl_index], x_grid_SI[cfg.c_index]]    #
    lims_rev = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.c_index], x_grid_SI[cfg.cl_index]]
    lims_im = [x_grid_SI[cfg.cl_index], x_grid_SI[cfg.c_index], y_grid_SI[0], y_grid_SI[-1]]    #

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(
        MULT * dyqy[:, cfg.cl_index:cfg.c_index],
        cmap=cmap,    #vmin=vmin,vmax=vmax,
        norm=norm,
        aspect=asp,
        extent=lims_im)
    levels = np.arange(vmin, vmax, (vmax - vmin) / 4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        #print 'levels ll = ', ll, levels[ll]
        levels[ll] *= cfg.divq_factor

    #
    levels = [0.0]
    if switch_kinetic_on:
        kc_str = 'k'
    else:
        kc_str = 'c'
    claby = r'-$\partial_y q_{y,RL,%s}$ %s' % (kc_str, cfg.divq_unit)
    CSy = ax.contour(MULT * dyqy[::-1, cfg.cl_index:cfg.c_index] * cfg.divq_factor,
                     levels=levels,
                     colors=('k'),
                     extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')

    #ax.set_title('dyqy at ' + str(time_col))

    #--- plot hlines
    lim = ax.get_ylim()
    if not cfg.lh_style == 'None':
        cf.plot_xvlines(ax,
                        x_grid_SI,
                        lim,
                        cfg.lineout_list,
                        linestyle=cfg.lh_style,
                        dashes=cfg.dashes)

    ax.set_xlabel(cfg.xlab_rel)
    ax.set_ylabel(cfg.ylab)
    #---

    c2 = fprl.custom_colorbar(fig, imy, ax, cax, claby)

    c2.add_lines(CSy)

    return imy


def plot_dBdt_bier(fig, ax, cax, path, time, cfg=cfg):

    #kohnew = load_qdata(path,time)
    kohnew = kbier.get_kinetic_E(path, fpre(path), time, xlim=73)
    x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
    y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
    time_col = kohnew[path]['time'] * cfg.tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8, 8)    # x,y inches
    lab_dict['lims'] = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.c_index], x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'], lab_dict['ylab'] = cfg.xlab, cfg.ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'    #'hot'
    var = 'q_RL'
    print ' shape = ', np.shape(kohnew[path][var + '_x']['data'])
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #data_x = np.transpose(kohnew[path][var +'_x']['data'])#q_SH_x[path][tt,:,:]
    data_y = np.transpose(kohnew[path][var + '_y']['data'])    #q_SH_y[path][tt,:,:]
    dxqy, dyqy = cf.get_grad(x_grid_SI, y_grid_SI, data_y)
    dyqy = np.transpose(dyqy)    #q_SH_y[path][tt,:,:]

    dyqy *= cfg.divq_factor * 0.1

    MULT = -1.0
    #vmin,vmax = np.min(MULT*dyqy),np.max(MULT*dyqy)
    #vmin,vmax = -1,2
    multtemp = 0.9
    vmin, vmax = np.min(MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp, np.max(
        MULT * dyqy[:, cfg.cl_index:cfg.c_index]) * multtemp
    print ' vmin, vmax = ', vmin, vmax
    print ' np.shape(dyqy) = ', np.shape(dyqy)
    #sys.exit()
    #vmin, vmax =  -0.00189426831011, 0.000207009514978
    lims_im = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.cl_index], x_grid_SI[cfg.c_index]]    #
    lims_rev = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[cfg.c_index], x_grid_SI[cfg.cl_index]]
    lims_im = [x_grid_SI[cfg.cl_index], x_grid_SI[cfg.c_index], y_grid_SI[0], y_grid_SI[-1]]    #

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(MULT * dyqy[:, cfg.cl_index:cfg.c_index],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    norm=norm,
                    aspect=asp,
                    extent=lims_im)
    levels = np.arange(vmin, vmax, (vmax - vmin) / 4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        levels[ll] *= cfg.divq_factor

    #
    levels = [0.0]
    claby = r'-$\partial_y q_{y,RL}$ '    #+ cfg.divq_unit
    CSy = ax.contour(MULT * dyqy[::-1, cfg.cl_index:cfg.c_index] * cfg.divq_factor,
                     levels=levels,
                     colors=('k'),
                     extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')

    #ax.set_title('dyqy at ' + str(time_col))

    #--- plot hlines
    lim = ax.get_ylim()
    if not cfg.lh_style == 'None':
        cf.plot_xvlines(ax,
                        x_grid_SI,
                        lim,
                        cfg.lineout_list,
                        linestyle=cfg.lh_style,
                        dashes=cfg.dashes)

    ax.set_xlabel(cfg.xlab_rel)
    ax.set_ylabel(cfg.ylab)
    c2 = fig.colorbar(imy, cax=cax, ax=ax, format='%1i', label=claby)
    c2.add_lines(CSy)

    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()
    return imy


def plot_ylineout_qRLy_custom(fig,
                              axy,
                              path_list,
                              var_amp='Te',
                              time='15',
                              mstyle_list=[None, None, None, 'x', '^', 'o'],
                              style_list=['-', '--', ':', '-', '--', ':'],
                              leg_on=True,
                              dict_list=[],
                              cfg=cfg):
    # --- init stuff
    if len(dict_list) != 0:
        cfg.color_lineout = dict_list['color_lineout']
        cfg.lineout_list = dict_list['lineout_list']

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
    save_name = '2D1_' + var_amp + 'amp'
    print('--- plot_ylineout_qRLy_custom ----')
    print('path_list = ', path_list)
    for pp in range(len(path_list)):
        var = var_amp

        path = path_list[pp]
        fprefix = path_list[pp].split('/')[-1]
        print('--- path = ', path, fprefix)
        ##-------------->>
        kohnew = load_qdata(path, time)
        x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
        y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
        time_col = kohnew[path]['time'] * cfg.tstep_factor
        asp = 'auto'

        dict_c, dict_k = repack_2D(path, time)
        T_data_c = dict_c['RL y']

        cmap = 'RdBu_r'    #'hot'
        var = 'q_RL'
        tt = int(time)

        #--- do everything in units of 10 eV/ps
        T_data = np.transpose(kohnew[path][var + '_y']['data'])    #q_SH_y[path][tt,:,:]
        #<<<<---------------------

        MULT = -1.0
        multtemp = 0.9

        dict_T = cf.load_dict(path, fprefix, 'Te', time)
        time_col = float(dict_T['time']) * cfg.tstep_factor
        print '\n --> time = ', time_col, ' ps  = ', dict_T['time'], ' tcol <-----\n'
        x_c_grid = dict_T['x_grid'] * cfg.xstep_factor
        y_c_grid = dict_T['y_grid'][1:-1] * cfg.xstep_factor

        min_data = np.min(T_data)
        dict = cf.calc_norms_2(var, min_data)
        norm_const = 1.0
        c_title = dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab = r'$q_{RL,y}$'
        var_name = r'$q_{RL}$'
        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]
        take_amp_on = True
        print(' shape qrL = ', np.shape(y_c_grid))
        print(' shape qrL = ', np.shape(T_data))
        for xl in range(len(cfg.lineout_list)):
            x_lineout = cfg.lineout_list[xl]
            cfg.ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            if take_amp_on:
                T_dev = T_data[cfg.lineout_list[xl], :] * norm_const - np.average(
                    T_data[cfg.lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data[cfg.lineout_list[xl], :]))

                T_dev_c = T_data_c[cfg.lineout_list[xl], :] * norm_const - np.average(
                    T_data_c[cfg.lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data_c[cfg.lineout_list[xl], :]))

            else:
                T_dev = T_data[cfg.lineout_list[xl], :] * norm_const
                T_dev_c = T_data_c[cfg.lineout_list[xl], :] * norm_const

            p2, = axy.plot(y_c_grid,
                           T_dev,
                           linestyle=lstyle,
                           c=cfg.color_lineout[xl],
                           marker=mstyle,
                           markersize=cfg.MSIZE,
                           markevery=3)
            p2c, = axy.plot(y_c_grid,
                            T_dev_c,
                            linestyle='--',
                            c=cfg.color_lineout[xl],
                            marker=mstyle,
                            markersize=cfg.MSIZE,
                            markevery=3)

            p2.set_linewidth = 2
            leg2_list.append(p2)
            lab_lineout_list.append(cfg.ylab_lineout)
        leg_list.append(p2)
        lab_list.append(final_lab)

    if take_amp_on:
        ydifflab = var_name + r'$ - \langle $' + var_name + r'$ \rangle$ ' + units
        ydifflab = r'$\delta q_{y,RL}$'
    else:
        ydifflab = final_lab
    axy.set_ylabel(ydifflab)
    axy.set_xlabel(cfg.ylab)
    if leg_on:
        leg_list = [p2c, p2]
        lab_list = [r'$\delta q_{y,RL,c}$', r'$\delta q_{y,RL,k}$']
        #axy.legend(leg_list,lab_list,numpoints=1)
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
