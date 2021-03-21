'''

    Te line out --- 
'''

import numpy as np, re, os, sys, getpass, matplotlib.pyplot as plt, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import MODULES.chfoil_module as cfoil
from chfoil_module import conv_factors_eos
from chfoil_module import cd5_switches
import MODULES.house_keeping as hk
import figure_prl_twocol as fprl
import matplotlib.ticker as ticker
from pylab import *
userid = getpass.getuser()

#path_tag = cfoil.retrieve_path_tag(path_list[0])

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
cd5 = cfoil.conv_factors_custom(norm_path)

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


#-----------------------------------------------------------------------
def plot_Te_xlineout(fig, ax, path_list, var_amp='Te', time='15', vert_lines='False', dict_list=[]):
    '''
        plot_Te_xlineout(fig,ax,path_list,var_amp='Te')
    '''

    if len(dict_list) != 0:
        color_lineout = dict_list['color_lineout']
        lineout_list = dict_list['lineout_list']

    #SI_on = True
    #norm_name = '/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/chfoil_default5_norm.txt'
    #cd5 = cfoil.conv_factors_cd5(SI_on,norm_name)
    # --- init stuff
    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data
    #yi = 20
    #xi = 60 # arbitrary x index to plot at
    #c_index = 200

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    save_name = '2D1_' + var_amp + 'amp'
    style_list = ['-', '--', ':', '-', '--', ':']
    mstyle_list = [None, None, None, 'x', '^', 'o']
    for pp in range(len(path_list)):
        #print ' path = ', path_list[pp]
        var = var_amp

        fname = cfoil.construct_fname(path_list[pp], fpre(path_list[pp]), var, time)
        dict_T = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), var, time)
        T_data = dict_T['mat']
        time_col = float(dict_T['time']) * tstep_factor
        x_c_grid = dict_T['x_grid'] * xstep_factor
        #x_c_grid = x_c_grid[:len(T_data[cl_index:c_index,0])]
        y_c_grid = dict_T['y_grid'] * xstep_factor
        print ' x_c_grid shape = ', np.shape(x_c_grid), np.shape(T_data[cl_index:c_index])

        ax.yaxis.set_ticks(np.array([0.6, 0.8, 1.0]))
        #ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))

        min_data = np.min(T_data[cl_index:c_index, :])
        dict = cfoil.calc_norms_2('Te keV', min_data)
        norm_const = dict['norm_const']
        c_title = dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
        final_lab = r'$' + final_lab + '$'
        lstyle, marker = cfoil.get_path_style(final_lab, h_lab, B_lab)

        b00 = re.search('=0$', B_lab)
        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]

        dashes = (4, 2)

        print 'lstyle = ', lstyle, mstyle
        #if h_lab == 'static':
        #   lstyle=':'

        if ~marker_on:
            mstyle = None
        #lstyle='-'
        yminax, ymaxax = ax.get_ylim()

        p1, = ax.plot(x_c_grid[cl_index:c_index],
                      T_data[cl_index:c_index, yi] * norm_const,
                      c='k',
                      linestyle=lstyle,
                      marker=mstyle,
                      markersize=MSIZE,
                      markevery=3,
                      label=final_lab)
        leg_list_x.append(p1)
        lab_list_x.append(final_lab)
        save_name = save_name + str(h_lab)
        if vert_lines:
            cfoil.plot_xvlines(ax,
                               x_c_grid, [],
                               lineout_list,
                               custom_cmap=color_lineout,
                               linestyle=lh_style,
                               dashes=dashes)
            #self.color_lineout
        '''
        for xl in range(len(lineout_list)):
            x_lineout = lineout_list[xl]
            #ylab_lineout = r'$%3.1f$' % x_grid_im[x_lineout]
            ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            
            if hlines_on:
                #ax.axvline(x_grid_im[x_lineout],yminax,ymaxax,
                ax.axvline(x_c_grid[x_lineout],yminax,ymaxax,
                            c=color_lineout[xl],
                            linewidth=0.9,
                            linestyle=lh_style,
                            dashes=dashes)
        '''

    ax.grid(color='0.5', linestyle='-')

    ax.set_xlabel(xlab)
    ax.set_ylabel(c_title)
    return p1


#plot_custom_xlineout_amp
#-----------------------------------------------------------------------
def plot_custom_xlineout_amp(fig,
                             ax,
                             path_list,
                             var_amp='Te',
                             time='15',
                             vert_lines='False',
                             style_list=['-', '--', ':', '-', '--', ':'],
                             mstyle_list=[None, None, None, 'x', '^', 'o'],
                             leg_dict=[]):
    '''
        p2 = plot_custom_xlineout(fig,ax,path_list,var_amp='Te')
    '''
    #SI_on = True
    #norm_name = '/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/chfoil_default5_norm.txt'
    #cd5 = cfoil.conv_factors_cd5(SI_on,norm_name)
    # --- init stuff
    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data
    #yi = 20
    #xi = 60 # arbitrary x index to plot at
    #c_index = 200

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    p_list = []
    save_name = '2D1_' + var_amp + 'amp'
    #style_list = ['-','--',':','-','--',':']
    #mstyle_list = [None,None,None,'x','^','o']
    for pp in range(len(path_list)):
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
        y_c_grid = dict_T['y_grid'] * xstep_factor
        min_data = np.min(T_data[cl_index:c_index, :])
        dict = cfoil.calc_norms_2(var, min_data, normal_class=cd5, forced_power=[0])
        norm_const = dict['norm_const']
        c_title = r'$\delta$ ' + dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
        final_lab = r'$' + final_lab + '$'
        #lstyle,marker = cfoil.get_path_style(final_lab,h_lab,B_lab)

        dashes = (4, 2)

        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]
        print 'lstyle = ', lstyle, mstyle
        yminax, ymaxax = ax.get_ylim()
        data_amp_2D = cfoil.get_U_dev_abs(T_data[cl_index:c_index, :])
        data_amp = np.max(data_amp_2D, axis=1)

        p1, = ax.semilogy(x_c_grid[cl_index:c_index],
                          np.abs(data_amp * norm_const),
                          c='k',
                          linestyle=lstyle,
                          marker=mstyle,
                          markersize=MSIZE,
                          markevery=3,
                          label=final_lab)

        if lstyle == '--':
            p1.set_dashes(dashes)
        p_list.append(p1)
    ax.grid(color='0.5', linestyle='-')
    if len(leg_dict) != 0:
        ax.legend(p_list, leg_dict)
    ax.set_xlabel(xlab)
    ax.set_ylabel(c_title)
    return p1


#-----------------------------------------------------------------------


def plot_custom_xlineout_amp_tevol(fig,
                                   ax,
                                   path_list,
                                   var_amp='Te',
                                   time_list=['15'],
                                   vert_lines='False',
                                   style_list=['-', '--', ':', '-', '--', ':'],
                                   mstyle_list=[None, None, None, 'x', '^', 'o'],
                                   leg_dict=[],
                                   axleg=[]):
    '''
        p2 = plot_custom_xlineout(fig,ax,path_list,var_amp='Te')
    '''
    # --- init stuff

    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    save_name = '2D1_' + var_amp + 'amp'
    for pp in range(len(path_list)):
        p_list = []
        leg_list = []
        c_list = cm.plasma(np.linspace(0, 1, len(time_list)))
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
            y_c_grid = dict_T['y_grid'] * xstep_factor
            min_data = np.min(T_data[cl_index:c_index, :])
            dict = cfoil.calc_norms_2(var, min_data, normal_class=cd5, forced_power=[0])
            norm_const = dict['norm_const']
            c_title = r'$\delta$ ' + dict['title']
            c_fmt = dict['c_fmt']
            var_name = dict['var_name']
            units = dict['units']

            final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
            final_lab = r'$' + final_lab + '$'

            dashes = (4, 2)

            lstyle = style_list[pp]
            mstyle = mstyle_list[pp]
            print 'lstyle = ', lstyle, mstyle
            yminax, ymaxax = ax.get_ylim()
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


#-----------------------------------------------------------------------
def plot_custom_xlineout(fig,
                         ax,
                         path_list,
                         var_amp='Te',
                         time='15',
                         vert_lines='False',
                         lstyle='-',
                         mstyle=None):
    '''
        p2 = plot_custom_xlineout(fig,ax,path_list,var_amp='Te')
    '''
    #SI_on = True
    #norm_name = '/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/chfoil_default5_norm.txt'
    #cd5 = cfoil.conv_factors_cd5(SI_on,norm_name)
    # --- init stuff
    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data
    #yi = 20
    #xi = 60 # arbitrary x index to plot at
    #c_index = 200

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    save_name = '2D1_' + var_amp + 'amp'
    style_list = ['-', '--', ':', '-', '--', ':']
    mstyle_list = [None, None, None, 'x', '^', 'o']
    for pp in range(len(path_list)):
        #print ' path = ', path_list[pp]
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
        #x_c_grid = x_c_grid[:len(T_data[cl_index:c_index,0])]
        y_c_grid = dict_T['y_grid'] * xstep_factor
        print ' x_c_grid shape = ', np.shape(x_c_grid), np.shape(T_data[cl_index:c_index])

        #ax.yaxis.set_ticks(np.array([0.6,0.8,1.0]))
        min_data = np.min(T_data[cl_index:c_index, :])
        dict = cfoil.calc_norms_2(var, min_data, normal_class=cd5, forced_power=[0])
        norm_const = dict['norm_const']
        c_title = dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
        final_lab = r'$' + final_lab + '$'
        #lstyle,marker = cfoil.get_path_style(final_lab,h_lab,B_lab)

        dashes = (4, 2)
        yminax, ymaxax = ax.get_ylim()

        p1, = ax.plot(x_c_grid[cl_index:c_index],
                      T_data[cl_index:c_index, yi] * norm_const,
                      c='k',
                      linestyle=lstyle,
                      marker=mstyle,
                      markersize=MSIZE,
                      markevery=3,
                      label=final_lab)

        if lstyle == '--':
            p1.set_dashes(dashes)

    ax.grid(color='0.5', linestyle='-')

    ax.set_xlabel(xlab)
    ax.set_ylabel(c_title)
    return p1


#-----------------------------------------------------------------------


def plot_Te_ylineout(fig, axy, path_list, var_amp='Te', time='15', dict_list=[]):
    '''
        Plots amplitude as a function of 
    '''

    if len(dict_list) != 0:
        color_lineout = dict_list['color_lineout']
        lineout_list = dict_list['lineout_list']
    #SI_on = True
    #norm_name = '/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/chfoil_default5_norm.txt'
    #cd5 = cfoil.conv_factors_cd5(SI_on,norm_name)
    # --- init stuff
    leg_list = []
    lab_list = []
    lim_data = 0.01
    min, max = 0.0, lim_data
    take_amp_on = True
    lab_type = 'B'

    n = 4    # number of ks

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    style_list = ['-', '--', ':', '-', '--', ':']
    mstyle_list = [None, None, None, 'x', '^', 'o']
    save_name = '2D1_' + var_amp + 'amp'
    for pp in range(len(path_list)):
        #print ' path = ', path_list[pp]
        var = var_amp

        fname = cfoil.construct_fname(path_list[pp], fpre(path_list[pp]), var, time)
        dict_T = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), var, time)
        T_data = dict_T['mat']
        time_col = float(dict_T['time']) * tstep_factor
        print '\n --> time = ', time_col, ' ps  = ', dict_T['time'], ' tcol <-----\n'
        x_c_grid = dict_T['x_grid'] * xstep_factor
        #x_c_grid = x_c_grid[:len(T_data[cl_index:c_index,0])]
        y_c_grid = dict_T['y_grid'] * xstep_factor

        min_data = np.min(T_data)
        dict = cfoil.calc_norms_2(var, min_data)
        norm_const = dict['norm_const']
        c_title = dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
        final_lab = r'$' + final_lab + '$'
        lstyle, marker = cfoil.get_path_style(final_lab, h_lab, B_lab)
        b00 = re.search('=0$', B_lab)
        #if B_lab == 'B' or B_lab== '0H':
        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]
        if lstyle == '--':
            dashes = (4, 2)
        else:
            dashes = (1, 0)
        print 'lstyle = ', lstyle, mstyle
        '''
        if B_lab == 'B' or b00:
           lstyle='-'
           mstyle = '^'
           
        else:
           lstyle='--'
           mstyle = 'x'
           dashes=(4,2)
        #if h_lab == 'static':
        #   lstyle=':'
        '''
        if not marker_on:
            mstyle = None
            print ' ---- no markesr on---- \n\n'

        for xl in range(len(lineout_list)):
            x_lineout = lineout_list[xl]
            ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            if take_amp_on:
                T_dev = T_data[lineout_list[xl], :] * norm_const - np.average(
                    T_data[lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data[lineout_list[xl], :]))
            else:
                T_dev = T_data[lineout_list[xl], :] * norm_const
            #p2, = axy.plot(y_c_grid,T_dev,c=color_lineout[xl],
            #                linestyle=lstyle,label=ylab_lineout)
            p2, = axy.plot(y_c_grid,
                           T_dev,
                           linestyle=lstyle,
                           c=color_lineout[xl],
                           marker=mstyle,
                           markersize=MSIZE,
                           markevery=3)
            if lstyle == '--':
                p2.set_dashes(dashes)
            p2.set_linewidth = 2
            leg2_list.append(p2)
            lab_lineout_list.append(ylab_lineout)
        leg_list.append(p2)
        if lab_type == 'B':
            if B_lab == 'B':
                lab_out = r'$B$'
            elif B_lab == r'no$ $B':
                lab_out = r'$no$ $B$'
                mstyle = 'x'
            else:
                lab_out = r'$' + B_lab + '$'
            lab_list.append(lab_out)
        else:
            lab_list.append(final_lab)

    if take_amp_on:
        ydifflab = var_name + r'$ - \langle $' + var_name + r'$ \rangle$ ' + units
    else:
        ydifflab = c_title
    axy.set_ylabel(ydifflab)
    axy.set_xlabel(ylab)
    axy.legend(leg_list, lab_list, loc='upper center', numpoints=1)
    axy.grid(color='0.5', linestyle='-')

    return p2


#-----------------------------------------------------------------------
def plot_ylineout_custom(fig,
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

    # ---- loop over paths
    leg2_list = []
    lab_lineout_list = []
    leg_list_x = []
    lab_list_x = []
    #style_list = ['-','--',':','-','--',':']
    #mstyle_list = [None,None,None,'x','^','o']
    save_name = '2D1_' + var_amp + 'amp'
    for pp in range(len(path_list)):
        var = var_amp

        fname = cfoil.construct_fname(path_list[pp], fpre(path_list[pp]), var, time)
        dict_T = cfoil.load_dict(path_list[pp], fpre(path_list[pp]), var, time)
        T_data = dict_T['mat']
        time_col = float(dict_T['time']) * tstep_factor
        print '\n --> time = ', time_col, ' ps  = ', dict_T['time'], ' tcol <-----\n'
        x_c_grid = dict_T['x_grid'] * xstep_factor
        y_c_grid = dict_T['y_grid'] * xstep_factor

        min_data = np.min(T_data)
        dict = cfoil.calc_norms_2(var, min_data)
        norm_const = dict['norm_const']
        c_title = dict['title']
        c_fmt = dict['c_fmt']
        var_name = dict['var_name']
        units = dict['units']

        final_lab, k_lab, h_lab, B_lab = cfoil.construct_label(fname)
        final_lab = r'$' + final_lab + '$'
        #b00 = re.search('=0$', B_lab)
        lstyle = style_list[pp]
        mstyle = mstyle_list[pp]
        for xl in range(len(lineout_list)):
            x_lineout = lineout_list[xl]
            ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            if take_amp_on:
                T_dev = T_data[lineout_list[xl], :] * norm_const - np.average(
                    T_data[lineout_list[xl], :] * norm_const) * np.ones(
                        len(T_data[lineout_list[xl], :]))
            else:
                T_dev = T_data[lineout_list[xl], :] * norm_const

            p2, = axy.plot(y_c_grid,
                           T_dev,
                           linestyle=lstyle,
                           c=color_lineout[xl],
                           marker=mstyle,
                           markersize=MSIZE,
                           markevery=3)
            #if lstyle=='--':
            #    p2.set_dashes(dashes)
            p2.set_linewidth = 2
            leg2_list.append(p2)
            lab_lineout_list.append(ylab_lineout)
        leg_list.append(p2)
        if lab_type == 'B':
            if B_lab == 'B':
                lab_out = r'$B$'
            elif B_lab == r'no$ $B':
                lab_out = r'$no$ $B$'
                mstyle = 'x'
            else:
                lab_out = r'$' + B_lab + '$'
            lab_list.append(lab_out)
        else:
            lab_list.append(final_lab)

    if take_amp_on:
        ydifflab = var_name + r'$ - \langle $' + var_name + r'$ \rangle$ ' + units
    else:
        ydifflab = c_title
    axy.set_ylabel(ydifflab)
    axy.set_xlabel(ylab)
    if leg_on:
        axy.legend(leg_list, lab_list, loc='upper center', numpoints=1)
    axy.grid(color='0.5', linestyle='-')

    return p2


#-----------------------------------------------------------------------
def fpre(path_in):
    return path_in.split('/')[-1]


#-----------------------------------------------------------------------
if __name__ == "__main__":

    fig, ax = fprl.newfig(1.0)
    plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})

    path1 = 'insert_test_path'
    path2 = 'insert_test_path'
    path_list = [path1, path2]
    #p1 = plot_Te_xlineout(fig,ax,path_list,var_amp='Te',time=loc_nspace.time_glb)

    if True:
        figy, axy = fprl.newfig(1.0)
        p2 = plot_Te_ylineout(figy, axy, path_list, var_amp='Te', time=loc_nspace.time_glb)

    plt.show()
