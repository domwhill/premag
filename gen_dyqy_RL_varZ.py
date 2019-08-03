'''
    plots the divergence of q_RL
'''
import numpy as np, re, os, sys, getpass,matplotlib.pyplot as plt
import getpass,site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')

import chfoil_module as cf
from chfoil_module import conv_factors_eos
from chfoil_module import cd5_switches
#import plot_comparison as pc
import figure_prl_twocol as fprl
import kinetic_ohmslaw_module_varZ as kohb
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
yi= cd5.yi

Bz_factor = cd5.Bz_ref
divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
lh_style = cd5.lh_style
dashes = cd5.dashes
MSIZE=4

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


def repack_2D(path,time):
    '''
        dict_c,dict_k = repack_2D(path,time)
    '''
    dict_qc,dict_qk = kohb.get_q_abs(path,time)#get_q_ratio(path,time)#get_q_individ(path,time)
    v_nx_k,v_ny_k,v_nx_c,v_ny_c = kohb.get_Nernst_abs(path,time)
    
    dict_out_c = {}
    dict_out_k = {}
    var_list = ['SH x','SH y','RL x','RL y','E x','E y','vN x','vN y']

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
    dict_out_k['Bz'] = dict_qk['Bz']*Bz_factor
    dict_out_k['x_grid'] = dict_qk['x_grid']*xstep_factor
    return dict_out_c, dict_out_k
#-----------------------------------------------------------------------

def load_qdata(path,time='10'):
    #-------- LOAD q data ----------------------------------------------
    kohnew = {}

    
    tt = int(time)
    string_tt_glb = '%2.2i' % int(time)
    #string_tt_glb = '%02i' %  int(time)
    time = string_tt_glb
    print ' ---------- tt ==== ', tt
    print '--------'
    print '\n\nkinetic model time = ', time
    dict = kohb.get_kinetic_heatflow_b(path,str(time))
    kohnew[path] = dict
        
    return kohnew


   

def extract_ohms(path,time,x_limit=73):
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
    dict_kinetic = kbier.get_kinetic_E(path,fprefix,time,xlim=73)
    

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

    func_l = lambda data: np.transpose(data)#[:,:]
    for key in dict.keys():
        dict[key]['data'] = func_l(dict[key]['data'])
    
    return dict
#-----------------------------------------------------------------------
def fpre(path):
    return path.split('/')[-1]



def plot_dyqy_RL(fig,ax,cax,path,time):

    kohnew = load_qdata(path,time)
    x_grid_SI = kohnew[path]['x_grid']*xstep_factor 
    y_grid_SI = kohnew[path]['y_grid']*xstep_factor       
    time_col = kohnew[path]['time']*tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8,8) # x,y inches
    lab_dict['lims'] = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'],lab_dict['ylab'] = xlab,ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'#'hot'
    var = 'q_RL'
    print ' shape = ', np.shape(kohnew[path][var + '_x']['data'])
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #data_x = np.transpose(kohnew[path][var +'_x']['data'])#q_SH_x[path][tt,:,:]
    data_y = np.transpose(kohnew[path][var + '_y']['data'])#q_SH_y[path][tt,:,:]
    dxqy,dyqy =  cf.get_grad(x_grid_SI,y_grid_SI,data_y)
    dyqy = np.transpose(dyqy)#q_SH_y[path][tt,:,:]
   

    dyqy *= divq_factor*0.1
    
    
    MULT = -1.0
    #vmin,vmax = np.min(MULT*dyqy),np.max(MULT*dyqy)
    #vmin,vmax = -1,2
    multtemp=0.9
    vmin,vmax = np.min(MULT*dyqy[:,cl_index:c_index])*multtemp,np.max(MULT*dyqy[:,cl_index:c_index])*multtemp
    print ' vmin, vmax = ', vmin, vmax
    print ' np.shape(dyqy) = ', np.shape(dyqy)
    #sys.exit()
    #vmin, vmax =  -0.00189426831011, 0.000207009514978
    lims_im = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[cl_index],x_grid_SI[c_index]]#
    lims_rev = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[cl_index]]
    lims_im = [x_grid_SI[cl_index],x_grid_SI[c_index],y_grid_SI[0],y_grid_SI[-1]]#

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(MULT*dyqy[:,cl_index:c_index],
                    cmap=cmap,vmin=vmin,vmax=vmax,
                    norm=norm,
                    aspect=asp,extent=lims_im)
    levels = np.arange(vmin,vmax,(vmax-vmin)/4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        #print 'levels ll = ', ll, levels[ll]
        levels[ll] *= divq_factor

    #
    levels = [0.0] 
    claby = r'-$\partial_y q_{y,RL}$ ' #+ divq_unit
    CSy = ax.contour(MULT*dyqy[::-1,cl_index:c_index]*divq_factor,levels=levels,colors=('k'),extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')
    
    #ax.set_title('dyqy at ' + str(time_col))
    
    #--- plot hlines
    lim = ax.get_ylim()
    if not lh_style=='None':
        cf.plot_xvlines(ax,x_grid_SI,lim,lineout_list,linestyle=lh_style,dashes=dashes)
    
    ax.set_xlabel(xlab_rel)
    ax.set_ylabel(ylab)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #---
    c2 = fig.colorbar(imy,cax=cax,ax=ax,
                        format= '%1i',
                        label=claby)
    c2.add_lines(CSy)
    
    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()
    #c2 = fig.colorbar(CSy,cax=cax,ax=ax,
    #                    label=claby)
    
    #c2.ax.plot([0,1],[0,0],'k')
    #print ' dir = ', dir(c2.ax)
    return imy


def plot_dyqy_RL_c(fig,ax,cax,path,time):

    kohnew = load_qdata(path,time)
    dict_c, dict_k = repack_2D(path,time)

    x_grid_SI = kohnew[path]['x_grid']*xstep_factor 
    y_grid_SI = kohnew[path]['y_grid']*xstep_factor       
    time_col = kohnew[path]['time']*tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8,8) # x,y inches
    lab_dict['lims'] = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'],lab_dict['ylab'] = xlab,ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'#'hot'
    var = 'q_RL'
    print ' shape = ', np.shape(kohnew[path][var + '_x']['data'])
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #data_x = np.transpose(kohnew[path][var +'_x']['data'])#q_SH_x[path][tt,:,:]
    data_y = np.transpose(kohnew[path][var + '_y']['data'])#q_SH_y[path][tt,:,:]
    
    data_yc = dict_c['RL y']
    print('shape classical = ', np.shape(data_yc),'shape kinetic qrl = ', np.shape(data_y))
    data_y = data_yc*1.0
    dxqy,dyqy =  cf.get_grad(x_grid_SI,y_grid_SI,data_y)
    dyqy = np.transpose(dyqy)#q_SH_y[path][tt,:,:]
   

    dyqy *= divq_factor*0.1
    
    
    MULT = -1.0
    #vmin,vmax = np.min(MULT*dyqy),np.max(MULT*dyqy)
    #vmin,vmax = -1,2
    multtemp=0.9
    vmin,vmax = np.min(MULT*dyqy[:,cl_index:c_index])*multtemp,np.max(MULT*dyqy[:,cl_index:c_index])*multtemp
    print ' vmin, vmax = ', vmin, vmax
    print ' np.shape(dyqy) = ', np.shape(dyqy)
    #sys.exit()
    #vmin, vmax =  -0.00189426831011, 0.000207009514978
    lims_im = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[cl_index],x_grid_SI[c_index]]#
    lims_rev = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[cl_index]]
    lims_im = [x_grid_SI[cl_index],x_grid_SI[c_index],y_grid_SI[0],y_grid_SI[-1]]#

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(MULT*dyqy[:,cl_index:c_index],
                    cmap=cmap,vmin=vmin,vmax=vmax,
                    norm=norm,
                    aspect=asp,extent=lims_im)
    levels = np.arange(vmin,vmax,(vmax-vmin)/4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        #print 'levels ll = ', ll, levels[ll]
        levels[ll] *= divq_factor

    #
    levels = [0.0] 
    claby = r'-$\partial_y q_{y,RL,c}$ ' #+ divq_unit
    CSy = ax.contour(MULT*dyqy[::-1,cl_index:c_index]*divq_factor,levels=levels,colors=('k'),extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')
    
    #ax.set_title('dyqy at ' + str(time_col))
    
    #--- plot hlines
    lim = ax.get_ylim()
    if not lh_style=='None':
        cf.plot_xvlines(ax,x_grid_SI,lim,lineout_list,linestyle=lh_style,dashes=dashes)
    
    ax.set_xlabel(xlab_rel)
    ax.set_ylabel(ylab)
    #---
    c2 = fig.colorbar(imy,cax=cax,ax=ax,
                        format= '%1i',
                        label=claby)
    c2.add_lines(CSy)
    
    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()

    return imy



def plot_dBdt_bier(fig,ax,cax,path,time):

    #kohnew = load_qdata(path,time)
    kohnew = kbier.get_kinetic_E(path,fpre(path),time,xlim=73)
    x_grid_SI = kohnew[path]['x_grid']*xstep_factor 
    y_grid_SI = kohnew[path]['y_grid']*xstep_factor       
    time_col = kohnew[path]['time']*tstep_factor
    asp = 'auto'
    lab_dict = {}
    lab_dict['figsize'] = (8,8) # x,y inches
    lab_dict['lims'] = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'],lab_dict['ylab'] = xlab,ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    cmap = 'RdBu_r'#'hot'
    var = 'q_RL'
    print ' shape = ', np.shape(kohnew[path][var + '_x']['data'])
    tt = int(time)

    #--- do everything in units of 10 eV/ps
    #data_x = np.transpose(kohnew[path][var +'_x']['data'])#q_SH_x[path][tt,:,:]
    data_y = np.transpose(kohnew[path][var + '_y']['data'])#q_SH_y[path][tt,:,:]
    dxqy,dyqy =  cf.get_grad(x_grid_SI,y_grid_SI,data_y)
    dyqy = np.transpose(dyqy)#q_SH_y[path][tt,:,:]
   

    dyqy *= divq_factor*0.1
    
    
    MULT = -1.0
    #vmin,vmax = np.min(MULT*dyqy),np.max(MULT*dyqy)
    #vmin,vmax = -1,2
    multtemp=0.9
    vmin,vmax = np.min(MULT*dyqy[:,cl_index:c_index])*multtemp,np.max(MULT*dyqy[:,cl_index:c_index])*multtemp
    print ' vmin, vmax = ', vmin, vmax
    print ' np.shape(dyqy) = ', np.shape(dyqy)
    #sys.exit()
    #vmin, vmax =  -0.00189426831011, 0.000207009514978
    lims_im = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[cl_index],x_grid_SI[c_index]]#
    lims_rev = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[cl_index]]
    lims_im = [x_grid_SI[cl_index],x_grid_SI[c_index],y_grid_SI[0],y_grid_SI[-1]]#

    norm = cf.MidPointNorm(0.0)
    imy = ax.imshow(MULT*dyqy[:,cl_index:c_index],
                    cmap=cmap,vmin=vmin,vmax=vmax,
                    norm=norm,
                    aspect=asp,extent=lims_im)
    levels = np.arange(vmin,vmax,(vmax-vmin)/4.0)
    levels = levels.tolist()

    #levels= [-0.001]
    for ll in range(len(levels)):
        #print 'levels ll = ', ll, levels[ll]
        levels[ll] *= divq_factor

    #
    levels = [0.0] 
    claby = r'-$\partial_y q_{y,RL}$ ' #+ divq_unit
    CSy = ax.contour(MULT*dyqy[::-1,cl_index:c_index]*divq_factor,levels=levels,colors=('k'),extent=lims_im)
    #ax.clabel(CSy,[0.0],inline=True,fmt='%1.1f')
    #ax.clabel(CSy,levels,inline=True,fmt='%3.4f')
    
    #ax.set_title('dyqy at ' + str(time_col))
    
    #--- plot hlines
    lim = ax.get_ylim()
    if not lh_style=='None':
        cf.plot_xvlines(ax,x_grid_SI,lim,lineout_list,linestyle=lh_style,dashes=dashes)
    
    ax.set_xlabel(xlab_rel)
    ax.set_ylabel(ylab)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #---
    c2 = fig.colorbar(imy,cax=cax,ax=ax,
                        format= '%1i',
                        label=claby)
    c2.add_lines(CSy)
    
    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()
    #c2 = fig.colorbar(CSy,cax=cax,ax=ax,
    #                    label=claby)
    
    #c2.ax.plot([0,1],[0,0],'k')
    #print ' dir = ', dir(c2.ax)
    return imy



def plot_ylineout_qRLy_custom(fig,axy,path_list,var_amp='Te',time='15',mstyle_list = [None,None,None,'x','^','o'],style_list = ['-','--',':','-','--',':'],leg_on=True, dict_list = []):
    # --- init stuff
    if len(dict_list)!=0:
        color_lineout = dict_list['color_lineout']
        lineout_list = dict_list['lineout_list']
    
    leg_list = []
    lab_list = []
    lim_data = 0.01
    min,max = 0.0, lim_data 
    take_amp_on = True
    lab_type = 'B'
    n =4 # number of ks    
    
    
    
    #===
    
    cmap = 'RdBu_r'#'hot'
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
        kohnew = load_qdata(path,time)
        x_grid_SI = kohnew[path]['x_grid']*xstep_factor 
        y_grid_SI = kohnew[path]['y_grid']*xstep_factor       
        time_col = kohnew[path]['time']*tstep_factor
        asp = 'auto'

        dict_c, dict_k = repack_2D(path,time)
        T_data_c = dict_c['RL y']
        
        cmap = 'RdBu_r'#'hot'
        var = 'q_RL'
        tt = int(time)

        #--- do everything in units of 10 eV/ps
        T_data = np.transpose(kohnew[path][var + '_y']['data'])#q_SH_y[path][tt,:,:]
        #<<<<---------------------

        MULT = -1.0
        multtemp=0.9
        
 
        
        dict_T = cf.load_dict(path,fprefix,'Te',time)
        time_col = float(dict_T['time'])*tstep_factor
        print '\n --> time = ',time_col, ' ps  = ', dict_T['time'], ' tcol <-----\n'
        x_c_grid =dict_T['x_grid']*xstep_factor
        y_c_grid = dict_T['y_grid'][1:-1]*xstep_factor
        
       
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
        print( ' shape qrL = ',np.shape(y_c_grid))
        print( ' shape qrL = ',np.shape(T_data))
        for xl in range(len(lineout_list)):
            x_lineout = lineout_list[xl]
            ylab_lineout = r'$%3.1f$' % x_c_grid[x_lineout]
            if take_amp_on:
                T_dev = T_data[lineout_list[xl],:]*norm_const - np.average(T_data[lineout_list[xl],:]*norm_const)*np.ones(len(T_data[lineout_list[xl],:]))

                T_dev_c = T_data_c[lineout_list[xl],:]*norm_const - np.average(T_data_c[lineout_list[xl],:]*norm_const)*np.ones(len(T_data_c[lineout_list[xl],:]))



            else:
                T_dev = T_data[lineout_list[xl],:]*norm_const
                T_dev_c = T_data_c[lineout_list[xl],:]*norm_const

            p2, = axy.plot(y_c_grid,T_dev,
                           linestyle=lstyle,
                           c=color_lineout[xl],
                           marker=mstyle,markersize=MSIZE, markevery=3)
            p2c, = axy.plot(y_c_grid,T_dev_c,
                           linestyle='--',
                           c=color_lineout[xl],
                           marker=mstyle,markersize=MSIZE, markevery=3)

            p2.set_linewidth=2
            leg2_list.append(p2)
            lab_lineout_list.append(ylab_lineout)
        leg_list.append(p2)
        lab_list.append(final_lab)
    
    if take_amp_on:
        ydifflab = var_name  + r'$ - \langle $' + var_name + r'$ \rangle$ ' + units
        ydifflab = r'$\delta q_{y,RL}$'
    else:
        ydifflab = final_lab
    axy.set_ylabel(ydifflab)
    axy.set_xlabel(ylab)
    if leg_on:
        leg_list = [p2c,p2]
        lab_list = [r'$\delta q_{y,RL,c}$', r'$\delta q_{y,RL,k}$']
        #axy.legend(leg_list,lab_list,numpoints=1)
    axy.grid(color='0.5',linestyle='-')

    return p2
    

#-----------------------------------------------------------------------
if __name__=="__main__":
    #-- generate new axis for color bar with gridspec here....
    fig, ax  = fprl.newfig(1.0)
    path, time = 'chfoil_default5_heat11_2D1--cx1','15'
    #matplotlib.rc('font', **{'family':"sans-serif"})

    params = {'text.latex.preamble': [r'\usepackage{siunitx}']}
    plt.rcParams.update(params)
    
    cax = []
    imy = plot_dyqy_RL(fig,ax,cax,path,time)
    plt.show()