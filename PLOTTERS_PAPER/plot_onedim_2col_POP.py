'''#
    This file plots the inital profiles then the evolved profiles for a 50T
    sim. in format for 2column paper
'''


import os,sys,site,getpass
userid = getpass.getuser()
site.addsitedir('/Users/'+ userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/MODULES')
import sys, os, re
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import house_keeping as hk
import chfoil_module as cf
import impact_norms as inorm
import figure_prl_twocol as fprl
import matplotlib as mpl

paths = hk.directory_paths()
norm_path = paths.norm_dir
cfg = cf.conv_factors_custom(norm_path)

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
#---------------

mult = 6
cfg.yi = 1


#--->
#--- generating  path_list
dim = '1D'
bz_in = 50
save_tag = '%s_%iT_' % (sys.argv[0].split('.')[0],bz_in)
s_len = 1
path1 = paths.get_path(s_len, 50, lambda_p=5, dim=dim)
path2 = paths.get_path(s_len, 400, lambda_p=5, dim=dim)
path_list = [path1]

save_path = paths.save_dir

xmin,xmax=-10.0,25.0
#-----------------------------------------------------------------------
    
def plot_ax(ax1,x_grid,data,norm_const,c='b',tlab='00',cfg=cfg):
    if len(np.shape(data))>1:
        if len(data[0,:])<3:
            x_gr= x_grid[:len(data[:,0])]
            ax1.plot(x_gr,data[:,0]*norm_const,c=c,label=r'\textit{$' + str(tlab) + '$}') 
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        else:
            x_gr= x_grid[:len(data[:,0])]
            ax1.plot(x_gr,data[:,cfg.yi]*norm_const,c=c,label=r'\textit{$' + str(tlab) + '$}') 
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
    else:
        x_gr= x_grid[:len(data)]
        ax1.plot(x_gr,data[:,cfg.yi]*norm_const,c=c,label=r'\textit{$' + str(tlab) + '$}') 

    return    
#-----------------------------------------------------------------------
if __name__=="__main__":
    var_list = ['Cx','n','Te','Bz','wt']
    t_list = []
    tc_list = []
    
    #-------------------------------------------------------------------
    fig = fprl.newfig_generic(1.4,0.9)
    plt.subplots_adjust(left=0.1,right=0.92,bottom=0.18,top=0.9,wspace=0.1)
    
    axn0 = fig.add_subplot(121)
    ax20 = axn0.twinx()
    axn1 = fig.add_subplot(122)
    ax21 = axn1.twinx()
    ax_list = [(axn0,ax20),(axn1,ax21)]
    if len(t_list)>5:
        loop = range(0,len(t_list),len(t_list)/5)
    else:
        loop = range(0,len(t_list))
    loop = [0,7,9,10,11]
    n = len(loop)# number of curves to plot
    ##color=iter(cm.plasma(np.linspace(0,1,n)))
        
    if not os.path.isdir(save_path):
        os.system('mkdir ' + save_path)
        print 'making directory ', save_path
    lstyle_list=['-','--',':']
    if False:
        var_extra = 'wt'
        var_extra_ylab = r'\omega \tau_{ei}'
        var_extra_lab = r'$\omega \tau_{ei}$'
    elif False:
        var_extra = 'Z'
        var_extra_lab = r'Z'
        var_extra_ylab = r'Z'
    else:
        var_extra = ''
        var_extra_lab = r''
        var_extra_ylab = r''
    
    tt_list= [0,15]
    for path in path_list:
        for itt, tt in enumerate(tt_list):
            axn = ax_list[itt][0]
            ax2 = ax_list[itt][1]
            #lstyle=lstyle_list[path_list.index(path)]
            lstyle=lstyle_list[path_list.index(path)]

            #--- load data
            time = '%02i' % tt
            fprefix = path.split('/')[-1]
            Te_dict = cf.load_dict(path,fprefix,'Te',time)
            n_dict = cf.load_dict(path,fprefix,'n',time)
            Cx_dict = cf.load_dict(path,fprefix,'Cx',time)
            Bz_dict = cf.load_dict(path,fprefix,'Bz',time)
        
            T_data = Te_dict['mat']
            n_data = n_dict['mat']
            Cx_data = Cx_dict['mat']
            Bz_data = Bz_dict['mat']
            time_col = Te_dict['time']
        
            x_grid_ccg = n_dict['x_grid']*cfg.xstep_factor
            power = 0#cf.extract_power(ne_dict['norm_const'])
            mod = ''#r'$ 10^{' +str(power)  + '} $'
        
            n_mult,T_mult,C_mult,B_mult = 1.0/40.0,2.5,20.0,1.2
        
            #n_str = get_norm_const(cfg.n0*(n_mult**-1))
            #T_str = get_norm_const((T_mult**-1)*2.0*Te_ref*1e-3)
            #C_str = get_norm_const((v_te*1e-3)*(C_mult**-1))
            #I_str = get_norm_const(B_const*(B_mult**-1))
        
        
            ne_lab = r'$n_e$[\SI{%1.1e}{\centi\meter^{-3}}]' % (cfg.n0*(n_mult**-1))
            Te_lab = r'$T_e$[\SI{%1.1f}{\kilo\electronvolt}]' % ((T_mult**-1)*2.0*cfg.T0*1e-3)
            Cx_lab = r'$C_x$[\SI{%1.1e}{\kilo\meter\per\second}]' % ((cfg.v_te*1e-3)*(C_mult**-1))
            Bz_lab = r'$B_z$ [\SI{%1.1e}{T}]' % ((cfg.Bz0**-1)*(B_mult**-1))

        
            ###n_crit_norm_oneplot = n_crit_norm*n_mult
            ne_norm = n_data[cfg.cl_index:cfg.c_index,cfg.yi]*n_mult
            Te_norm = T_data[cfg.cl_index:cfg.c_index,cfg.yi]*T_mult
            C_norm = Cx_data[cfg.cl_index:cfg.c_index,cfg.yi]*C_mult
            B_norm = Bz_data[cfg.cl_index:cfg.c_index,cfg.yi]*B_mult
            #axI = axn.twinx()
            b1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],B_norm,'k',linestyle=lstyle,label='$B_z$')
            n1, = ax2.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],ne_norm,c='g',linestyle=lstyle,label=r'$n_e$')
            t1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],Te_norm,c='r',linestyle=lstyle,label='$T_e$')
            Cx1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],C_norm,c='b',linestyle=lstyle,label='$C_x$')
            
            if var_extra != '':
                wt_dict = cf.load_dict(path,fprefix,var_extra,time)
                wt_data = wt_dict['mat']
                wt_lab = var_extra_lab
                wt_norm = wt_data[cfg.cl_index:cfg.c_index,cfg.yi]
                wt1, = axn.plot(x_grid_ccg[cfg.cl_index:cfg.c_index],wt_norm,c='k',linestyle=lstyle,label='$\omega \tau$')
            ylim = axn.get_ylim()
            axn.set_ylim(ylim[0],ylim[1]*1.2)
            axn.set_xlim(xmin,xmax)
            time_si = r'%i $\si{ps}$' % (time_col*cfg.tstep_factor)
            if itt==1:
                dx_mult = 0.38
                dy_mult = 0.28
            else:
                dx_mult = 0.3
                dy_mult = 0.24
            cf.annotate_time(
                    axn, lett=time_si,
                    dx_mult=dx_mult,dy_mult=dy_mult,
                    loc='top', fontsize=8
                    )

        # axn.axhline(n_crit_norm_oneplot,c='k',linestyle='--')
    y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)

    yn1,yn2 = axn1.get_ylim()
    ya1,ya2 = ax21.get_ylim()

    axn0.set_ylim(yn1,yn2)
    ax20.set_ylim(ya1,ya2)
        
    ax_all = [axn0, ax20, axn1, ax21]
    map(lambda ax: ax.tick_params(direction='in'),ax_all)
    if var_extra != '':
        p_list = [n1,t1,Cx1,b1,wt1]
        leg_lab_list = [ne_lab,Te_lab,Cx_lab,Bz_lab,wt_lab]
    else:
        p_list = [n1,t1,Cx1,b1]
        leg_lab_list = [ne_lab,Te_lab,Cx_lab,Bz_lab]

    leg =axn0.legend(p_list,leg_lab_list,ncol=1,fontsize=7,loc='upper left',
            framealpha=1.0,frameon=False)

    axn0.set_xlabel(r'$x-x_{abl}$ [$\si{\micro\meter}$ ]')
    axn1.set_xlabel(r'$x-x_{abl}$ [$\si{\micro\meter}$ ]')
    var_extra_ylab = ',\,' + var_extra_ylab
    axn0.set_ylabel(r'$T_e,\, C_x,\, B_z' + var_extra_ylab + '$')
    axn1.set_yticks([])
    ax20.set_yticks([])

    ax21.set_ylabel(r'$n_e$')
    axn0.set_xlim(xmin,xmax)
    axn1.set_xlim(xmin,xmax)
    
    
    #plt.show()
    save_suffix = 'one_dim_2col' + time
    savename = save_path + save_suffix + '.png'
    print 'saving fig as ', savename
    print 'copy to terminal: \n open -a preview ' + savename
    
    plt.savefig(savename ,dpi=400)
    
    axn.cla()
    ax2.cla()