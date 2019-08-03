'''

    21/3/2019 - DWH
    This script plots the ratio between the heat flow angles and theta_k/theta_c
'''
import numpy as np, sys, os, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import matplotlib.pyplot as plt
#---> kinetic/classical transport post processing module
import kinetic_ohmslaw_module_varZ as kohb

import figure as fprl
import matplotlib.gridspec as GS
import matplotlib.ticker as ticker
from pylab import *
import chfoil_module as cf


path_pre = '/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/PREMAG/2D_RUNS/'
path_0T = path_pre + 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_0T_E'
path_50T = path_pre + 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E'
path_400T = path_pre + 'r5_v40_vmax20_vnonuni10_400T_Zb'


path = path_50T
save_path = '/Users/' + userid + '/Dropbox/York/UPDATE/PICS_TEMP/'
#---------------------------------->
SI_on = cf.cd5_switches.SI_on
save_on = cf.cd5_switches.save_on
hlines_on = cf.cd5_switches.hlines_on
grid_on = cf.cd5_switches.grid_on
horizontal_on = cf.cd5_switches.horizontal_on
separate_plots_on = cf.cd5_switches.separate_plots_on

norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' +norm_name + 'norm.log' 
[T0,n0,Z0,Bz0] = np.loadtxt(log_file)

cd5 = cf.conv_factors_custom(norm_path,Z0,Ar=6.51)
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

divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
lh_style = cd5.lh_style
dashes = cd5.dashes
#-----------------------------?


fig = fprl.newfig_generic_3yscale(width=1.0,scale_width=1.0,scale_ratio=0.5,clearon=True)
[ax1,ax2,ax3] = [fig.add_subplot(2,3,1),fig.add_subplot(2,3,6),fig.add_subplot(2,3,2)]
[ax4,ax5,ax6]= [fig.add_subplot(2,3,4),fig.add_subplot(2,3,3),fig.add_subplot(2,3,5)]

#[ax1,ax2,ax3] = [fig.add_subplot(1,3,1),fig.add_subplot(1,3,2),fig.add_subplot(1,3,3)]


ax_list = [ax1,ax2,ax3]
map(fprl.fix_ax,[ax1,ax2,ax3])
#plt.subplots_adjust(hspace=0.1,left=0.2,right=0.95,wspace=0.68)
fig.subplots_adjust(hspace=0.5,left=0.15,right=0.96,wspace=0.5,bottom=0.2,top=0.9)

t_list = [3,6]
fcmap = cm.plasma
c_list = iter(fcmap(np.linspace(0,1,len(t_list))))
log_on = True
lab_list = []
plot_list = []


dict_fo = cf.load_dict(path,path.split('/')[-1],'fo','00')
x_grid_SI = dict_fo['x_grid']*xstep_factor
y_grid_SI = dict_fo['y_grid']*xstep_factor


for tt in t_list:
    c= next(c_list)
    print 'tt = ', tt
    time = '%02i' % tt

    dict_ratio = kohb.get_q_ratio(path,time)
    data_x2d_SH,data_y2d_SH,prelab_SH = dict_ratio['SH x'],dict_ratio['SH y'],'\perp'
    data_x2d_RL,data_y2d_RL,prelab_RL = dict_ratio['RL x'],dict_ratio['RL y'],'RL'
    data_x2d_vN,data_y2d_vN = kohb.get_Nernst_ratio(path,time)
    
    
    #rat_vN = 
    #rat_RL = 
    
    
    prelab_vN = 'v_N'
    #prelab = 'v_{N'
    data_xSH = cf.get_avgx(data_x2d_SH)
    data_ySH = cf.get_avgx(data_y2d_SH)
    data_xRL = cf.get_avgx(data_x2d_RL)
    data_yRL = cf.get_avgx(data_y2d_RL)
    data_xvN = cf.get_avgx(data_x2d_vN)
    data_yvN = cf.get_avgx(data_y2d_vN)
    
    
    theta_SH = cf.get_avgx(data_y2d_SH/data_x2d_SH)
    
    #theta_SH = cf.get_avgx(data_y2d_SH/data_x2d_SH)
    #p1, =ax1.plot(x_grid_SI[:x_lim],data_x,c=c)

    #- This one would plot the angles of stuff...
    #p1, =ax1.semilogy(x_grid_SI,data_ySH/data_xSH,c=c)
    #p2, =ax2.semilogy(x_grid_SI,data_yvN/data_xvN,c=c)  
    #p3, =ax3.semilogy(x_grid_SI,data_yRL/data_xRL,c=c)
    
    p1, =ax1.plot(x_grid_SI,data_xSH,c=c)
    p2, =ax2.plot(x_grid_SI,data_xvN,c=c)
    p3, =ax3.plot(x_grid_SI,data_yRL,c=c)
    
    p4, =ax4.plot(x_grid_SI,data_ySH,c=c)    
    p5, =ax5.plot(x_grid_SI,data_yvN,c=c)
    p6, =ax6.plot(x_grid_SI,data_xRL,c=c)    

    #p3, =ax3.semilogy(x_grid_SI,data_yRL,c=c)

    
    #ax3.set_ylim(0.1,5.0)
    #--- legend stuff
    dict_wt = cf.load_dict(path,path.split('/')[-1],'wt',time)
    t_l = dict_wt['time']*tstep_factor
    t_lab = '%3i \si{ps}' % t_l
    t_lab = '%3i' % t_l
    
    plot_list.append(p1)
    lab_list.append(t_lab)

#xtit = lambda ax, name: ax.set_ylabel('$ %s,x,K}/ %s,x,B}$' % (name,name))
ytit = lambda ax, name: ax.set_ylabel('$\Theta_{%s,K}/ \Theta_{%s,B}$' % (name,name))
ytit = lambda ax, name: ax.set_ylabel('$\frac{q_{y,RL,K}}{q_{y,\perp,K}}\frac{q_{y,\perp,c}}{q_{y,RL,c}}$')
ytit = lambda ax, name: ax.set_ylabel(r'$(q_{y,RL,K}/q_{y,\perp,K})(q_{y,\perp,c}/q_{y,RL,c})$')
ytit = lambda ax, name: ax.set_ylabel('$%s$' % (name))

def clearxaxis(ax):
    ax.set_xlabel('')
    ax.set_xticklabels([])        
    
map(ytit,[ax1,ax2,ax3],[prelab_SH,prelab_vN,prelab_RL])
#map(clearxaxis,[ax1,ax2,ax3])
map(lambda axy: axy.set_xlabel(xlab_rel), [ax1,ax2,ax3])
map(lambda axy: axy.xaxis.set_major_locator(MaxNLocator(4,prune='both')), [ax1,ax2,ax3])
#map(lambda axy: axy.set_xlabel(xlab), [ax1,ax2,ax3,ax4,ax5,ax6])

xlim = ax1.get_xlim()
map(lambda axy: axy.set_xlim(-5.0,20.0),ax_list)
map(lambda axy: axy.set_xlim(-5.0,20.0),[ax4,ax5,ax6])

# - qSHx, vNx
def sort_y(axy):
    #for axis in [ax.xaxis, ax.yaxis]:
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    axy.yaxis.set_major_formatter(formatter)

    #axy.set_ylim(0.1,2.0)
    #axy.set_ylabel('')
    #axy.set_ylim(0.5e-2,1.5e2)
    #axy.set_yticks([0.5,1.0])


def sort_y_ylim(axy):
    #axy.set_ylim(0.5e-2,1.5e2)
    axy.set_yticks([0.1,1.0,10.0])
    
    #axy.set_ylabel('')
    #axy.yaxis.labelpad= 4
    #axy.yaxis.set_ticklabels([0.2,1.0,2.0])

if False:
    ax5.legend(plot_list,lab_list,
            bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
            title=r't $\si{ps}$')
savename =  save_path + 'theta_50T_both.pdf'
print('--- saving as : ',savename)
fig.savefig(savename)
plt.show()
    
