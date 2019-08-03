'''
test plot
'''

import numpy as np
import matplotlib.pyplot as plt
import sys, re, os, getpass,site
userid = getpass.getuser()
site.addsitedir('/Users/'+ userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import chfoil_module as cf
#import q_SH_Te_tsteps as q_mod
import matplotlib as mpl
from pylab import *
import figure_prl_twocol as fprl
import kinetic_ohmslaw_module_1D_varZ as q_mod
from matplotlib.legend_handler import HandlerBase

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12

norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' +norm_name + 'norm.log' 
[T0,n0,Z0,Bz0] = np.loadtxt(log_file)

cd5 = cf.conv_factors_custom(norm_path,Z0,Ar=6.51)
print 'Z0 = ', Z0
print 'T0 = ', T0
print 'n0 = ', n0
print ' Bz0 = ', Bz0


cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index,c_index = 0,-1
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab_rel
ylab = cd5.ylab
Bz0  = (m_e/(q_e*tau_ei))
Bz_lab = r'$B_z$[\si{T}]'
Bz_factor = Bz0
#<--------------
# inputs

                                                                                                                    #r5_v40_Z_FEOS_trackncrit_D - 400T
#r5_v40_Z_FEOS_trackncrit_50T_FF - 50T
#r5_v40_Z_FEOS_trackncrit_1D_0T - 0T
#-------------------------------------------->
path_pre = '../PREMAG/'
path = path_pre + 'r5_v40_Z_FEOS_trackncrit_50T_FF'
path_0T = path_pre + 'r5_v40_Z_FEOS_trackncrit_1D_0T'
path_400T = path_pre + 'r5_v40_Z_FEOS_trackncrit_D'

path_0T = path_pre + 'r5_v40_Z_FEOS_trackncrit_1D_0T_2'
path_50T = path_pre + 'r5_v40_Z_FEOS_trackncrit_50T_FF'
path_400T = path_pre + 'r5_v40_Z_FEOS_MODNE_1D_400T_hires'#'r5_v40_Z_FEOS_trackncrit_D'
Nernst_on = False
plot_zero = True

time = '14'
iy = 1
xmin,xmax = -5.0,20.0
thresh_N_ratio = 5e-2
save_path = '/Users/' + userid + '/Dropbox/York/UPDATE/PICS_POP/'
save_name = save_path + 'qperp_ratio_1d.png'
#--------------------------------------------->
#--- functions
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,
                       x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0,y0+width], [0.7*height,0.7*height],
                           linestyle=orig_handle[1], color='k')
        l2 = plt.Line2D([x0,y0+width], [0.3*height,0.3*height], 
                           color=orig_handle[0])
        return [l1, l2]

#def get_ratio():
#--------------- 0T loading--->
dict_qc,dict_qk = q_mod.get_q_abs(path_0T,time)#get_q_ratio(path,time)#get_q_individ(path,time)
v_nx_k,v_ny_k,v_nx_c,v_ny_c = q_mod.get_Nernst_abs(path_0T,time)
x_grid = dict_qk['x_grid']*xstep_factor
#<------------------ SH + RL--
Te_0 = dict_qk['Te'][:,iy]
Bz_0 = dict_qk['Bz'][:,iy]*Bz_factor
#rat_0T = q_mod.get_ratio_lim(dict_qc['SH x'][:,0],dict_qk['SH x'][:,iy],vmax=100.0)
if Nernst_on:
    
    rat_0T = q_mod.get_ratio_lim(v_nx_c[:,iy],v_nx_k[:,iy],vmax=100.0,thresh=thresh_N_ratio)
    print(' shapes = ',np.shape(rat_0T), np.shape(v_nx_k),np.shape(v_nx_c))
    #sys.exit()
else:
    rat_0T = q_mod.get_ratio_lim(dict_qc['SH x'][:,0],dict_qk['SH x'][:,iy],vmax=100.0)

#<--------------- 0T loading---



#--------------- 50T loading--->
dict_qc,dict_qk = q_mod.get_q_abs(path,time)#get_q_ratio(path,time)#get_q_individ(path,time)
v_nx_k,v_ny_k,v_nx_c,v_ny_c = q_mod.get_Nernst_abs(path,time)
x_grid = dict_qk['x_grid']*xstep_factor

#--------> debug

print(np.shape(dict_qc['SH x']), np.shape(dict_qk['SH x']), np.shape(dict_qk['tot x']), np.shape(dict_qk['tot x']), np.shape(dict_qc['tot x']))
print(np.shape(dict_qc['SH x']), np.shape(dict_qk['SH x']), np.shape(dict_qk['tot x']), np.shape(dict_qk['tot x']), np.shape(dict_qc['tot x']))
print(np.shape(x_grid))
print('v_nx = ', np.shape(v_nx_k),np.shape(v_nx_k),np.shape(v_ny_k),np.shape(v_nx_c),np.shape(v_ny_c))

print(' Bz = ', np.shape(dict_qk['Bz']), 'np.shape(Te = ', np.shape(dict_qk['Te']))

#<------------------ SH + RL--
Te_50 = dict_qk['Te'][:,iy]
Bz_50 = dict_qk['Bz'][:,iy]*Bz_factor
if Nernst_on:
    
    v_nx_c50 = v_nx_c[:,iy]
    v_nx_k50 = v_nx_k[:,iy]
    
    rat_50T = q_mod.get_ratio_lim(v_nx_c[:,iy],v_nx_k[:,iy],vmax=100.0,thresh=thresh_N_ratio)
else:
    rat_50T = q_mod.get_ratio_lim(dict_qc['SH x'][:,0],dict_qk['SH x'][:,iy],vmax=100.0)
#--------------- 50T loading--->
dict_qc,dict_qk = q_mod.get_q_abs(path_400T,time)#get_q_ratio(path,time)#get_q_individ(path,time)
v_nx_k,v_ny_k,v_nx_c,v_ny_c = q_mod.get_Nernst_abs(path_400T,time)
x_grid = dict_qk['x_grid']*xstep_factor


#--------> debug
# print(np.shape(dict_qc['SH x']), np.shape(dict_qk['SH x']), np.shape(dict_qk['tot x']), np.shape(dict_qk['tot x']), np.shape(dict_qc['tot x']))
# print(np.shape(dict_qc['SH x']), np.shape(dict_qk['SH x']), np.shape(dict_qk['tot x']), np.shape(dict_qk['tot x']), np.shape(dict_qc['tot x']))
# print(np.shape(x_grid))
# print('v_nx = ', np.shape(v_nx_k),np.shape(v_nx_k),np.shape(v_ny_k),np.shape(v_nx_c),np.shape(v_ny_c))
# 
# print(' Bz = ', np.shape(dict_qk['Bz']), 'np.shape(Te = ', np.shape(dict_qk['Te']))
#<------------------ SH + RL--
Te_400 = dict_qk['Te'][:,iy]
Bz_400 = dict_qk['Bz'][:,iy]*Bz_factor
if Nernst_on:
    rat_400T = q_mod.get_ratio_lim(v_nx_c[:,iy],v_nx_k[:,iy],vmax=100.0,thresh=thresh_N_ratio)
    v_nx_c400 = v_nx_c[:,iy]
    v_nx_k400 = v_nx_k[:,iy]
else:
    rat_400T = q_mod.get_ratio_lim(dict_qc['SH x'][:,0],dict_qk['SH x'][:,iy],vmax=100.0)

##-------- Some multiplication factors to sort out axes normalisation
multT = 2.0



#----> plotting
print( ' shape ratio = ', np.shape(rat_50T),np.shape(dict_qc['SH x'][:,0]),np.shape(dict_qk['SH x'][:,iy]))
#fig = plt.figure()
fig = fprl.newfig_generic_twinx(1.0)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.tick_params(which='both',direction='in')
ax2.tick_params(which='both',direction='in')


vmin = np.min(rat_50T[np.abs(x_grid - xmin)<=xmax])*0.8
vmax = np.max(rat_50T[np.abs(x_grid - xmin)<=xmax])*1.2

#p0, = ax.plot(x_grid,rat_0T,c='g')
if Nernst_on:
    p50, = ax.plot(x_grid,v_nx_k50,c='r')
    p50, = ax.plot(x_grid,v_nx_c50,c='r',linestyle='--')

    p400, = ax.plot(x_grid,v_nx_k400,c='b')
    p400, = ax.plot(x_grid,v_nx_c400,c='b',linestyle='--')
else:
    p0, = ax.plot(x_grid,rat_0T,c='g')
    p50, = ax.plot(x_grid,rat_50T,c='r')
    p400, = ax.plot(x_grid,rat_400T,c='b')


#ax.plot(x_grid,dict_qc['SH x'],c='r',linestyle='--')
#ax.plot(x_grid,dict_qk['RL y'][:,iy],c='k')
#ax2.plot(x_grid,dict_qk['Bz'][:,iy],c='k',linestyle=':')


#pT0, = ax2.plot(x_grid,Te_0,c='k',linestyle='-.')
#pT50, = ax2.plot(x_grid,Te_50,c='k',linestyle='--')
#pT400, = ax2.plot(x_grid,Te_400,c='k',linestyle=':')
#pT0, = ax2.semilogy(x_grid,Bz_0,c='k',linestyle='-.')
if Nernst_on:

    pT50, = ax2.plot(x_grid,Bz_50,c='k',linestyle='--')
    pT400, = ax2.plot(x_grid,Bz_400,c='k',linestyle=':')
else:
    if plot_zero:
        pT0, = ax2.plot(x_grid,Te_0*multT,c='k',linestyle='-.')
    pT50, = ax2.plot(x_grid,Te_50*multT,c='k',linestyle='--')
    pT400, = ax2.plot(x_grid,Te_400*multT,c='k',linestyle=':')


#ax3 = ax.twinx()
#ax3.plot(x_grid,Bz_0,c='k',linestyle='-.')
#ax3.plot(x_grid,Bz_400,c='k',linestyle=':')
#ax3.plot(x_grid,Bz_50,c='k',linestyle='--')

#plist = [(p0,pT0),(p50,pT50),(p400,pT400)]
#plist = [('r','-.'),(p50,pT50),(p400,pT400)]
if plot_zero:
    leg_list = ['$\SI{0}{T}$','$\SI{50}{T}$','$\SI{400}{T}$'] 
    plt.legend([("g",'-.'),("r","--"),("b",":")], leg_list,
               handler_map={tuple: AnyObjectHandler()},
                frameon=False)


else:
    leg_list = ['$\SI{50}{T}$','$\SI{400}{T}$'] 
    plt.legend([("r","--"),("b",":")], leg_list,
               handler_map={tuple: AnyObjectHandler()},
               frameon=False)

#ax.legend(plist,leg_list,
#          handler_map={tuple: AnyObjectHandler()})


#ax2.plot(x_grid,dict_qk['Z2ni'][:,iy],c='gray',linestyle=':')
#ax2.plot(x_grid,dict_qc['RL y'],c='k',linestyle='--')
#ax2.plot(x_grid,v_nx_k,c='b')
#ax2.plot(x_grid,v_nx_c,c='b',linestyle='--')




Te_lab = r'$T_e$[\SI{%1.1f}{\kilo\electronvolt}]' % ((multT**-1)*2.0*T0*1e-3)
ax.set_xlabel(xlab)
if Nernst_on:
    ax.set_ylabel(r'$v_{x,N,k}, v_{x,N,c}$')
    ax2.set_ylabel(Bz_lab)
    ax.set_ylim(-0.1,0.01)#(vmin,vmax)
    ax2.set_ylim(0.0,1000.0)
    
else:
    ax.set_ylabel(r'$q_{x,\perp,k}/q_{x,\perp,c}$')
    ax2.set_ylabel(Te_lab)
    ax.set_ylim(0.0,2.0)
    
ax.set_xlim(xmin,xmax)
plt.savefig(save_name,dpi=600)
print( 'saving as: ', save_name)
print( ' copy and paste: open -a preview ' + save_name)
#ax2.plot(x_grid,dict_qk['SH x'][:,1],c='r',linestyle='--')
os.system('open -a preview ' + save_name)
#plt.show()
