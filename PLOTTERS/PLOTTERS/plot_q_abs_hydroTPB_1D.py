'''
test plot
'''

import numpy as np
import matplotlib.pyplot as plt
import sys, re, os, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import chfoil_module as cf
#import q_SH_Te_tsteps as q_mod
import matplotlib as mpl
from pylab import *
import figure as fprl
import kinetic_ohmslaw_module_1D_varZ as q_mod
from matplotlib.legend_handler import HandlerBase

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12

norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + 'norm.log'
[T0, n0, Z0, Bz0] = np.loadtxt(log_file)

cd5 = cf.conv_factors_custom(norm_path, Z0, Ar=6.51)
print 'Z0 = ', Z0
print 'T0 = ', T0
print 'n0 = ', n0
print ' Bz0 = ', Bz0

cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index, c_index = 0, -1
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab_rel
ylab = cd5.ylab
Bz0 = (m_e / (q_e * tau_ei))
Bz_lab = r'$B_z$[\si{T}]'
Bz_factor = Bz0
P_factor = 2.0 / 3.0
P_lab = r'$P$'
#<--------------
# inputs

#r5_v40_Z_FEOS_trackncrit_D - 400T
#r5_v40_Z_FEOS_trackncrit_50T_FF - 50T
#r5_v40_Z_FEOS_trackncrit_1D_0T - 0T
#-------------------------------------------->
path_pre = '../PREMAG/'
path_0T = path_pre + 'r5_v40_Z_FEOS_trackncrit_1D_0T'
path_50T = path_pre + 'r5_v40_Z_FEOS_trackncrit_50T_FF'
path_400T = path_pre + 'r5_v40_Z_FEOS_trackncrit_D'
Nernst_on = True
time = '13'
iy = 1
xmin, xmax = -5.0, 20.0
thresh_N_ratio = 5e-2
var = 'vN x'
var_list = ['SH x', 'SH y', 'RL x', 'RL y', 'E x', 'E y', 'vN x', 'vN y']
lab_list = [
    r'$q_{\perp,x}$', r'$q_{\perp,y}$', r'$q_{RL,x}$', r'$q_{RL,y}$', r'$q_{E,x}$', r'$q_{E,y}$',
    r'$q_{E,x}$', r'$v_{N,x}$', r'$v_{N,y}$'
]
ylab = lab_list[var_list.index(var)]


#--------------------------------------------->
#--- functions
class AnyObjectHandler(HandlerBase):

    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.7 * height, 0.7 * height],
                        linestyle=orig_handle[1],
                        color='k')
        l2 = plt.Line2D([x0, y0 + width], [0.3 * height, 0.3 * height], color=orig_handle[0])
        return [l1, l2]


#-----------------------------------------------------------------------
def repack(path, time):
    '''
        dict_c,dict_k = repack(path,time)
    '''
    dict_qc, dict_qk = q_mod.get_q_abs(path,
                                       time)    #get_q_ratio(path,time)#get_q_individ(path,time)
    v_nx_k, v_ny_k, v_nx_c, v_ny_c = q_mod.get_Nernst_abs(path, time)

    dict_out_c = {}
    dict_out_k = {}
    for var in var_list:
        if var[0] != 'v':
            dict_out_k[var] = dict_qk[var][:, 0]
        elif var == 'vN x':
            dict_out_k[var] = v_nx_k[:, iy]
        elif var == 'vN y':
            dict_out_k[var] = v_ny_k[:, iy]

    for var in var_list:
        if var[0] != 'v':
            dict_out_c[var] = dict_qc[var][:, 0]
        elif var == 'vN x':
            dict_out_c[var] = v_nx_c[:, iy]
        elif var == 'vN y':
            dict_out_c[var] = v_ny_c[:, iy]

    dict_out_k['U'] = dict_qk['U'][:, iy] * P_factor
    dict_out_k['Te'] = dict_qk['Te'][:, iy]
    dict_out_k['Bz'] = dict_qk['Bz'][:, iy] * Bz_factor
    dict_out_k['x_grid'] = dict_qk['x_grid'] * xstep_factor

    return dict_out_c, dict_out_k


#----> data loading
dict_0T_c, dict_0T_k = repack(path_0T, time)
dict_50T_c, dict_50T_k = repack(path_50T, time)
dict_400T_c, dict_400T_k = repack(path_400T, time)

x_grid = dict_0T_k['x_grid']
#q_mod.get_ratio_lim(dict_qc['SH x'][:,0],dict_qk['SH x'][:,iy],vmax=100.0)
#<--------------- 0T loading---
#----> plotting
#print( ' shape ratio = ', np.shape(rat_50T),np.shape(dict_qc['SH x'][:,0]),np.shape(dict_qk['SH x'][:,iy]))
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax3 = ax.twinx()

data_k = dict_50T_k[var]

vmin = np.min(data_k[np.abs(x_grid - xmin) <= xmax]) * 1.8
vmax = np.max(data_k[np.abs(x_grid - xmin) <= xmax]) * 1.2

#p0, = ax.plot(x_grid,rat_0T,c='g')

p50, = ax.plot(x_grid, dict_50T_k[var], c='r')
p50, = ax.plot(x_grid, dict_50T_c[var], c='r', linestyle='--')

p400, = ax.plot(x_grid, dict_400T_k[var], c='b')
p400, = ax.plot(x_grid, dict_400T_c[var], c='b', linestyle='--')

#ax.plot(x_grid,dict_qc['SH x'],c='r',linestyle='--')
#ax.plot(x_grid,dict_qk['RL y'][:,iy],c='k')
#ax2.plot(x_grid,dict_qk['Bz'][:,iy],c='k',linestyle=':')

#pT0, = ax2.plot(x_grid,Te_0,c='k',linestyle='-.')
#pT50, = ax2.plot(x_grid,Te_50,c='k',linestyle='--')
#pT400, = ax2.plot(x_grid,Te_400,c='k',linestyle=':')
#pT0, = ax2.semilogy(x_grid,Bz_0,c='k',linestyle='-.')
pT50, = ax2.plot(x_grid, dict_50T_k['Bz'], c='k', linestyle='--')
pT400, = ax2.plot(x_grid, dict_400T_k['Bz'], c='k', linestyle=':')
pT50, = ax3.plot(x_grid, dict_50T_k['U'], c='gray', linestyle='--')
pT400, = ax3.plot(x_grid, dict_400T_k['U'], c='gray', linestyle=':')
pT50, = ax3.plot(x_grid, dict_50T_k['Te'], c='b', linestyle='--')
pT400, = ax3.plot(x_grid, dict_400T_k['Te'], c='b', linestyle=':')

leg_list = ['$\SI{50}{T}$', '$\SI{400}{T}$']
plt.legend([("r", "--"), ("b", ":")], leg_list, handler_map={tuple: AnyObjectHandler()})

Te_lab = r'$T_e$[\SI{%1.1f}{\kilo\electronvolt}]' % (2.0 * T0 * 1e-3)
ax.set_xlabel(xlab)
ax.set_ylabel(ylab)
ax2.set_ylabel(Bz_lab)

#ax.set_ylim(-0.1,0.01)#(vmin,vmax)
ax.set_ylim(vmin, vmax)
ax2.set_ylim(0.0, 1000.0)
ax.set_xlim(xmin, xmax)

#ax2.plot(x_grid,dict_qk['SH x'][:,1],c='r',linestyle='--')
plt.show()
