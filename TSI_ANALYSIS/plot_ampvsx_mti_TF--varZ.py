'''
    16/3/19 - this script plots dT at different distances away from ablation surface
    
     + the Biermann + the div.q_RL - i.e. where dyqyrl intersect and have same sign? We have TS instability...
'''

import numpy as np, sys, os, getpass, site
sys.path.extend(["./"])

import matplotlib.pyplot as plt
import PLOTTERS.gen_Te_lineout as TEL
import gen_dyqy_RL_varZ as dyqy
import gen_dBdt_bier_varZ as gdBdt
import figure as fprl
import matplotlib.gridspec as GS
import matplotlib.ticker as ticker
from pylab import *
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
#-------------------------------------------------#
save_name = 'modeinv'
paths = hk.directory_paths()
src_dir = paths.src_dir    #'/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS'
data_dir = paths.data_dir_2D
save_path = '/OUTPUT/AMP/'
loc_nspace = hk.load_names()
path1 = loc_nspace.single_B
path1_noB = loc_nspace.single_noB
path1_speckle_tc5 = loc_nspace.speckle_B
path1_speckle_tc5_noB = loc_nspace.speckle_noB

save_path = loc_nspace.save_path
time_list = ['03', '05', '06']
time_glb = time_list[-1]
cwdpath = hk.directory_paths().norm_dir
cd5 = cf.conv_factors_custom(cwdpath)
cl_index = cd5.cl_index
c_index = cd5.c_index
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
dict_list = {}
# custom colours and lineout indices
dict_list['color_lineout'] = cd5.color_lineout    # custom colours for lines
dict_list['lineout_list'] = [110, 120, 130]    #[90,140, 160]##

if len(sys.argv) > 1:
    save_path = sys.argv[-1]
    if save_path[-1] != '/':
        save_path = save_path + '/'

    if not os.path.isdir(save_path):
        print 'NO DIRECTORY FOR SAVING IN WITH NAME: ', save_path
        sys.exit()
else:
    ##save_path = './pics'
    if save_path[-1] != '/':
        save_path = save_path + '/'
    print ' save_path = %s ' % save_path
save_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/RESULTS/'
if not os.path.exists(save_path):
    os.system('mkdir ' + save_path)

print ' SAVE_PATH = ', save_path
fname = 'Te_amp_premag_' + time_glb + '_dyqy'

savename = save_path + fname + 'v2'

#----- input data dirs ----
path0 = data_dir + 'LT1/r5_v40_TF_50T_dy01_2D1_vmax20'
path1 = data_dir + 'LT1/r5_v40_TF_50T_dy01_2D1_vmax20_100T'
path2 = data_dir + 'LT1/r5_v40_vmax20_vnonuni10_400T_Zb'

#path2 = 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_0T_E'
#path_list = [path1,path2]
path_list = [path0, path1, path2]
#leg_list = [r'$50\,\si{T}$']
leg_list = ['$50\,\si{T}$', '$100\,\si{T}$', '$400\,\si{T}$']

#-----------

plt.rcParams.update({
    'lines.linewidth': 1,
    'legend.fontsize': 8,
    'axes.titlesize': 12,
    'axes.linewidth': 1,
    'lines.linewidth': 1.0,
    'axes.labelsize': 12,
    'axes.labelpad': 1,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'ytick.minor.width': 0.7,
    'ytick.minor.size': 1,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12
})
#-----> new rc params

fig = fprl.newfig_generic(1.0)

gs1 = GS.GridSpec(3, 3, width_ratios=[1.4, 0.7, 0.04])
gs2 = GS.GridSpec(3, 3, width_ratios=[1.4, 0.7, 0.04])

#gs1.update(wspace=0.01, hspace=0.05)
gs1.update(wspace=0.5, hspace=0.07, left=0.15, right=0.92)
gs2.update(wspace=0.05, hspace=0.07, left=0.15, right=0.92)

ax_yTe = plt.subplot(gs1[:, 0])
ax_dyqy = plt.subplot(gs2[2, 1])
ax_Bz = plt.subplot(gs2[1, 1])
ax_xTe = plt.subplot(gs2[0, 1])
cax_Bz = plt.subplot(gs2[1, 2])
#cax_Bz = plt.subplot(gs2[1,2])

cax_dyqy = plt.subplot(gs2[2, 2])

if False:
    ax_yTe = fig.add_subplot(1, 2, 1)
    ax_xTe = fig.add_subplot(3, 2, 2)
    ax_Bz = fig.add_subplot(3, 2, 4)

    ax_dyqy = fig.add_subplot(3, 2, 6)
    ax_caxBz = fig.add_subplot(3, 3, 6)
    ax_caxdyqy = fig.add_subplot(3, 3, 9)
ax_list = [ax_yTe, ax_xTe, ax_Bz, ax_dyqy]
##plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})
style_list = ['-', '--', ':', '-', '--', ':']
mstyle_list = [None, None, None, 'x', '^', 'o']

p1 = TEL.plot_Te_xlineout(fig, ax_xTe, path_list, var_amp='Te', time=time_glb, dict_list=dict_list)
axwt = ax_xTe.twinx()
p1wt = TEL.plot_custom_xlineout(fig, axwt, path_list, var_amp='U', time=time_glb, lstyle='--')
#p1 = TEL.plot_Te_xlineout(fig,axwt,path_list,var_amp='Te',time=time_glb,dict_list=dict_list)

#p1wt = TEL.plot_custom_xlineout(fig,axwt,path_list,var_amp='wt',time=time_glb,lstyle=':')
#p1wt = TEL.plot_custom_xlineout_amp(fig,axwt,path_list,var_amp='U',time=time_glb,
#                                    style_list=style_list,mstyle_list=mstyle_list)

p1wt = TEL.plot_custom_xlineout_amp_tevol(fig,
                                          ax_dyqy,
                                          path_list,
                                          var_amp='Te',
                                          time_list=time_list,
                                          style_list=style_list,
                                          mstyle_list=mstyle_list,
                                          leg_dict=dict_list,
                                          axleg=cax_dyqy)
cax_dyqy.get_xaxis().set_visible(False)
cax_dyqy.axes.get_yaxis().set_visible(False)
for item in [cax_dyqy]:
    item.patch.set_visible(False)
cax_dyqy.axis('off')

p2 = TEL.plot_ylineout_custom(fig,
                              ax_yTe,
                              path_list,
                              var_amp='Te',
                              time=time_list[-1],
                              style_list=style_list,
                              mstyle_list=mstyle_list,
                              dict_list=dict_list)

mstyle_list = ['x', 'x', 'o', None, None, None]
if False:
    p2 = TEL.plot_ylineout_custom(fig,
                                  ax_yTe,
                                  path_list,
                                  var_amp='Te',
                                  time=time_list[1],
                                  style_list=style_list,
                                  mstyle_list=mstyle_list,
                                  dict_list=dict_list)

#ax_xTe.ticklabel_format(axis='x','%1.1f')

#ax_xTe.yaxis.set_major_locator(MaxNLocator(4,prune='upper '))#sets the max number of ticks on yaxis
ax2_list = [ax_Bz, ax_dyqy, axwt]
for ax in ax2_list:
    ax.yaxis.set_major_locator(MaxNLocator(4,
                                           prune='both'))    #sets the max number of ticks on yaxis

ax_yTe.yaxis.set_major_locator(MaxNLocator(4))    #sets the max number of ticks on yaxis
ax_yTe.xaxis.set_major_locator(MaxNLocator(3))    #sets the max number of ticks on yaxis

#ax_yTe.set_ylim(-0.1,0.1)
yte_l, yte_u = ax_xTe.get_ylim()
#ax_xTe.set_ylim(yte_l,1.2)
#ax_yTe.set_ylim(-2,2)

# plot Bz ---
#--> disable plot dyqy
#path,fprefix,var,time = path1,path1,'Bz',time_glb
#im  = Bz_im.get_var_imshow_custom(fig,ax_Bz,cax_Bz,path,fprefix,var,time,xi_max = 354)
path, time = path1, time_glb
#---> enable plot dyqy of other
#im = gdBdt.plot_dBdt_bier(fig,ax_Bz,cax_Bz,path,time)
im = dyqy.plot_dyqy_RL(fig, ax_Bz, cax_Bz, path, time)

ax_xTe.set_xticklabels([])
ax_Bz.set_xticklabels([])
ax_xTe.set_xlabel('')
ax_Bz.set_xlabel('')
#------>
ax_list = [ax_Bz, ax_xTe, ax_dyqy, cax_Bz, ax_yTe]

#<-------

fprefix = path.split('/')[-1]
dict_temp = cf.load_dict(path, path.split('/')[-1], 'Te', time_glb)
x_grid, y_grid = dict_temp['x_grid'], dict_temp['y_grid']
x_grid_SI = x_grid[cl_index:c_index] * xstep_factor

print 'limits = ', x_grid_SI[0], x_grid_SI[-1]
xax_min, xax_max = ax_Bz.get_xlim()
xax_min, xax_max = -5.0, 20.0

#xax_min,xax_max = -20.0,x_grid_SI[-1]

yax_min, yax_max = ax_Bz.get_ylim()

#xBz_min,xBz_max = ax_xTe.get_xlim()
ax_xTe.set_xlim(xax_min, xax_max)
axwt.set_xlim(xax_min, xax_max)
ax_Bz.set_xlim(xax_min, xax_max)
ax_dyqy.set_xlim(xax_min, xax_max)
ax_yTe.set_xlim(yax_min, yax_max)

multx, multy = 1.2, 1.2
fsize = 12
fprl.annotate_axis(ax_yTe, lett=r'$(a)$', dx_mult=1.4, dy_mult=0.0, fontsize=fsize)
fprl.annotate_axis(ax_xTe, lett=r'$(b)$', dx_mult=multx, dy_mult=0.0, fontsize=fsize)
fprl.annotate_axis(ax_Bz, lett=r'$(c)$', dx_mult=multx, dy_mult=multy, fontsize=fsize)
fprl.annotate_axis(ax_dyqy, lett=r'$(d)$', dx_mult=multx, dy_mult=multy, fontsize=fsize)

#ax_xTe.get_shared_x_axes().join(ax_dyqy, ax_Bz,ax_xTe)
#ax_Bz.get_shared_x_axes().join(ax_dyqy, ax_Bz,ax_xTe)
#ax_dyqy.get_shared_x_axes().join(ax_dyqy, ax_Bz,ax_xTe)
#ax_Bz.get_shared_x_axes().join(ax_xTe)
#ax_dyqy.get_shared_x_axes().join(ax_xTe)

for ax in ax_list:

    ax.tick_params(which='both', direction='in')
fprl.savefig(savename)
print ' saving as --- ', savename
os.system('open -a preview ' + savename + '.pdf')
#plt.show()
