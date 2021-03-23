'''
    16/3/19 - this script plots dT at different distances away from ablation surface
    
     + the Biermann + the div.q_RL - i.e. where dyqyrl intersect and have same sign? We have TS instability...
'''

import sys
import os
import numpy as np
from pylab import *
import matplotlib.gridspec as GS

sys.path.extend(["./"])
import PLOTTERS.gen_Te_lineout as TEL
import PLOTTERS.gen_dyqy_RL_varZ as dyqy
import MODULES.figure as fprl
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
#-------------------------------------------------#
# load in initial modules
save_name = 'modeinv'
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir

#--- User input --->------------>>>>>
scale_length = 1
lambda_p = 5    #perturbation length in mfp, can only be 5 or 10
bz_in = [-1, 50, 400]
bz_in = [0, 400]

#-- initialise zero size arrays
path_list = []
leg_list = []


def get_bz_lab(bz):
    if int(bz) == -1:
        lab = 'no B'
    elif int(bz) == -2:
        lab = r'$100\% \,\, 50\,\si{T}$'
    else:
        lab = '$%i\,\si{T}$' % (bz)
    return lab


for ip in range(len(bz_in)):

    path_loc = paths.get_path(scale_length, bz_in[ip], lambda_p)
    path_list.append(path_loc)
    leg_list.append(get_bz_lab(bz_in[ip]))

print('--> path_list = ', path_list)
print('--> leg_list = ', leg_list)
#<-------------------------------------------
t_list, tc_list = cf.get_t_list(path_list[1], var='Te')
time = t_list[-1]
time_list = ['03', t_list[-1]]
#time_list = ['03','15']

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
dict_list['lineout_list'] = [100, 120, 140]    #[64,67,70]#[110,120,130]#[90,140, 160]##

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
if not os.path.exists(save_path):
    os.system('mkdir ' + save_path)

print ' SAVE_PATH = ', save_path
fname = 'Te_amp_premag_' + time_glb + '_dyqy'

savename = save_path + fname + 'v2'

#----- input data dirs ----

cd5.leg_list = leg_list
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

gs1.update(wspace=0.5, hspace=0.07, left=0.14, right=0.91)
gs2.update(wspace=0.05, hspace=0.07, left=0.14, right=0.91)

ax_a = plt.subplot(gs1[:, 0])
ax_d = plt.subplot(gs2[2, 1])
ax_c = plt.subplot(gs2[1, 1])
ax_b = plt.subplot(gs2[0, 1])
cax_c = plt.subplot(gs2[1, 2])
cax_d = plt.subplot(gs2[2, 2])

ax_list = [ax_a, ax_b, ax_c, ax_d]
##plt.rcParams.update({'font.sans-serif': 'Arial', 'font.family': 'sans-serif'})
style_list = ['--', '-', ':', '-', '--', ':']
mstyle_list = [None, None, None, 'x', '^', 'o']

p1 = TEL.plot_Te_xlineout(fig, ax_b, path_list, var_amp='Te', time=time_glb, dict_list=dict_list)
ax_b2 = ax_b.twinx()
p1wt = TEL.plot_custom_xlineout(fig, ax_b2, path_list, var_amp='wt', time=time_glb, lstyle='--')

shade_of_grey = lambda alpha: [alpha, alpha, alpha]

p1wt = TEL.plot_custom_xlineout_amp_tevol(fig,
                                          ax_d,
                                          path_list,
                                          var_amp='Te',
                                          time_list=time_list,
                                          style_list=style_list,
                                          mstyle_list=mstyle_list,
                                          cmap=[shade_of_grey(0.5),
                                                shade_of_grey(0.0)],
                                          leg_dict=dict_list,
                                          axleg=cax_d)
cax_d.get_xaxis().set_visible(False)
cax_d.axes.get_yaxis().set_visible(False)
for item in [cax_d]:
    item.patch.set_visible(False)
cax_d.axis('off')

p2 = TEL.plot_ylineout_custom(fig,
                              ax_a,
                              path_list,
                              var_amp='Te',
                              time=time_list[-1],
                              style_list=style_list,
                              mstyle_list=mstyle_list,
                              dict_list=dict_list,
                              cfg=cd5)

mstyle_list = ['x', 'x', 'o', None, None, None]

#ax_b.yaxis.set_major_locator(MaxNLocator(4,prune='upper '))#sets the max number of ticks on yaxis
ax2_list = [ax_c, ax_d, ax_b2]
for ax in ax2_list:
    ax.yaxis.set_major_locator(MaxNLocator(4,
                                           prune='both'))    #sets the max number of ticks on yaxis

ax_a.yaxis.set_major_locator(MaxNLocator(4))    #sets the max number of ticks on yaxis
ax_a.xaxis.set_major_locator(MaxNLocator(3))    #sets the max number of ticks on yaxis

#ax_a.set_ylim(-0.1,0.1)
yte_l, yte_u = ax_b.get_ylim()

# plot Bz ---
#--> disable plot dyqy
path, time = path_list[1], time_glb
#---> enable plot dyqy of other
#im = gdBdt.plot_dBdt_bier(fig,ax_c,cax_c,path,time)
im = dyqy.plot_dyqy_RL_c(fig, ax_c, cax_c, path, time, switch_kinetic_on=False)
#plot_dyqy_RL_c(fig,ax,cax,path,time,cfg=cfg,**kwargs)

ax_b.set_xticklabels([])
ax_c.set_xticklabels([])
ax_b.set_xlabel('')
ax_c.set_xlabel('')
#------>
ax_list = [ax_c, ax_b, ax_d, cax_c, ax_a]

#<-------

fprefix = path.split('/')[-1]
dict_temp = cf.load_dict(path, path.split('/')[-1], 'Te', time_glb)
x_grid, y_grid = dict_temp['x_grid'], dict_temp['y_grid']
x_grid_SI = x_grid[cl_index:c_index] * xstep_factor

xax_min, xax_max = ax_c.get_xlim()
xax_min, xax_max = -5.0, 20.0

yax_min, yax_max = ax_c.get_ylim()

ax_b.set_xlim(xax_min, xax_max)
ax_b2.set_xlim(xax_min, xax_max)
ax_c.set_xlim(xax_min, xax_max)
ax_d.set_xlim(xax_min, xax_max)
ax_a.set_xlim(yax_min, yax_max)

multx, multy = 1.2, 1.2
fsize = 12
fprl.annotate_axis(ax_a, lett=r'$(a)$', dx_mult=1.4, dy_mult=0.0, fontsize=fsize)
fprl.annotate_axis(ax_b, lett=r'$(b)$', dx_mult=multx, dy_mult=0.0, fontsize=fsize)
fprl.annotate_axis(ax_c, lett=r'$(c)$', dx_mult=multx, dy_mult=multy, fontsize=fsize)
fprl.annotate_axis(ax_d, lett=r'$(d)$', dx_mult=multx, dy_mult=multy, fontsize=fsize)

for ax in ax_list:

    ax.tick_params(which='both', direction='in')
fprl.savefig(savename)
print ' saving as --- ', savename
os.system('open -a preview ' + savename + '.pdf')
#plt.show()
