'''
    16/3/19 - this script plots dT at different distances away from ablation surface
    
     + the Biermann + the div.q_RL - i.e. where dyqyrl intersect and have same sign? We have TS instability...
'''

import numpy as np, sys, os, getpass, site
sys.path.extend(["./"])

import matplotlib.pyplot as plt
import PLOTTERS.gen_Te_lineout as TEL
import PLOTTERS.gen_dyqy_RL_varZ as dyqy
import PLOTTERS.gen_dBdt_bier_varZ as gdBdt
import MODULES.figure as fprl
import matplotlib.gridspec as GS
from pylab import *
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
#-------------------------------------------------#
# load in initial modules
save_name = 'modeinv'
paths = hk.directory_paths()
src_dir = paths.src_dir    #'/Users/dominichill/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS'
data_dir = paths.data_dir_2D
save_path = paths.save_dir

#--- User input --->------------>>>>>
scale_length = 2
lambda_p = 5    #perturbation length in mfp, can only be 5 or 10
bz_in = [-1, 50, 400]
bz_in = [50, 100, 400]
xax_min, xax_max = -5.0, 40.0

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
def get_tmax(path_list):
    t_out = []

    for path in path_list:
        t_list, tc_list = cf.get_t_list(path, var='Te')
        print(' t_list = ', t_list)
        t_out.append(int(t_list[-1]))
    print('t_out = ', t_out)

    return '%02i' % np.min(t_out)


time_glb = get_tmax(path_list)
print('time_glb = ', time_glb)
time_list = ['02', time_glb]
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
fname = 'BZ_amp_premag_' + time_glb + str(scale_length) + sys.argv[0].split('.')[0]

savename = save_path + fname

#----- input data dirs ----
#path0 = data_dir + 'LT1/r5_v40_TF_50T_dy01_2D1_vmax20'
#path1 = data_dir + 'LT1/r5_v40_TF_50T_dy01_2D1_vmax20_100T'
#path2 = data_dir +  'LT1/r5_v40_vmax20_vnonuni10_400T_Zb'
#path2 = data_dir + 'LT1/r5_v40_TF_50T_dy01_2D1_vmax20_100T'

#leg_list = ['$50\,\si{T}$','$100\,\si{T}$','$400\,\si{T}$']
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

#gs1.update(wspace=0.01, hspace=0.05)
gs1.update(wspace=0.5, hspace=0.07, left=0.1, right=0.9)
gs2.update(wspace=0.05, hspace=0.07, left=0.1, right=0.9)

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

style_list = [':', '--', '-', '-', '--', ':']
mstyle_list = [None, None, None, 'x', '^', 'o']

p1 = TEL.plot_Te_xlineout(fig, ax_xTe, path_list, var_amp='Te', time=time_glb, dict_list=dict_list)
axwt = ax_xTe.twinx()
p1wt = TEL.plot_custom_xlineout(fig, axwt, path_list, var_amp='U', time=time_glb, lstyle='--')
axwt.grid('off')
ax_xTe.grid('off')

shade_of_grey = lambda alpha: [alpha, alpha, alpha]
p1wt = TEL.plot_custom_xlineout_amp_tevol(fig,
                                          ax_dyqy,
                                          path_list,
                                          var_amp='Bz',
                                          time_list=time_list,
                                          style_list=style_list,
                                          mstyle_list=mstyle_list,
                                          leg_dict=dict_list,
                                          axleg=cax_dyqy,
                                          cmap=[shade_of_grey(0.5),
                                                shade_of_grey(0.0)],
                                          log_on=True)
ax_dyqy.set_ylim(1e-7, 0.1)

cax_dyqy.get_xaxis().set_visible(False)
cax_dyqy.axes.get_yaxis().set_visible(False)
for item in [cax_dyqy]:
    item.patch.set_visible(False)
cax_dyqy.axis('off')
cd5.xmin, cd5.xmax = xax_min, xax_max
p2 = TEL.plot_ylineout_custom(fig,
                              ax_yTe,
                              path_list,
                              var_amp='Bz',
                              time=time_list[-1],
                              style_list=style_list,
                              mstyle_list=mstyle_list,
                              dict_list=dict_list,
                              cfg=cd5)

# plot Bz ---
#--> disable plot dyqy
path, time = path_list[-1], time_glb
#---> enable plot dyqy of other
im = gdBdt.plot_dBdt_bier(fig, ax_Bz, cax_Bz, path, time, xlim=[-5.0, 40.0])
#im  = dyqy.plot_dyqy_RL(fig,ax_Bz,cax_Bz,path,time)

#ax_xTe.yaxis.set_major_locator(MaxNLocator(4,prune='upper '))#sets the max number of ticks on yaxis
ax2_list = [ax_Bz, axwt]
for ax in ax2_list:
    ax.yaxis.set_major_locator(MaxNLocator(4,
                                           prune='both'))    #sets the max number of ticks on yaxis

#ax_dyqy.yaxis.set_major_locator(MaxNLocator(2,prune='upper'))#sets the max number of ticks on yaxis
ax_dyqy.yaxis.set_ticks([1e-7, 1e-4])

ax_yTe.yaxis.set_major_locator(MaxNLocator(4))    #sets the max number of ticks on yaxis
ax_yTe.xaxis.set_major_locator(MaxNLocator(3))    #sets the max number of ticks on yaxis

yte_l, yte_u = ax_xTe.get_ylim()

# ax tick labels formatting
ax_xTe.set_xticklabels([])
ax_Bz.set_xticklabels([])
ax_xTe.set_xlabel('')
ax_Bz.set_xlabel('')

ax_yTe.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
ax_Bz.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))

#------>
ax_list = [ax_Bz, ax_xTe, ax_dyqy, cax_Bz, ax_yTe]

#<-------

fprefix = path.split('/')[-1]
dict_temp = cf.load_dict(path, path.split('/')[-1], 'Te', time_glb)
x_grid, y_grid = dict_temp['x_grid'], dict_temp['y_grid']
x_grid_SI = x_grid[cl_index:c_index] * xstep_factor

yax_min, yax_max = ax_Bz.get_ylim()

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

for ax in ax_list:

    ax.tick_params(which='both', direction='in')
fprl.savefig(savename)
print ' saving as --- ', savename + '.pdf'
os.system('open -a preview ' + savename + '.pdf')
#plt.show()
