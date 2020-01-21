
'''
27/08/2019
- Rehash of plot_q_lineout.py


'''
import numpy as np, sys, os, getpass, site, re
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/MODULES')
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/PLOTTERS')
import matplotlib.pyplot as plt
import figure_prl_twocol as fprl
import matplotlib.gridspec as GS
import matplotlib.ticker as ticker
from pylab import *
import kinetic_ohmslaw_module_varZ as q_mod
import gen_Te_lineout as TEL
import chfoil_module as cf
import house_keeping as hk
import tsi_module as tsi



#---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
#-----
# functions
fpre = lambda path_in: path_in.split('/')[-1]
def b_lab(bz_in):
    if int(bz_in) == -1:
        lab = 'no B'
    else:
        lab = r'$%i\,\si{T}$' % (bz_in)
    return lab
#---> file inputs
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir
norm_dir = paths.norm_dir
log_file = norm_dir + 'norm.log'
[T_ref,n_ref,Z_ref,Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir,Z_ref,Ar=6.51)


Te_factor = 2.0*T_ref*1e-3
#--->  getting path ->
lambda_p = 5
time_list = ['03','05','06']

path_list = []

#--- generating  path_list
scale_len = 1
save_tag = 'LT%i_' % (scale_len)
bz_list = [-1,50.0,100.0]
lab_list = []
for ip in range(len(bz_list)):    
    path_loc = paths.get_path(scale_len,bz_list[ip],lambda_p)
    path_list.append(path_loc)
    lab_list.append(b_lab(bz_list[ip]))

#aesthetics
cmap = cm.viridis
style_list = [':','--','-','-','--',':']
mstyle_list = [None,None,None,'x','^','o']


#<<<--- transport inputs
time = '07'
iy = 20
xmin,xmax = -2.0,20.0
thresh_N_ratio = 5e-2
var = 'tot y'#'RL y'#sys.argv[-1]
save_name = '%sdT_%s.png' % (save_path,save_tag)

fig = fprl.newfig_generic_twinx(1.0,scale_width=1.0,scale_ratio=1.0)
ax_t = fig.add_subplot(111)
ax_Te = ax_t.twinx()

#------------------->>>

p1dt = TEL.plot_custom_xlineout_amp_tevol(fig,ax_t,path_list,var_amp='Te',time_list=time_list,
                                    style_list=style_list,mstyle_list=mstyle_list,
                                    cmap=['r','g','b','k'],
                                    axleg=ax_t,
                                    leg_dict={'title':cd5.tlab,'some_other thing':1.0})
#leg = ax_t.legend(p1dt, leg_list,loc = 'bottom right')

var = 'Te'
p_list = []
for pp in range(len(path_list)):
    time = time_list[-1]
    dict_T = cf.load_dict(path_list[pp],fpre(path_list[pp]),var,time)
    
    p2, = ax_Te.plot(dict_T['x_grid']*cd5.xstep_factor,dict_T['mat'][:,iy]*(2.0*T_ref/1e3),
               c='k',linestyle=style_list[pp])
    p_list.append(p2)
ax_Te.set_ylabel(r'$T_e\,[\si{keV}]$')
# legend
leg2 = ax_Te.legend(p_list,lab_list)
leg2.get_frame().set_linewidth(0.0)
leg2.get_frame().set_facecolor('none')

ax_t.set_xlim(xmin,xmax)
ax_Te.tick_params(axis='both',which='both',direction='in')
ax_t.tick_params(axis='both',which='both',direction='in')





#plt.show()

plt.savefig(save_name,dpi=600)
print( 'saving as: ', save_name)
print( ' copy and paste: open -a preview ' + save_name)

#plt.show()
#plt.close()
