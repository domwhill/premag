
'''
27/08/2019
- Rehash of plot_q_lineout.py
 Plots Phase shift of kinetic/classical heat flows

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
import pdb

import chfoil_module as cf
import house_keeping as hk
import tsi_module as tsi



#---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
#-----



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
time_list = ['03','05','07']
path_list = []
bz_in = 400.0# selected magnetic field [Tesla]

#--- generating  path_list
s_list = [2]
save_tag = '%s_%iT_' % (sys.argv[0].split('.')[0],bz_in)
for ip in range(len(s_list)):    
    path_loc = paths.get_path(s_list[ip],bz_in,lambda_p)
    path_list.append(path_loc)
    save_tag = save_tag + str(s_list[ip])

#aesthetics
xmin,xmax = -10,30.0 # max/min of x axis on plot
cmap = cm.viridis
style_list = ['--','--',':','-','--',':']
mstyle_list = [None,None,None,'x','^','o']


#<<<--- transport inputs
time = '07'
iy = 20
xmin,xmax = 0.0,40.0
thresh_N_ratio = 5e-2
if len(sys.argv)>1:
    var = sys.argv[1]
else:
    var = 'tot y'
if len(sys.argv)>2:
    save_path = save_path  + sys.argv[-1] + '/'
    
save_name = '%sdeltaqy_phase_%s_%s_%s' % (save_path,save_tag,re.sub(' ','', var),time)
phase_tag = '_phase_'


plot_amp_var_on = True
plot_div_var_on = False # this should overide above to take the divergence of the heat flow.
var_list = ['tot x','tot y',
            'SH x','SH y',
            'RL x','RL y',
            'E x','E y',
            'vN x','vN y']

lab_list = [r'$q_{x}$',    r'$q_{y}$',
            r'$q_{\perp,x}$',r'$q_{\perp,y}$',
            r'$q_{RL,x}$',r'$q_{RL,y}$',
            r'$q_{E,x}$', r'$q_{E,y}$', 
            r'$v_{N,x}$',r'$v_{N,y}$']

var_plot_list = ['tot x','tot y'
                'SH x','SH y',
                'RL x','RL y',
                'vN x','vN y']
var_dict = {'tot x':r'$q_{x}$', 
            'tot y': r'$q_{y}$',
            'SH x': r'$q_{\perp,x}$',
            'SH y': r'$q_{\perp,y}$',
            'RL x': r'$q_{RL,x}$',
            'RL y': r'$q_{RL,y}$',
            'E x': r'$q_{E,x}$', 
            'E y': r'$q_{E,y}$', 
            'vN x': r'$v_{N,x}$',
            'vN ,': r'$v_{N,y}$'}

ylab = lab_list[var_list.index(var)]

add_str = ''
if plot_amp_var_on:
    add_str = add_str + r'\delta ' 
if plot_div_var_on:
    add_str = r'-\partial_y' + add_str

ylab = r'$' + add_str + ylab[1:]

multQ = 1e-18
qrl_norm = (cd5.n0*1e6)*m_e*(cd5.v_te**3)*multQ
q_unit = r'$\,[\SI{%1.1e}{W/m^2}]$' % (multQ**-1) # kg s^-3 => kg m^-1 s^-2   ###nn    
ylab = ylab +  q_unit                   # kg m^2 s^{-2}/(m^2 s) = W/m^2
xlab = r'y [$\si{\micro\meter}$]'


def get_amp(data):

    return data - np.average(data)


def get_amp1(grid, data):
    
    if plot_amp_var_on:
        data_out =  data - np.average(data)

    elif plot_div_var_on:
        data_out = np.zeros((np.shape(data)))
        dy = (grid[2]-grid[0])
        data_out[1:-1] = -(data[2:] - data[:-2])/(dy)
        data_out[0] = data_out[1]
        data_out[-1] = data_out[-2]
    
    else:
        data_out = data.copy()   
    
    return data_out


def wavelength(grid, data):
    l2 = grid[np.where(data == np.max(data))] \
         - grid[np.where(data == np.min(data))]
    return l2*2.0


def get_phase(grid, data):
    l = wavelength(grid,data)
    max_val  = grid[np.where(data==np.max(data))]
    min_val = grid[np.where(data == np.min(data))]
    frac = max_val/l
    frac2 = min_val/l
    phase = frac*(2.0*np.pi)
    phase2 = np.max([frac2*(2.0*np.pi) - np.pi, frac2*(2.0*np.pi) + np.pi])

    return phase[0]


def get_phase_array(grid,data):
    """
        phase_arr = get_phase_array(grid,data,axis=0)
    """
    
    nx = len(data.take(0,axis=1))
    phase_arr = np.zeros((nx))
    for ix in range(nx):
        data_amp = get_amp1(grid, data.take(ix,axis=0))
        phase_arr[ix] = get_phase(
                            grid, data_amp
                            )
    return phase_arr


def conv_phase(phase):
    return (phase - np.pi)/(2.0*np.pi)
    
slice_array_y = lambda array: array[iy,:]

fig1 = fprl.newfig_generic_2yscale(1.0)#(1.1,scale_width=1.5,scale_ratio=0.5)#plt.figure()
ax2 = fig1.add_subplot(111)
ax5 = ax2.twinx()




#------------------->>>
for ipp in range(len(path_list)):
    pp = path_list[ipp]
    save_name_npy_c = pp.split('/')[-1] + phase_tag +  '_dict_c.npy'
    save_name_npy_k = pp.split('/')[-1] +  phase_tag +  '_dict_k.npy'

    if os.path.isfile(save_name_npy_c) and os.path.isfile(save_name_npy_k):
        dict_c = np.load(save_name_npy_c,allow_pickle=True).any()
        dict_k = np.load(save_name_npy_k,allow_pickle=True).any()
    else:
        dict_c,dict_k = q_mod.repack_2D(pp,time,recompute_q_c = True)
    x_grid = dict_k['x_grid']
    y_grid = dict_k['y_grid']
    data_k = dict_k[var]
    ipp = path_list.index(pp)
    
    # te_amp
    ax_amp = ax2.twinx()
    nx,ny = np.shape(dict_k['Te'])
    pT50, = ax2.plot(x_grid*cd5.xstep_factor,dict_k['Te'][:,ny//2]*Te_factor,c='k',
                    linestyle='-')

    te_amp  = cf.get_U_dev_abs(dict_k['Te'])
    pT_amp, = ax_amp.semilogy(
                        x_grid*cd5.xstep_factor,np.max(np.abs(te_amp*Te_factor), axis=1),
                        c='k',
                        linestyle='-'
                    )
    
    #--->
    var = 'tot y'
    color = 'b'
    c_list = ['b','r','g']
    var_list = ['tot y', 'SH y', 'RL y']
    p_list = []
    lab_list = []
    for ic, var in enumerate(var_list):
        color = c_list[ic]
        phase_arr_c = get_phase_array(y_grid,dict_c[var])
        phase_arr_k = get_phase_array(y_grid,dict_k[var])
        
        p1, = ax5.plot(x_grid*cd5.xstep_factor, conv_phase(phase_arr_k),
                c=color,
                linestyle='-'
                )
        ax5.plot(
                x_grid*cd5.xstep_factor, conv_phase(phase_arr_c),
                c=color,
                linestyle='--'
                )
        p_list.append(p1)
        lab_list.append(var_dict[var])
    if not os.path.isfile(save_name_npy_c) and not os.path.isfile(save_name_npy_k):
        np.save(save_name_npy_c, dict_c, allow_pickle=True)
        np.save(save_name_npy_k, dict_k, allow_pickle=True)


    ax2.set_ylabel(r'$T_e$ [eV]')
    #----->
    ax5.set_xlim(xmin,xmax)    
    ax_amp.set_xlim(xmin,xmax)
    ax2.set_xlim(xmin,xmax)    
    ax2.ticklabel_format(style='sci',scilimits=(-2,2))
    ax_amp.set_yticks([])
    
    ax5.set_ylabel(r' phase [$\pi$]')
    ax5.set_ylim(-1.0,1.0)
    ax5.legend(p_list,lab_list)


fig1.savefig(save_name + '.png', dpi=600)
#fig2.savefig(save_name + '_2.png', dpi=600)
plt.show()
print( '\n copy and paste: open -a preview ' + save_name + '.png')
#plt.show()
#plt.close()