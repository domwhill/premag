
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
bz_in = 0.0# selected magnetic field [Tesla]

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
path = path_list[-1]
t_list, tc_list = cf.get_t_list(path,var='Te')
time = '07'
#'07'
iy = 20
xmin,xmax = -10.0,60.0
thresh_N_ratio = 5e-2
var = sys.argv[1]

if len(sys.argv)>2:
    save_path = save_path  + sys.argv[-1] + '/'
save_name = '%sdeltaqy_wt_%s_%s_%s' % (save_path,save_tag,re.sub(' ','', var),time)


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


slice_array_y = lambda array: array[iy,:]


fig1 = fprl.newfig_generic_2yscale(1.0)#(1.1,scale_width=1.5,scale_ratio=0.5)#plt.figure()
gs1 = GS.GridSpec(1,2)
ax = plt.subplot(gs1[0])
ax2 = plt.subplot(gs1[1])
fig1.subplots_adjust(left=0.15, right=0.9, wspace = 0.55, top=0.88, bottom=0.2)


fig2 = fprl.newfig_generic_2yscale(1.0)#plt.figure()
ax3 = fig2.add_subplot(121)
ax4 = fig2.add_subplot(122)
ax5 = ax2.twinx()
fig2.subplots_adjust(left=0.15, right=0.9, wspace = 0.55, top=0.88, bottom=0.2)

x_list = np.array([110,120])

dict_c_lt1,dict_k_lt1 = q_mod.repack_2D(path_list[0],time)
nx,ny = np.shape(dict_k_lt1['Te'])

data_pos = dict_k_lt1['Te'][x_list,ny//2]

#------------------->>>
for ipp in range(len(path_list)):
    pp = path_list[ipp]

    dict_c,dict_k = q_mod.repack_2D(pp,time,recompute_q_c = True)
    np.save(pp.split('/')[-1] + 'dict_c.npy',dict_c,allow_pickle=True)
    np.save(pp.split('/')[-1] + 'dict_k.npy',dict_k,allow_pickle=True)
    x_grid = dict_k['x_grid']
    y_grid = dict_k['y_grid']
    data_k = dict_k[var]
    ipp = path_list.index(pp)
    

    # get ratios from particular point
    nx,ny = np.shape(dict_k['Te'])
    pT50, = ax2.plot(
                    x_grid*cd5.xstep_factor,dict_k['Te'][:,ny//2]*Te_factor,
                    c='k',
                    linestyle='-'
                    )
    
    # wt
    axwt = ax2.twinx()
    pT50, = axwt.plot(
                    x_grid*cd5.xstep_factor,dict_k['wt'][:,ny//2],c='k',
                    linestyle=style_list[ipp]
                    )
    # te_amp
    ax_amp = ax2.twinx()
    te_amp  = cf.get_U_dev_abs(dict_k['Te'])
    pT_amp, = ax_amp.semilogy(
                        x_grid*cd5.xstep_factor,np.max(np.abs(te_amp*Te_factor), axis=1),
                        c='k',
                        linestyle='-'
                    )
    

    axwt.set_ylabel(r'$\omega \tau$')
    ax2.set_ylabel(r'$T_e$ [eV]')
    #----->
    c_list = cmap(np.linspace(0.0, 1.0, len(x_list)))
    ix_list = cf.data_intersect(x_grid,dict_k['Te'][:,ny//2], data_pos)

    
    for indx,ix in enumerate(ix_list):
        colour = c_list[indx]
        
        pc, = ax.plot(
                        y_grid*cd5.xstep_factor,get_amp1(y_grid,dict_c[var][ix,:]),
                        c=colour,
                        linestyle='--'
                        )
    
        pk, = ax.plot(
                        y_grid*cd5.xstep_factor,get_amp1(y_grid,dict_k[var][ix,:]),
                        c=colour,
                        linestyle='-'
                        )
        
        # Plot dT lineout 
        # convert from eV -> keV

        pk, = ax3.plot(
                      y_grid*cd5.xstep_factor, get_amp(dict_k['Te'][ix,:])*Te_factor,
                      c=colour,
                      linestyle='--'
                      )
        # Plot wte y lineout
        pk, = ax4.plot(
                        y_grid*cd5.xstep_factor,get_amp(dict_k['wt'][ix,:]),
                        c=colour,
                        linestyle='-'
                        )
        
        # Mark points y amp is sampled at on Te                
        ax2.scatter(
                    x_grid[ix]*cd5.xstep_factor,dict_k['Te'][ix,ny//2]*Te_factor,
                    c=colour,marker='x'
                    )
                       
                    
    ax3.set_ylabel(r'$\delta T_e$ [eV]')
    ax4.set_ylabel(r'$\delta \omega \tau_{ei}$')

    ax.set_xlim(y_grid[2]*cd5.xstep_factor,y_grid[-2]*cd5.xstep_factor)
    ax2.set_xlim(xmin,xmax)
    ax.ticklabel_format(style='sci',scilimits=(-2,2))
    
    ax2.ticklabel_format(style='sci',scilimits=(-2,2))
    ax2.set_yticks([])
    ax3.ticklabel_format(style='sci',scilimits=(-2,2))
    ax4.ticklabel_format(style='sci',scilimits=(-2,2))

    ax.set_ylabel(ylab)
    ax.set_xlabel(xlab)
    ax.set_xlabel(xlab)



#fig1.savefig(save_name + '_1.png', dpi=600)
#fig2.savefig(save_name + '_2.png', dpi=600)
plt.show()
print( 'saving as: ', save_name + '.png')
print( '\n copy and paste: open -a preview ' + save_name + '_1.png')
print( '\n copy and paste: open -a preview ' + save_name + '_2.png')

#plt.show()
#plt.close()