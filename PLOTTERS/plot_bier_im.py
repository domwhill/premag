'''

    This script plots image plots of dBdt due to biermann

'''
import numpy as np, sys, os, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/'+ userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/MODULES')
import matplotlib.pyplot as plt
#---> kinetic/classical transport post processing module
import kinetic_ohmslaw_module_varZ as kohb
import figure_prl_twocol as fprl
from pylab import *
import chfoil_module as cf
import house_keeping as hk
import pdb

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12


#---> file inputs
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir
norm_dir = paths.norm_dir
log_file = norm_dir + 'norm.log'
[Te_ref,n_ref,Z_ref,Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir,Z_ref,Ar=6.51)
scale_length = 1
lambda_p = 5

bz_in = 400
path = paths.get_path(scale_length, bz_in, lambda_p)

save_name = '%s/biermann_im_LT%i_bz%i_lp%i.png' % (save_path,scale_length,bz_in,lambda_p)

#----------------
cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index,c_index = 0,-1
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
v_te = cd5.v_te
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab_rel
ylab = cd5.ylab
Bz_ref  = (m_e/(q_e*tau_ei))

ylab_bier = r'$\partial_t \mathbf{B}|_{B}$'
ylab_bier = r'$\partial_t \mathbf{B}|_{B}$ [$\si{Ts^{-1}}$]'

xmin,xmax = 0.0,20.0
t_list = [6]
t_list, tc_list = cf.get_t_list(path,var='Te')
time = t_list[-1]


#-----------------------------------------------------------------------------------------------------
def custom_im(fig, ax,xgrid,ygrid,data,lab,**kwargs):
    vmin = kwargs.pop('vmin',None)
    vmax = kwargs.pop('vmax', None)

    xmin =  0.0
    xmax = 20.0
    ymin,ymax = np.min(ygrid),np.max(ygrid)
    bool_t = (xgrid>=xmin)*(xgrid<xmax)
    treat_data = lambda a: np.sign(data[bool_t,:])*np.log10(np.abs(a[bool_t,:]))
    treat_data = lambda a: a[bool_t]

    im = ax.imshow(
                treat_data(data),
                aspect='auto',extent=[ymin,ymax,xmax,xmin],
                vmin=vmin,vmax=vmax,
                cmap = 'RdBu_r'
                )
    #cax = ax
    fprl.custom_colorbar(fig, im, ax, cax='', lab=lab)
    #plt.colorbar(im,ax=ax,label = lab)
    return im

#-----------------------------------------------------------------------------------------------------
def fpre(path):
    return path.split('/')[-1]


#-----------------------------------------------------------------------------------------------------


log_on = True
lab_list = []
plot_list = []


dict_fo = cf.load_dict(path,path.split('/')[-1],'fo','00')
x_grid_SI = dict_fo['x_grid'][1:-1]*xstep_factor
y_grid_SI = dict_fo['y_grid']*xstep_factor

classic_biermann, kinetic_biermann = kohb.get_kinetic_and_classicB_c(path, fpre(path), time)

fig = fprl.newfig_generic_2yscale(1.0)
ax1,ax2 = fig.add_subplot(121), fig.add_subplot(122)
im = custom_im(
        fig,
        ax1,x_grid_SI, y_grid_SI, classic_biermann, ylab_bier
        )

custom_im(
        fig,
        ax2,x_grid_SI, y_grid_SI, kinetic_biermann, ylab_bier
        )

plt.show()
#------------> plotting

