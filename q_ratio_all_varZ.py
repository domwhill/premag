'''
    q spitzer - modified for variable Z
    
'''
#-----------------------------------------------------------------------
import matplotlib.pyplot as plt
import re,os,sys,site,getpass, numpy as np
userid = getpass.getuser()
path_pre = '../'
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')
site.addsitedir(path_pre)

##site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/')

#from scipy import interpolate as SI
import chfoil_module as cf
import q_SH_Te_tsteps as q_mod
#import kinetic_biermann_compare_lineouts as bier
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter
from pylab import *

import figure_prl_twocol as fprl
from chfoil_module import conv_factors_cd5
from chfoil_module import cd5_switches
import figure as fprl
import getpass,site
userid = getpass.getuser()



nv = 100
#Z = 3.5
ones_nv = np.ones((nv))
vmax = 10.0
#xmin,xmax,nx = 0.0,142.662,40
#ymin,ymax,ny = 0.0,35.0,10
vmin,vmax,nv = 0.0,10.0,100

vb_grid = np.arange(vmin,vmax,(vmax-vmin)/(nv+1))
v_grid = 0.5*(vb_grid[1:] + vb_grid[:-1])
dv = v_grid[1] - v_grid[0]
v2dv = (v_grid**2)*dv

#v_grid = np.arange(vmin,vmax,(vmax-vmin)/nv)
#---- get vel grid---
#ic = np.arange(nv)
#ib = np.arange(-0.5,nv+0.5)
#vb_grid = cf.interp_data(ic,v_grid,ib)
dvc = vb_grid[1:] - vb_grid[:-1]
Y1_weight = 4.0*np.pi/3.0
Y0_weight = 4.0*np.pi
#---------------------


path0 = 'chfoil_default5_heat11_2D1--cx1'
path0_nob = 'chfoil_default5_heat11_2D0_noB--cx1'
path0_stat = 'chfoil_default5_heat11_2D0_static--cx1'
path1 = 'chfoil_default5_heat11_2D1--cx1'
path1_nob = 'chfoil_default5_heat11_2D1_noB--cx1'
path3 = 'chfoil_default5_heat11_2D3--cx1'
path3_noB = 'chfoil_default5_heat11_2D3_noB--cx1'
path3_stat = 'chfoil_default5_heat11_2D3_static--cx1'
#path_list = [path3,path3_noB,path3_stat,path1,path1_nob]
path_1pp = 'chfoil_default5_heat11_2D1_1pp--cx1'
path_1pp_95lmfp = 'chfoil_default5_heat11_2D1_1pp_95lmfp--cx1'
#path_list = [path1]
path = 'cd5_eos7_2D1_regrid96_B0_5_static'
##path_pre = '../'
path = 'eos7_B1_static_1pp'
hallparam = '0'
#path = path_pre + 'eos7_B'+ hallparam + '_static_1pp_v2'
path = path_pre + 'eos7_B'+ hallparam + '_static_1pp'

#path = path_1pp
path = 'p400nFL_5v22_hv3_2D1_1pp_v2'
path = 'p4_4_v29_hv1_0_3vo2_2D1'
path = 'p400_noFL_5_v37_hv4_5vo2_t2_c'
path = 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E'
time = '10'
colormap='RdBu_r'
colormap='plasma'
x_lim = 220
save_path = '/Users/' +userid + '/Dropbox/York/Pre-magnetised/RESULTS/' + path.split('/')[-1]  + '/'
if not os.path.exists(save_path):
    os.system('mkdir ' + save_path)
    
path_tag = cf.retrieve_path_tag(path)
save_path = './pics/r_5v38/'
save_path = '/Users/' +userid + '/Dropbox/York/Pre-magnetised/RESULTS/' + path.split('/')[-1]  + '/'
if not os.path.exists(save_path):
    os.system('mkdir ' + save_path)
#path_tag = '1pp_B_' + hallparam#cf.retrieve_path_tag(path)
#path_tag = 'r_5_v38_Z_FEOS_MODNE_0_04vo2_2D1_noconvol'

#-----------------------------------------------------------------------
log_on = True
'''
SI_on = True
#--- inputs
v_te = 1.93e7 # ms^-1
tau_ei = 7.44e-14 # s
nu_ei = 1.0/tau_ei # s^-1
lambda_mfp = 1.437e-6 # m
#----------
if SI_on:
    xstep_factor = lambda_mfp*1e6
    tstep_factor = tau_ei*1e12
    xlab = r'$x$ [ $\mu m$ ]'
    ylab = r'$y$ [ $\mu m$ ]'
    leg_title = r'time [ $ps$ ]'

else:
    xstep_factor = 1.0
    tstep_factor = 1.0
    xlab = r'\textit{x [ $\lambda_{mfp}$ ]}'
    ylab = r'\textit{x [ $\lambda_{mfp}$ ]}'
    leg_title = r'\text{$time$ [t$_{col}$ ]}'
'''
SI_on = cd5_switches.SI_on
save_on = cd5_switches.save_on
hlines_on = cd5_switches.hlines_on
grid_on = cd5_switches.grid_on
horizontal_on = cd5_switches.horizontal_on
separate_plots_on = cd5_switches.separate_plots_on
'''
cwdpath = os.getcwd() +'/'+ path_pre
print cwdpath
#cd5 = cf.conv_factors_custom(cwdpath)
cd5 = cf.conv_factors_premag()
Z0 = cd5.Z
T0 = cd5.T0
n0 = cd5.n0
'''
norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' +norm_name + 'norm.log' 
[T0,n0,Z0,Bz0] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_path,Z0,Ar=6.51)

print 'Z0 = ', Z0
print 'T0 = ', T0
print 'n0 = ', n0
cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index,c_index = 0,500
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab
xlab_rel = cd5.xlab_rel
ylab = cd5.ylab
t_unit = cd5.t_unit
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
lh_style = cd5.lh_style
dashes = cd5.dashes

ylab = cd5.ylab
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
yi= cd5.yi

divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
#ptag = cf.retrieve_path_tag(path)
ptag = path_tag
fcmap = fprl.plotting_params.lineouts_cmap
#--- saving ----
ratios_savename = 'qSHqRLvN_ratios_' + ptag + 'DDFIW'


#-----------------------------------------------------------------------
def extend_grid_xy_to_vxy(nv,ny,nx,grid_xy):
    '''
        
    '''
    grid_vyx = np.zeros((nv,ny,nx))
    for iv in range(nv):
        grid_vyx[iv,:,:] = grid_xy
    return grid_vyx

def extend_grid_v_to_vxy(nv,ny,nx,grid_v):
    '''
        
    '''
    grid_vyx = np.zeros((nv,ny,nx))
    for iv in range(nv):
        grid_vyx[iv,:,:] = grid_v[iv]
    return grid_vyx


def extend_grid_y_to_vxy(nv,ny,nx,grid_y):
    '''
        
    '''
    grid_vyx = np.zeros((nv,ny,nx))
    for iy in range(ny):
        grid_vyx[:,iy,:] = grid_y[iy]
    return grid_vyx

def extend_grid_x_to_vxy(nv,ny,nx,grid_x):
    '''
        
    '''
    grid_vyx = np.zeros((nv,ny,nx))
    for ix in range(nx):
        grid_vyx[:,:,ix] = grid_x[ix]
    return grid_vyx
#-----------------------------------------------------------------------
def get_omega(v_grid,ni,Bz):
    nu_ei = (Z**2)*ni*(v_grid**-3)
    omega = Bz*(nu_ei**-1)
    return omega
    
def get_v2dv(v_grid):
    dv = v_grid[1] - v_grid[0]
    v2dv = (v_grid**2)*dv
    return v2dv

def get_v2dv_vyx(v_grid):
    dv = v_grid[1] - v_grid[0]
    v2dv = (v_grid**2)*dv
    v2dv_vyx = extend_grid_v_to_vxy(nv,ny,nx,grid_v)
    return v2dv_vyx
    
def maxw_dist(v,vte):
    f_om = ((np.pi*(vte**2))**(-1.5)) * np.exp(-(v/vte)**2) ;
    return f_om
    
def get_v_mom_m(v_grid,omega,F0,m):
    '''
        int_0^\inf dv V^{m+2}F0/(1 + omega^2*V^6)
    '''
    
    #ones_nv
    v2dv = get_v2dv(v_grid)
    int = ((v_grid**m)*v2dv*F0)/(ones_nv + (omega**2)*v_grid**6)
    mom = 4.0*np.pi*np.sum(int)
    return mom
    
def get_v_mom_m_n(v_grid,omega,F0,m,n):
    '''
        V^m_n
    '''
    mom_m = get_v_mom_m(v_grid,omega,F0,m)
    mom_n = get_v_mom_m(v_grid,omega,F0,n)
    mom_m_n = mom_m/mom_n
    return mom_m_n
    

def get_delta(v_grid,omega,F0):
    
    mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    delta = 1.0 + (omega**2)*(mom_8_5**2)
    ##print 'DELTA = ', delta
    return delta


#-----------------------------------------------------------------------    
def get_fo_mom_int(m,rho_vyx,omega_vyx,fo):
    prefactor = (4.0*np.pi/3.0)*rho_vyx # rho = 1/Z2ni
    # omega  = Bz/Z2ni = Bz*rho
    
    
    nv,ny,nx = np.shape(fo)
    #print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    ones_vyx = np.ones((nv,ny,nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv,ny,nx,v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv,ny,nx,v_grid)
    
    int = ((v_grid_vyx**m)*v2dv_vyx*fo)/(ones_vyx + (omega_vyx**2)*v_grid_vyx**6)
    mom = np.sum(prefactor*int,axis=0)
    return mom
#-----------------------------------------------------------------------    

def get_fo_mom_int2(m,rho_vyx,omega_vyx,fo):
    prefactor = (4.0*np.pi/3.0)*rho_vyx # rho = 1/Z2ni
    
    nv,ny,nx = np.shape(fo)
    ##print ' int 2 = '
    ##print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    ones_vyx = np.ones((nv,ny,nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv,ny,nx,v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv,ny,nx,v_grid)
    
    int = ((v_grid_vyx**m)*v2dv_vyx*fo)/((ones_vyx + (omega_vyx**2)*(v_grid_vyx**6))**2)
    mom = np.sum(prefactor*int,axis=0)
    return mom
#-----------------------------------------------------------------------

def get_dfodv_int(n,rho_vyx,omega_vyx,fo):
    # omega  = Bz/Z2ni = Bz*rho
    #print ' ---- dfodv --- int -----'
    v_mom_nm3 =  get_fo_mom_int(n-3,rho_vyx,omega_vyx,fo)
    v_mom2_np3 =  get_fo_mom_int2(n+3,rho_vyx,omega_vyx,fo)
    mom_out = -1.0*(n*v_mom_nm3 - 6.0*(omega_vyx[0,:,:]**2)*v_mom2_np3) # unclear whether this minus one should actually be there or not
    return mom_out
    
def get_I1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''
    
    #print ' ---- I1 -----'
    m = 5
    mom_x = get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_x)
    mom_y = get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_y)
    return mom_x,mom_y

def get_I2_scalar(rho_vyx,omega_vyx,fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' I2_scalar'
    m = 6
    mom = get_dfodv_int(m,rho_vyx,omega_vyx,fo)
    return mom

def get_I4_scalar(rho_vyx,omega_vyx,fo):
    '''
        I4 = get_I4_scalar(rho,omega,fo)
    '''
    #print ' I4_scalar'
    
    m = 9
    mom = rho_vyx[0,:,:]*get_dfodv_int(m,rho_vyx,omega_vyx,fo)
    return mom
    

def get_I3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y):
    '''
        I3 =  get_I3_vec(rho,omega,gradfo)
    '''
    #print ' ---- I3 -----'
    m = 8
    mom_x = rho_vyx[0,:,:]*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_x)
    mom_y = rho_vyx[0,:,:]*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_y)
    return mom_x,mom_y

def get_K1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''
    
    #print ' ---- K1 -----'
    m = 7
    mom_x = 0.5*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_x)
    mom_y = 0.5*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_y)
    return mom_x,mom_y

def get_K2_scalar(rho_vyx,omega_vyx,fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
    m = 8
    mom = 0.5*get_dfodv_int(m,rho_vyx,omega_vyx,fo)
    return mom

def get_K3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''
    
    #print ' ---- K3 -----'
    m = 10
    mom_x = 0.5*rho_vyx[0,:,:]*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_x)
    mom_y = 0.5*rho_vyx[0,:,:]*get_fo_mom_int(m,rho_vyx,omega_vyx,gradfo_y)
    return mom_x,mom_y


def get_K4_scalar(rho_vyx,omega_vyx,fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
    m = 11
    mom = 0.5*rho_vyx[0,:,:]*get_dfodv_int(m,rho_vyx,omega_vyx,fo)
    return mom


def get_v_N(grid,rho_vyx,Bz_vyx,fo):
    
    omega_vyx = Bz_vyx*rho_vyx
    gradfo_x,gradfo_y = cf.get_grad_3d(grid,fo)
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x,I1_y =  get_I1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    I3_x,I3_y =  get_I3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx,omega_vyx,fo)
    I4 = get_I4_scalar(rho_vyx,omega_vyx,fo)
    
    omega_yx = omega_vyx[0,:,:]
    Bz_yx = Bz_vyx[0,:,:]
    Eta = (I2 + (omega_yx**2)*((I4**2)/I2))
    v_nx = (1.0/Eta)*(I3_x - (I4/I2)*I1_x)
    v_ny = (1.0/Eta)*(I3_y - (I4/I2)*I1_y)
    
    
    return v_nx,v_ny

def get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo):
    '''
        alpha_perp_kin = get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo)
    '''
    
    omega_vyx = Bz_vyx*rho_vyx
    gradfo_x,gradfo_y = cf.get_grad_3d(grid,fo)
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x,I1_y =  get_I1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    I3_x,I3_y =  get_I3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx,omega_vyx,fo)
    I4 = get_I4_scalar(rho_vyx,omega_vyx,fo)
    
    omega_yx = omega_vyx[0,:,:]
    Bz_yx = Bz_vyx[0,:,:]
    Eta = (I2 + (omega_yx**2)*((I4**2)/I2))
    alpha_perp_kin  = -(1.0/Eta)
    
    return alpha_perp_kin

def data_var(data,var):
    dict = {}
    dict['data'] = data
    dict['name'] = var
    return dict
def get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo):
    '''
        dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    '''
    
    omega_vyx = Bz_vyx*rho_vyx
    gradfo_x,gradfo_y = cf.get_grad_3d(grid,fo)
    #print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x,I1_y =  get_I1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    I3_x,I3_y =  get_I3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    #print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx,omega_vyx,fo)
    I4 = get_I4_scalar(rho_vyx,omega_vyx,fo)
    #========
    K2 = get_K2_scalar(rho_vyx,omega_vyx,fo)
    K4 = get_K4_scalar(rho_vyx,omega_vyx,fo)    
    K1_x,K1_y =  get_K1_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    K3_x,K3_y =  get_K3_vec(rho_vyx,omega_vyx,gradfo_x,gradfo_y)
    
    
    
    omega_yx = omega_vyx[0,:,:]
    Bz_yx = Bz_vyx[0,:,:]
    Eta = (I2 + (omega_yx**2)*((I4**2)/I2))
    v_nx = -(1.0/Eta)*(I3_x - (I4/I2)*I1_x)
    v_ny = -(1.0/Eta)*(I3_y - (I4/I2)*I1_y)
    
    #--- Spitzer
    q_SH_x = K1_x - (K2/Eta)*I1_x - (K2/Eta)*(I4/I2)*(Bz_yx**2)*I3_x + (K4/Eta)*(Bz_yx**2)*(I3_x - I1_x*(I4/I2))
    q_SH_y = K1_y - (K2/Eta)*I1_y - (K2/Eta)*(I4/I2)*(Bz_yx**2)*I3_y + (K4/Eta)*(Bz_yx**2)*(I3_y - I1_y*(I4/I2))
    
    #---- RL
    kappa_w_dxT = - (K2/Eta)*(I3_x - I1_x*(I4/I2)) + K3_x - (K4/Eta)*I1_x - (Bz_yx**2)*(K4/Eta)*(I4/I2)*I3_x
    kappa_w_dyT = - (K2/Eta)*(I3_y - I1_y*(I4/I2)) + K3_y - (K4/Eta)*I1_y - (Bz_yx**2)*(K4/Eta)*(I4/I2)*I3_y
    q_RL_x = -Bz_yx*kappa_w_dyT
    q_RL_y = Bz_yx*kappa_w_dxT
    
    #--- Thermo electric perp
    beta_perp = K2/Eta + (Bz_yx**2)*(K4*I4)/(Eta*I2)
    q_TE_x = beta_perp*jx
    q_TE_y = beta_perp*jy
    
    #--- Thermo electric wedge
    beta_wedge = (-(K2*I4/(Eta*I2)) + K4/Eta)
    q_E_x = -Bz_yx*beta_wedge*jy
    q_E_y = Bz_yx*beta_wedge*jx
    
    
    dict ={}
    dict['q_SH_x'] = data_var(-1.0*q_SH_x,'q_x SH')
    dict['q_SH_y'] = data_var(-1.0*q_SH_y,'q_y SH')
    dict['q_RL_x'] = data_var(-1.0*q_RL_x,'q_x RL')
    dict['q_RL_y'] = data_var(-1.0*q_RL_y,'q_y RL')
    dict['q_TE_x'] = data_var(-1.0*q_TE_x,'q_x thermoelectric')
    dict['q_TE_y'] = data_var(-1.0*q_TE_y,'q_y thermoelectric')
    dict['q_E_x'] =  data_var(-1.0*q_E_x,'q_x Ettinghausen')
    dict['q_E_y'] =  data_var(-1.0*q_E_y,'q_y Ettinghausen')
    
    return dict

#-----------------------------------------------------------------------    

def get_alpha_perp(w,rho,ne,Te,F0,v_grid):
    '''
    w = eBz/m_e (gyrofreq)
    rZZni = rho
    vte**2 = 2.0*Te
    alpha_perp = get_alpha_perp(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho
    vte = (2.0*Te)**0.5

    omega  = w*rZZni 

    v_mom_5 = get_v_mom_m(v_grid,omega,F0,5)
    delta = get_delta(v_grid,omega,F0)
    #print 'shapes delta vmom5 ; ',np.shape(delta), np.shape(v_mom_5)
    alpha  = (1.5)*((vte*vte)/(ne*rZZni))*((v_mom_5*delta)**-1)
       
    return alpha

#-----------------------------------------------------------------------        
def get_alpha_wedge(w,rho,ne,Te,F0,v_grid):
    '''
        alpha_wedge = K *omega*( (1.5 V_8_5/(V_5*Delta)) - 1.0)
        alpha_wedge = get_alpha_wedge(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho  
    vte = (2.0*Te)**0.5

    omega  = w*rZZni
    v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_5 = get_v_mom_m(v_grid,omega,F0,5)
    delta = get_delta(v_grid,omega,F0)
    
    alpha  =  (omega/(rZZni*ne))*(((1.5)*(vte**2)*(v_mom_8_5/(v_mom_5*delta))) - 1.0)
    return alpha
#-----------------------------------------------------------------------
def get_beta_perp(w,rho,ne,Te,F0,v_grid):
    '''
        beta_perp = get_beta_perp(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho  
    vte = (2.0*Te)**0.5
    omega  = w*rZZni
    delta = get_delta(v_grid,omega,F0)
    v_mom_7_5 = get_v_mom_m_n(v_grid,omega,F0,7,5)
    v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_10_8 = get_v_mom_m_n(v_grid,omega,F0,10,8)
    
    beta = ( v_mom_7_5 + ((omega*v_mom_8_5)**2)*v_mom_10_8)/delta - 2.5
    return beta


#-----------------------------------------------------------------------
def get_beta_wedge(w,rho,ne,Te,F0,v_grid):
    '''
        beta_wedge = get_beta_wedge(w,rho,ne,Te,F0,v_grid)
        w = B-field Normalised
        rho = 1/(Z^2ni)
    '''
    rZZni = rho  
    vte = (2.0*Te)**0.5
    omega  = w*rZZni
 
    v_mom_7_5 = get_v_mom_m_n(v_grid,omega,F0,7,5)
    v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_10_8 = get_v_mom_m_n(v_grid,omega,F0,10,8)
    delta = get_delta(v_grid,omega,F0)
    beta = (1.0/(vte*vte))*omega*v_mom_8_5*(v_mom_10_8 - v_mom_7_5)/delta
    return beta


#-----------------------------------------------------------------------
def get_kappa_perp(w,rho,ne,Te,F0,v_grid):
    rZZni = rho  
    vte = (2.0*Te)**0.5

    omega  = w*rZZni
    v_mom_7_5 = get_v_mom_m_n(v_grid,omega,F0,7,5)
    v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_10_5 = get_v_mom_m_n(v_grid,omega,F0,10,5)
    v_mom_10_8 = get_v_mom_m_n(v_grid,omega,F0,10,8)

    v_mom_5 = get_v_mom_m(v_grid,omega,F0,5)
    v_mom_7 = get_v_mom_m(v_grid,omega,F0,7)
    v_mom_9 = get_v_mom_m(v_grid,omega,F0,9)


    delta = get_delta(v_grid,omega,F0)
    
    prefactor = (1.0/3.0)*((ne*rZZni)/(vte**4))
    
    kappa = prefactor* (v_mom_9-((v_mom_7_5 +(omega**2)*v_mom_8_5*v_mom_10_5)*v_mom_7/delta) +
            ((omega**2)*v_mom_8_5*v_mom_10_5*(v_mom_10_8 - v_mom_7_5)*v_mom_5/delta))
    return kappa

#-----------------------------------------------------------------------
def get_kappa_wedge(w,rho,ne,Te,F0,v_grid):
    rZZni = rho  
    vte = (2.0*Te)**0.5

    omega  = w*rZZni
    v_mom_7_5 = get_v_mom_m_n(v_grid,omega,F0,7,5)
    #v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_10_5 = get_v_mom_m_n(v_grid,omega,F0,10,5)
    #v_mom_10_8 = get_v_mom_m_n(v_grid,omega,F0,10,8)

    #v_mom_5 = get_v_mom_m(v_grid,omega,F0,5)
    v_mom_7 = get_v_mom_m(v_grid,omega,F0,7)
    v_mom_8 = get_v_mom_m(v_grid,omega,F0,8)
    #v_mom_9 = get_v_mom_m(v_grid,omega,F0,9)
    v_mom_12 = get_v_mom_m(v_grid,omega,F0,12)


    delta = get_delta(v_grid,omega,F0)
      
  
    kappa = (1.0/3.0)*((ne*rZZni)/(vte**4))*omega*(v_mom_12
             -(2.0*v_mom_10_5*v_mom_7/delta)
     +(((v_mom_7_5**2)-(omega*v_mom_10_5)**2)*v_mom_8/delta)) 

    return kappa
    
#-----------------------------------------------------------------------

def get_v_N_classical(Z2ni,ne,Te,w,dxT,dyT,jx=0.0,jy=0.0):
    
    #Te = 0.5
    vte = (2.0*Te)**0.5
    rho = Z2ni**-1
    
    ##print ' ######################################'
    ##print 'inputs : vte: %3.4f \t vmax: %3.4f \t \n \t w = %3.4f \n\t ne =  %3.4f \t rho %3.4f ' % (vte,vmax, w, ne ,rho)
    ##print ' ######################################'
    ##print ' v_grid ===== ', v_grid
    # ------
    
    F0 = maxw_dist(v_grid,vte)
    
    beta_wedge_c = get_beta_wedge(w,rho,ne,Te,F0,v_grid)
    beta_wedge = beta_wedge_c
    
    v_N_x = - beta_wedge*dxT/w
    v_N_y = - beta_wedge*dyT/w
    
    
    return v_N_x,v_N_y

#-----------------------------------------------------------------------

def get_vN_from_path(path,fprefix,time):
    '''
        v_nx,v_ny = get_vN_from_path(path,fprefix,time)
    '''

    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']

    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    #dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    
    
    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    grid = dict_fo
    nv,ny,nx = np.shape(fo)
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    rA = (Z2ni)**-1
    
    #-------

    dxT,dyT = get_gradT(x_grid,y_grid,Te)
    ##print 'qx_c = ', np.shape(qx_c), 'qy_c = ', np.shape(qy_c), np.shape(Bz_data), np.shape(T_data),np.shape(n_data)
    
    ny,nx = np.shape(qx_c)#len(T_data[:,0])-1,len(T_data[0,:])-1
    v_nx = np.zeros((ny,nx))
    v_ny = np.zeros((ny,nx))
    v_nx_hf = np.zeros((ny,nx))
    v_ny_hf = np.zeros((ny,nx))
        
    for ix in range(1,len(Te[0,:])-1):
       for iy in range(1,len(Te[:,0])-1):
            Z2ni_loc = Z2ni[iy,ix]
            v_nx[iy-1,ix-1],v_ny[iy-1,ix-1] = get_v_N_classical(Z2ni_loc,ne[iy,ix],Te[iy,ix],Bz[iy-1,ix-1],dxT[iy-1,ix-1],dyT[iy-1,ix-1])
            v_nx_hf[iy-1,ix-1] = qx_c[iy-1,ix-1]/(2.5*n[iy,ix]*T[iy,ix])
            v_ny_hf[iy-1,ix-1] = qy_c[iy-1,ix-1]/(2.5*n[iy,ix]*T[iy,ix])
    ##print 'shapes = ',np.shape(Bz),np.shape(v_nx)
    v_nx = np.where(Bz[1:-1,1:-1]==0.0,0.0,v_nx)
    v_ny = np.where(Bz[1:-1,1:-1]==0.0,0.0,v_ny)
    return v_nx,v_ny, v_nx_hf,v_ny_hf

#-----------------------------------------------------------------------
def get_q_SH(ne,Te,w,dxT,dyT):
    
    #Te = 0.5
    vte = (2.0*Te)**0.5
    w = 0.1
    #ne = 1.0
    rho = ne
    
    #print ' ######################################'
    #print 'inputs : vte: %3.4f \t vmax: %3.4f \t \n \t w = %3.4f \n\t ne =  %3.4f \t rho %3.4f ' % (vte,vmax, w, ne ,rho)
    #print ' ######################################'
    
    # ------
    F0 = maxw_dist(v_grid,vte)
    kappa_perp_c = get_kappa_perp(w,rho,ne,Te,F0,v_grid)
    ##print 'kp'
    kappa_wedge_c = get_kappa_wedge(w,rho,ne,Te,F0,v_grid)
    
    kappa_perp = kappa_perp_c*(Te**2.5)
    kappa_wedge = kappa_wedge_c*(Te**2.5)
    q_SH_x = - kappa_perp*dxT
    q_SH_y = - kappa_perp*dyT
    q_RL_x = + kappa_wedge*dyT
    q_RL_y = - kappa_wedge*dxT
    return q_SH_x

def get_gradT(x_grid,y_grid,T_data):
    '''
        ONLY FOR CC cells - centred differencing
    '''
    
    if len(np.shape(T_data))==1:
        dx = x_grid[2:] - x_grid[:-2]
        dxT = (T_data[2:] - T_data[:-2])/dx
        return dxT
    else:
        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        dxT = (T_data[:,2:] - T_data[:,:-2])/dx
        dyT = (T_data[2:,:] - T_data[:-2,:])/dy
        return dxT, dyT

def vmom_fo(m,v_grid,fo):
    '''
     mom = vmom_fo(m,v_grid,fo)
     fo can be 1D array or 3D [v,y,x]
    '''
    #dvc = get_dvc(v_grid)
    if len(np.shape(fo))>1:
        vmom = np.einsum('ijk,i->ijk',fo,Y0_weight*(v_grid**(m+2))*dvc)
    else:
        vmom = Y0_weight*fo*(v_grid**(m+2))*dvc
    mom = np.sum(vmom,axis=0)
    return mom
    
    
def convert_var_to_str(var):
    for name in globals():
        if eval(name)==var:
            return name

def plot_2D(ax,data,label,middle=None,colormap='jet',limits=None,xlabel=None,ylabel=None,clabel=None):
    if middle != None:
        try:
            norm = cf.MidPointNorm(middle)
        except ValueError:
            norm = None
    else:
        norm = None
    
    if limits:
        try:
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm,extent=limits)
            plt.colorbar(im,ax=ax,aspect='auto',norm=norm,label=clabel)
        except ValueError:
            norm=None
            #print 'setting norm to None;'
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm,extent=limits)
            plt.colorbar(im,ax=ax,aspect='auto',norm=norm,label=clabel)
            
    else:
        try:
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm)
            plt.colorbar(im,ax=ax,aspect='auto',norm=norm,label=clabel)
        except ValueError:
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm)
            plt.colorbar(im,ax=ax,aspect='auto',norm=norm,label=clabel)

    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
        
    #print ' plotted ---- ', label
    
#-----------------------------------------------------------------------
def vmom_f1(m,v_grid,f1x_c,f1y_c):
    '''
    
    '''
    #dvc = get_dvc(v_grid)
    
    if len(np.shape(f1x_c))>1:
        momx = np.einsum('ijk,i->ijk',f1x_c,Y1_weight*(v_grid**(m+3))*dvc)
        momy = np.einsum('ijk,i->ijk',f1y_c,Y1_weight*(v_grid**(m+3))*dvc)
    else:
        momx = Y1_weight*f1x_c*(v_grid**(m+3))*dvc     
        momy = Y1_weight*f1y_c*(v_grid**(m+3))*dvc        
    vmomx = np.sum(momx,axis=0)
    vmomy = np.sum(momy,axis=0)
    return vmomx,vmomy

def get_kinetic_heatflow_b(path,time):
    '''
        dict = get_kinetic_heatflow_b(path,time)
        var_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
        for var in range(len(dict.keys())):
            var_name = var_list[var]
            ax = fig.add_subplot(ny_plots,nx_plots,var+1)
            plot_2D(ax,dict[var_name]['data'],dict[var_name]['name'])
    
    '''
    
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']

    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])



    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    
    
    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    grid = dict_fo
    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    dict['x_grid']= dict_fo['x_grid']
    dict['y_grid']= dict_fo['y_grid']
    dict['v_grid']= dict_fo['v_grid']
    dict['time'] = dict_fo['time']
    return dict




def get_kinetic_heatflow(path,time):
    '''
       dict_x,dict_y = get_kinetic_heatflow(path,time)
    '''
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne

    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    
    
    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    
    
    dict_fxX = cf.load_dict(path,fpre(path),'fxX',time)
    fxX = dict_fxX['mat']

    dict_fxY = cf.load_dict(path,fpre(path),'fxY',time)
    fxY = dict_fxY['mat']

    dict_fyX = cf.load_dict(path,fpre(path),'fyX',time)
    fyX = dict_fyX['mat']

    dict_fyY = cf.load_dict(path,fpre(path),'fyY',time)
    fyY = dict_fyY['mat']

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    
    f1x_c = 0.5*(fxX[:,:,1:] + fxX[:,:,:-1])
    f1y_c = 0.5*(fyY[:,1:,:] + fyY[:,:-1,:])
    
    qx_c = 0.5*(qxX[:,1:] + qxX[:,:-1])
    qy_c = 0.5*(qyY[1:,:] + qyY[:-1,:])

    vmomx_0,vmomy_0 = vmom_f1(0,v_grid,f1x_c,f1y_c)
    vmomx_3,vmomy_3 = vmom_f1(3,v_grid,f1x_c,f1y_c)
    vmomx_5,vmomy_5 = vmom_f1(5,v_grid,f1x_c,f1y_c)
    
    jx = -vmomx_0
    jy = -vmomy_0
    
    mom0_3 = vmom_fo(3,v_grid,fo)
    mom0_5 = vmom_fo(5,v_grid,fo)
    mom0_7 = vmom_fo(7,v_grid,fo)

    ##print 'np.shape(fxX) = ', np.shape(fxX), np.shape(fyY), np.shape(mom0_3), np.shape(mom0_5)
    ##print ' shapes = ',np.shape(vmomx_3),np.shape(vmomy_3),np.shape(vmomx_5),np.shape(vmomy_5)
    ny,nx = np.shape(vmomy_5)
    ##print ' shape(Z2ni) = ', np.shape(Z2ni)
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = prof_Z*ne
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    ##print 'np.sort(Z2ni) = ', np.sort(Z2ni)[0,0], np.sort(fo)
    ##print ' np.shape(Z2ni) = ', np.shape(Z2ni)
    rA = (Z2ni*2.0)**-1

    vqx = rA*(vmomx_5 - (4.0/3.0)*(mom0_5/mom0_3)*(vmomx_3))
    vqy = rA*(vmomy_5 - (4.0/3.0)*(mom0_5/mom0_3)*(vmomy_3))

    vqBz_x = -1.0*vqy*Bz
    vqBz_y = -1.0*-1.0*vqx*Bz
    #---------------
    #print np.shape(mom0_7), np.shape(mom0_7)[:], np.shape(mom0_7)[:][0] 
    ny,nx = np.shape(mom0_7)
    
    x_grid = cf.trim_array(x_grid,nx,[1])
    y_grid = cf.trim_array(y_grid,ny,[1])
    
    ##print ' np.shape(x_Grid) = ', np.shape(x_grid),np.shape(y_grid), np.shape(mom0_7)
    dymom7,dxmom7 = cf.get_grad(y_grid,x_grid,mom0_7)
    dymom5,dxmom5 = cf.get_grad(y_grid,x_grid,mom0_5)
    dymom7,dxmom7 = -1.0*dymom7,-1.0*dxmom7
    dymom5,dxmom5 = -1.0*dymom5,-1.0*dxmom5
    
    ny,nx = np.shape(dxmom7)
    eta_star = Z2ni
    eta = Z2ni/(2.0*mom0_3)#Z2ni/(2vmom_0)
    #--- trim arrays
    ##print ' qx_c = ', np.shape(qx_c), np.shape(qy_c), np.shape(vqBz_y)
    array_list = [rA,eta,jx,jy,mom0_7,mom0_5,mom0_3,vqBz_x,vqBz_y]
    ##print ' nx = ',nx, 'ny = ', ny, np.shape(dxmom7)
    qx_c =  cf.trim_array(np.transpose(qx_c),ny,nx)
    qy_c = cf.trim_array(np.transpose(qy_c),ny,nx)
    ##print ' qx_c = ', np.shape(qx_c), np.shape(qy_c)
    #cf.trim_array(array,ny,nx)
    rA,eta,jx,jy,mom0_7,mom0_5,mom0_3,vqBz_x,vqBz_y = map(lambda array: cf.trim_array(array,ny,nx),array_list)
    
    resist_x = rA*(8.0/3.0)*eta*jx*mom0_5
    resist_y = rA*(8.0/3.0)*eta*jy*mom0_5
    
    
    grad1_x = rA*(1.0/3.0)*dxmom7
    grad1_y = rA*(1.0/3.0)*dymom7
    
    grad2_x = -(4.0/9.0)*rA*(dxmom5)*(mom0_5/mom0_3)
    grad2_y = -(4.0/9.0)*rA*(dymom5)*(mom0_5/mom0_3)
    #---------------
    
    qx_tot = grad1_x + grad2_x + resist_x + vqBz_x
    qy_tot = grad1_y + grad2_y + resist_y + vqBz_y
    dict_x = {}
    dict_x['x_grid'] = x_grid
    dict_x['grad1_x'] = grad1_x
    dict_x['grad2_x'] = grad2_x
    dict_x['resist_x'] = resist_x
    dict_x['vqBz_x'] = vqBz_x
    
    dict_y = {} 
    dict_y['y_grid'] = y_grid
    dict_y['grad1_y'] = grad1_y
    dict_y['grad2_y'] = grad2_y
    dict_y['resist_y'] = resist_y
    dict_y['vqBz_y'] = vqBz_y
    
    return dict_x,dict_y
def get_avgx(mat,ax=1):
    '''
        ax  = [n,..,2,1,0]
        ax = 1 x axis
        ax = 0 y axis
    '''
    
    avg = np.average(mat,axis=ax)
    return avg
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
def get_Nernst_ratio(path,time,sep_var=False):
    '''
    data_x,data_y = get_Nernst_ratio(path,time)
    '''
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor

    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #------ NERNST -----------------------------------------------------
    v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    v_nx_classical,v_ny_classical,v_nx_hf,v_ny_hf = q_mod.get_vN_from_path(path,fpre(path),time)
    if np.any(np.isnan(v_nx_classical)):
        print ' shape v_nx_classical = ', np.shape(v_nx_classical)
        print ' np.where(np.isnan(v_nx_classical))', np.where(np.isnan(v_nx_classical))
        print ' np.where(np.isnan(v_ny_classical))', np.where(np.isnan(v_ny_classical))
        
        sys.exit()
    v_nx_kinetic = np.transpose(v_nx)[cl_index:c_index,:]
    v_ny_kinetic = np.transpose(v_ny)[cl_index:c_index,:]
    
    #data_y = v_ny_kinetic/v_ny_classical[cl_index:c_index,:]
    #data_y = np.where(v_ny_classical[cl_index:c_index,:]==0.0,1.0,data_y)
    #vmin,vmax = -1.5,1.5
    #data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)
    
    
    
    if sep_var:
        print '--- Nernst separate ---'
        data_x = get_combined_dict(v_nx_classical[cl_index:c_index,:],v_nx_kinetic)
        data_y = get_combined_dict(v_ny_classical[cl_index:c_index,:],v_ny_kinetic)
    else:
        print ' --- nernst ratios ------'
        data_x = get_ratio_lim(v_nx_classical[cl_index:c_index,:],v_nx_kinetic,vmax=1.5)
        data_y = get_ratio_lim(v_ny_classical[cl_index:c_index,:],v_ny_kinetic,vmax=1.5)        
    #data_x = v_nx_kinetic/v_nx_classical[cl_index:c_index,:]
    #data_x = np.where(v_nx_classical[cl_index:c_index,:]==0.0,1.0,data_x)
    #vmin,vmax = -1.0,1.0
    #data_x = np.where(np.abs(data_x)>vmax,np.sign(data_x)*vmax,data_x)
    return data_x,data_y

#-----------------------------------------------------------------------
def get_Nernst_abs(path,time,sep_var=False):
    '''
    v_nx_kinetic,v_ny_kinetic,v_nx_classical[cl_index:c_index,:], v_ny_classical[cl_index:c_index,:] = get_Nernst_abs(path,time)
    '''
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor

    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #------ NERNST -----------------------------------------------------
    v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    v_nx_classical,v_ny_classical,v_nx_hf,v_ny_hf = q_mod.get_vN_from_path(path,fpre(path),time)
    if np.any(np.isnan(v_nx_classical)):
        print ' shape v_nx_classical = ', np.shape(v_nx_classical)
        print ' np.where(np.isnan(v_nx_classical))', np.where(np.isnan(v_nx_classical))
        print ' np.where(np.isnan(v_ny_classical))', np.where(np.isnan(v_ny_classical))
        
        sys.exit()
    v_nx_kinetic = np.transpose(v_nx)[cl_index:c_index,:]
    v_ny_kinetic = np.transpose(v_ny)[cl_index:c_index,:]
    
    #data_y = v_ny_kinetic/v_ny_classical[cl_index:c_index,:]
    #data_y = np.where(v_ny_classical[cl_index:c_index,:]==0.0,1.0,data_y)
    #vmin,vmax = -1.5,1.5
    #data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)
    
    
    
    if sep_var:
        print '--- Nernst separate ---'
        data_x = get_combined_dict(v_nx_classical[cl_index:c_index,:],v_nx_kinetic)
        data_y = get_combined_dict(v_ny_classical[cl_index:c_index,:],v_ny_kinetic)
    else:
        print ' --- nernst ratios ------'
        data_x = get_ratio_lim(v_nx_classical[cl_index:c_index,:],v_nx_kinetic,vmax=1.5)
        data_y = get_ratio_lim(v_ny_classical[cl_index:c_index,:],v_ny_kinetic,vmax=1.5)        
    #data_x = v_nx_kinetic/v_nx_classical[cl_index:c_index,:]
    #data_x = np.where(v_nx_classical[cl_index:c_index,:]==0.0,1.0,data_x)
    #vmin,vmax = -1.0,1.0
    #data_x = np.where(np.abs(data_x)>vmax,np.sign(data_x)*vmax,data_x)
    return v_nx_kinetic,v_ny_kinetic,v_nx_classical[cl_index:c_index,:], v_ny_classical[cl_index:c_index,:]


#-----------------------------------------------------------------------
def get_ratio_lim(q_yc,q_yk,vmax=1.5):
    print ' np.max(q_yc) min = ', np.min(np.abs(q_yc)), np.max(np.abs(q_yc))
    ratio_factor = np.abs(q_yk/q_yc)
    ratio_factor = np.where(np.isnan(q_yc),0.0,ratio_factor)
    ratio_factor = np.where(np.abs(q_yc)<=(1e-10)*np.max(np.abs(q_yc)),0.0,ratio_factor)
    #ratio_factor = np.where(np.abs(ratio_factor)>=10.0,np.sign(ratio_factor)*10.0,ratio_factor)
    data_y = ratio_factor
    #data_y = (q_yc - q_yk)/np.abs(q_yc)

    #data_y = np.where(q_yc==0.0,1.0,data_y)
    #data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)
    return data_y

#-----------------------------------------------------------------------
def get_combined_dict(q_yc,q_yk):
    comb_dict = {}
    comb_dict['classical'] = q_yc
    comb_dict['kinetic'] = q_yk
    return comb_dict
#-----------------------------------------------------------------------

def get_q_ratio(path,time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor

    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    #Z2ni = ne
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #heat flow ratios
    
    #--- kinetic
    dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    name_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[cl_index:c_index]
    
    q_KSH_x,q_KSH_y,q_KRL_x,q_KRL_y,q_KTE_x,q_KTE_y,q_KE_x,q_KE_y =dict['q_SH_x']['data'],dict['q_SH_y']['data'],dict['q_RL_x']['data'],dict['q_RL_y']['data'],dict['q_TE_x']['data'],dict['q_TE_y']['data'],dict['q_E_x']['data'],dict['q_E_y']['data']
    #--- classical

    
    q_dir = path + '/q_dir/' 
    q_SH_x = np.load(q_dir + 'q_SH_x.txt.npy')
    q_SH_y = np.load(q_dir + 'q_SH_y.txt.npy')
    q_RL_x = np.load(q_dir + 'q_RL_x.txt.npy')
    q_RL_y = np.load(q_dir + 'q_RL_y.txt.npy')
    q_E_x = np.load(q_dir + 'q_E_x.txt.npy')
    q_E_y = np.load(q_dir + 'q_E_y.txt.npy')
    q_xX = np.load(q_dir + 'q_xX.txt.npy')
    q_yY = np.load(q_dir + 'q_yY.txt.npy')
    #nt,ny,nx
    q_x_VFP = (q_xX[:,1:,:] + q_xX[:,:-1,:])*0.5
    q_y_VFP = (q_yY[:,:,1:] + q_yY[:,:,:-1])*0.5
    q_x_B = q_RL_x + q_E_x + q_SH_x
    q_y_B = q_RL_y + q_E_y + q_SH_y
    
    if np.isnan(np.sum(q_E_x)):
        print ' nans in Q_E_X'
        print ' SUM q_E_x = ', np.sum(q_E_x)
        print ' SUM q_RL_x = ', np.sum(q_RL_x)
        print ' SUM q_SH_x = ', np.sum(q_SH_x)
        print np.where(np.isnan(q_E_x))
    
    if np.isnan(np.sum(q_E_y)):
        print ' nans in Q_E_y'
        print ' SUM q_E_y = ', np.sum(q_E_y)
        print ' SUM q_RL_y = ', np.sum(q_RL_y)
        print ' SUM q_SH_y = ', np.sum(q_SH_y)
        print ' np.where(q_E_x)', np.where(np.isnan(q_E_y))
        

    
    #------ RL
    
    q_yk = q_KRL_y
    q_yc = q_RL_y[time_int,cl_index:c_index,:]   
    rat_RL_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int,cl_index:c_index,:]   
    rat_RL_x = get_ratio_lim(q_xc,q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int,cl_index:c_index,:]    
    rat_SH_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int,cl_index:c_index,:]    
    rat_SH_x = get_ratio_lim(q_xc,q_xk)
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int,cl_index:c_index,:]
    q_xc = q_x_B[time_int,cl_index:c_index,:]    
    rat_tot_x = get_ratio_lim(q_xc,q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int,cl_index:c_index,:]
    q_yc = q_y_B[time_int,cl_index:c_index,:]    
    rat_tot_y = get_ratio_lim(q_yc,q_yk)
    
    dict_ratio = {}
    dict_ratio['SH x'] = rat_SH_x#q_SH_x[time_int,cl_index:c_index,:]
    dict_ratio['SH y'] = rat_SH_y#q_SH_y[time_int,cl_index:c_index,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y
    
    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y
    return dict_ratio

#-----------------------------------------------------------------------

def get_q_abs(path,time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor

    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z= cf.load_dict(path,fpre(path),'Z',time)
    Z_prof = dict_Z['mat']
    
    
    Z2ni = ne*Z_prof
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #heat flow ratios
    
    #--- kinetic
    dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    name_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[cl_index:c_index]
    
    q_KSH_x,q_KSH_y,q_KRL_x,q_KRL_y,q_KTE_x,q_KTE_y,q_KE_x,q_KE_y =dict['q_SH_x']['data'],dict['q_SH_y']['data'],dict['q_RL_x']['data'],dict['q_RL_y']['data'],dict['q_TE_x']['data'],dict['q_TE_y']['data'],dict['q_E_x']['data'],dict['q_E_y']['data']
    #--- classical

    
    q_dir = path + '/q_dir/' 
    q_SH_x = np.load(q_dir + 'q_SH_x.txt.npy')
    q_SH_y = np.load(q_dir + 'q_SH_y.txt.npy')
    q_RL_x = np.load(q_dir + 'q_RL_x.txt.npy')
    q_RL_y = np.load(q_dir + 'q_RL_y.txt.npy')
    q_E_x = np.load(q_dir + 'q_E_x.txt.npy')
    q_E_y = np.load(q_dir + 'q_E_y.txt.npy')
    q_xX = np.load(q_dir + 'q_xX.txt.npy')
    q_yY = np.load(q_dir + 'q_yY.txt.npy')
    #nt,ny,nx
    q_x_VFP = (q_xX[:,1:,:] + q_xX[:,:-1,:])*0.5
    q_y_VFP = (q_yY[:,:,1:] + q_yY[:,:,:-1])*0.5
    q_x_B = q_RL_x + q_E_x + q_SH_x
    q_y_B = q_RL_y + q_E_y + q_SH_y
    
    if np.isnan(np.sum(q_E_x)):
        print ' nans in Q_E_X'
        print ' SUM q_E_x = ', np.sum(q_E_x)
        print ' SUM q_RL_x = ', np.sum(q_RL_x)
        print ' SUM q_SH_x = ', np.sum(q_SH_x)
        print np.where(np.isnan(q_E_x))
    
    if np.isnan(np.sum(q_E_y)):
        print ' nans in Q_E_y'
        print ' SUM q_E_y = ', np.sum(q_E_y)
        print ' SUM q_RL_y = ', np.sum(q_RL_y)
        print ' SUM q_SH_y = ', np.sum(q_SH_y)
        print ' np.where(q_E_x)', np.where(np.isnan(q_E_y))
        

    
    #------ RL
    
    q_yk = q_KRL_y
    q_yc = q_RL_y[time_int,cl_index:c_index,:]   
    rat_RL_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int,cl_index:c_index,:]   
    rat_RL_x = get_ratio_lim(q_xc,q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int,cl_index:c_index,:]    
    rat_SH_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int,cl_index:c_index,:]    
    rat_SH_x = get_ratio_lim(q_xc,q_xk)
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int,cl_index:c_index,:]
    q_xc = q_x_B[time_int,cl_index:c_index,:]    
    rat_tot_x = get_ratio_lim(q_xc,q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int,cl_index:c_index,:]
    q_yc = q_y_B[time_int,cl_index:c_index,:]    
    rat_tot_y = get_ratio_lim(q_yc,q_yk)
    
    
    v_nx_k,v_ny_k,v_nx_c, v_ny_c = get_Nernst_abs(path,time)
    dict_ratio_k = {}
    dict_ratio_k['SH x'] = q_KSH_x#q_SH_x[time_int,cl_index:c_index,:]
    dict_ratio_k['SH y'] = q_KSH_y#q_SH_y[time_int,cl_index:c_index,:]
    dict_ratio_k['RL x'] = q_KRL_x
    dict_ratio_k['RL y'] = q_KRL_y
    dict_ratio_k['vN x'] = v_nx_k
    dict_ratio_k['vN y'] = v_ny_k
    
    dict_ratio_c = {}
    dict_ratio_c['SH x'] = q_SH_x[time_int,cl_index:c_index,:]#q_SH_x[time_int,cl_index:c_index,:]
    dict_ratio_c['SH y'] = q_SH_y[time_int,cl_index:c_index,:]#q_SH_y[time_int,cl_index:c_index,:]
    dict_ratio_c['RL x'] = q_RL_x[time_int,cl_index:c_index,:]
    dict_ratio_c['RL y'] = q_RL_y[time_int,cl_index:c_index,:]
    dict_ratio_c['vN x'] = v_nx_c
    dict_ratio_c['vN y'] = v_ny_c

    
    return dict_ratio_k, dict_ratio_c


def get_q_individ(path,time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])

    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor

    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    nv,ny,nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    
    #Z2ni = np.transpose(cf.trim_array(Z2ni,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #heat flow ratios
    
    #--- kinetic
    dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    name_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[cl_index:c_index]
    
    q_KSH_x,q_KSH_y,q_KRL_x,q_KRL_y,q_KTE_x,q_KTE_y,q_KE_x,q_KE_y =dict['q_SH_x']['data'],dict['q_SH_y']['data'],dict['q_RL_x']['data'],dict['q_RL_y']['data'],dict['q_TE_x']['data'],dict['q_TE_y']['data'],dict['q_E_x']['data'],dict['q_E_y']['data']
    #--- classical

    
    q_dir = path + '/q_dir/' 
    q_SH_x = np.load(q_dir + 'q_SH_x.txt.npy')
    q_SH_y = np.load(q_dir + 'q_SH_y.txt.npy')
    q_RL_x = np.load(q_dir + 'q_RL_x.txt.npy')
    q_RL_y = np.load(q_dir + 'q_RL_y.txt.npy')
    q_E_x = np.load(q_dir + 'q_E_x.txt.npy')
    q_E_y = np.load(q_dir + 'q_E_y.txt.npy')
    q_xX = np.load(q_dir + 'q_xX.txt.npy')
    q_yY = np.load(q_dir + 'q_yY.txt.npy')
    #nt,ny,nx
    q_x_VFP = (q_xX[:,1:,:] + q_xX[:,:-1,:])*0.5
    q_y_VFP = (q_yY[:,:,1:] + q_yY[:,:,:-1])*0.5
    q_x_B = q_RL_x + q_E_x + q_SH_x
    q_y_B = q_RL_y + q_E_y + q_SH_y
    
    if np.isnan(np.sum(q_E_x)):
        print ' nans in Q_E_X'
        print ' SUM q_E_x = ', np.sum(q_E_x)
        print ' SUM q_RL_x = ', np.sum(q_RL_x)
        print ' SUM q_SH_x = ', np.sum(q_SH_x)
        print np.where(np.isnan(q_E_x))
    
    if np.isnan(np.sum(q_E_y)):
        print ' nans in Q_E_y'
        print ' SUM q_E_y = ', np.sum(q_E_y)
        print ' SUM q_RL_y = ', np.sum(q_RL_y)
        print ' SUM q_SH_y = ', np.sum(q_SH_y)
        print ' np.where(q_E_x)', np.where(np.isnan(q_E_y))
        

    
    #------ RL
    q_yk = q_KRL_y
    q_yc = q_RL_y[time_int,cl_index:c_index,:]   
    rat_RL_y = get_combined_dict(q_yc,q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int,cl_index:c_index,:]   
    rat_RL_x = get_combined_dict(q_xc,q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int,cl_index:c_index,:]    
    rat_SH_y = get_combined_dict(q_yc,q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int,cl_index:c_index,:]    
    rat_SH_x = get_combined_dict(q_xc,q_xk)
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int,cl_index:c_index,:]
    q_xc = q_x_B[time_int,cl_index:c_index,:]
    rat_tot_x = get_combined_dict(q_xc,q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int,cl_index:c_index,:]
    q_yc = q_y_B[time_int,cl_index:c_index,:]    
    rat_tot_y = get_combined_dict(q_yc,q_yk)
    
    
    
    dict_ratio = {}
    dict_ratio['SH x'] = rat_SH_x#q_SH_x[time_int,cl_index:c_index,:]
    dict_ratio['SH y'] = rat_SH_y#q_SH_y[time_int,cl_index:c_index,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y
    
    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y
    return dict_ratio
    
def get_alpha_perp_path(path,time):
    '''
     alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time)
    '''
    #--- load in data ....
    time_int = int(time)
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']
    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
    jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])
    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    #Z2ni = ne
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor
    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']


    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    nv,ny,nx = np.shape(fo)

    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    
    # trim data to correct size...

    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))

    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    alpha_perp_kinetic = get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo)
    
    
    #---- perform the classical calculation
    T_data = np.transpose(cf.trim_array(Te,nx,ny))
    n_data = np.transpose(cf.trim_array(ne,nx,ny))
    Bz_data = Bz
    #---------------------------------------
    alpha_perp_classical = np.zeros((np.shape(T_data)))
    #------
    #print ' shapes = ', np.shape(T_data),np.shape(n_data), np.shape(Bz_data), np.shape(alpha_perp_kinetic)
    #for ix in range(1,len(T_data[0,:])-1):
    #   for iy in range(1,len(T_data[:,0])-1):
    for ix in range(nx):
        for iy in range(ny):
            Z2ni_loc = n_data[iy,ix]
            Te_loc = T_data[iy,ix]
            ne_loc = n_data[iy,ix]
            w_loc = Bz_data[iy,ix]
            #v_nx[iy-1,ix-1],v_ny[iy-1,ix-1] = get_v_N_classical(Z2ni,n_data[iy,ix],T_data[iy,ix],Bz_data[iy-1,ix-1],dxT[iy-1,ix-1],dyT[iy-1,ix-1])
            #v_nx_hf[iy-1,ix-1] = q_xc[iy-1,ix-1]/(2.5*n_data[iy,ix]*T_data[iy,ix])
            #v_ny_hf[iy-1,ix-1] = q_yc[iy-1,ix-1]/(2.5*n_data[iy,ix]*T_data[iy,ix])

            vte = (2.0*Te_loc)**0.5
            rho = Z2ni_loc**-1
            F0 = maxw_dist(v_grid,vte)
            alpha_perp_classical[iy,ix] = get_alpha_perp(w_loc,rho,ne_loc,Te_loc,F0,v_grid)
    #----------
    #dict[name]['data'] = np.transpose(dict[name]['data'])[cl_index:c_index]
    return alpha_perp_classical,alpha_perp_kinetic

def fpre(pp):
    return pp.split('/')[-1]

if __name__ == "__main__":
    '''
        dxT = T_data[1:] - T_data[:-1]/(x_grid[1:] - x_grid[:-1])
    '''
    
    
    dict_Bz = cf.load_dict(path,fpre(path),'Bz',time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path,fpre(path),'jxX',time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path,fpre(path),'jyY',time)
    jyY = dict_jyY['mat']
    #print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    try:
        jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
        jx_c = 0.5*(jxX[1:,:] + jxX[:-1,:])
    except IndexError:
        jy_c = 0.5*(jyY[:,1:] + jyY[:,:-1])
        jx_c = 0.5*(jxX[1:] + jxX[:-1])
        
    dict_qxX = cf.load_dict(path,fpre(path),'qxX',time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path,fpre(path),'qyY',time)
    qyY = dict_qyY['mat']
    #print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    #qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
    #qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    try:
        qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
        qx_c = 0.5*(qxX[1:,:] + qxX[:-1,:])
    except IndexError:
        qy_c = 0.5*(qyY[:,1:] + qyY[:,:-1])
        qx_c = 0.5*(qxX[1:] + qxX[:-1])


    dict_wt = cf.load_dict(path,fpre(path),'wt',time)
    wt = dict_wt['mat']

    
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    
    dict_te = cf.load_dict(path,fpre(path),'Te',time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid*xstep_factor
    y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time']*tstep_factor
    
    
    dict_fo = cf.load_dict(path,fpre(path),'fo',time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    #print ' -------------------- SHAPES --------------------------------'
    
    grid = dict_fo
    nv,ny,nx = np.shape(fo)
    dict_ne = cf.load_dict(path,fpre(path),'n',time)
    ne = dict_ne['mat']
    dict_Z = cf.load_dict(path,fpre(path),'Z',time)
    prof_Z = dict_Z['mat']
    nv,ny,nx = np.shape(fo)

    Z2ni = np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    # #print ' np.shapae(jx_c) = ', np.shape(jx_c),np.shape(jy_c)
    # #print ' np.shape(Z2ni) = ', np.shape(Z2ni),' np.shape(Bz) = ', np.shape(Bz)
    
    rA = (Z2ni)**-1
    
    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid,fo)
    rho = rA
    omega = Bz*rA
    rho_vyx =  extend_grid_xy_to_vxy(nv,ny,nx,rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv,ny,nx,omega)
    #------ NERNST -----------------------------------------------------
    v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo) # kinetic Nernst
    v_nx_classical,v_ny_classical,v_nx_hf,v_ny_hf = q_mod.get_vN_from_path(path,fpre(path),time)

    dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    
    #print ' dict keys = ', dict.keys()
    if False:
        #fig = plt.figure()
        fig = fprl.newfig_generic_2yscale(width,scale_width=1.0,scale_ratio=1.0,clearon=True)

        ny_plots, nx_plots = 1,2
        var_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
        
        var_list = ['q_SH_x','q_SH_y']
        q_SH_x,q_SH_y,q_RL_x,q_RL_y,q_TE_x,q_TE_y,q_E_x,q_E_y =dict['q_SH_x']['data'],dict['q_SH_y']['data'],dict['q_RL_x']['data'],dict['q_RL_y']['data'],dict['q_TE_x']['data'],dict['q_TE_y']['data'],dict['q_E_x']['data'],dict['q_E_y']['data']
        qx_sum = q_SH_x + q_RL_x + q_TE_x + q_E_x
        qy_sum = q_SH_y + q_RL_y + q_TE_y + q_E_y

    print ' np.shape( y_grid_SI) = ', np.shape(y_grid_SI),np.shape(x_grid_SI),c_index,cl_index   
    limits = [y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[cl_index]]
    lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    lab_dict['lims'] = limits#[y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'],lab_dict['ylab'] = xlab,ylab
    lab_dict['cbar_title'] = '   '
    lab_dict['title'] = '   '
    lab_dict['middle'] = 0.0

    left_cbar = 0.05
    bottom_cbar = 0.05
    right_cbar = 0.8
    top_cbar = 0.9
    hspace_cbar = 0
    wspace_cbar = 0
    #---- 2D plots of Nernst and ting
    if False:     
        plt.subplots_adjust(left=left_cbar, bottom=bottom_cbar, right=right_cbar, top=top_cbar, wspace=wspace_cbar, hspace=hspace_cbar)
    
        #-------------------------------------------------------------------
        var = 'v_nx_k'
        save_name = save_path + '/'+ path_tag + '_' + var +  '_' + str(int(time_col)) + '.pdf'
        lab_dict['cbar_title'] = '$v_{N,x}$ [$v_{T}$]'
        lab_dict['save_name'] = save_name
        data = np.transpose(v_nx)[cl_index:c_index,:]
        min_vx,max_vx = np.min(data),np.max(data)
        lab_dict['vmin'] = min_vx
        lab_dict['vmax'] = max_vx

        var = 'v_nx_ratio'
        save_name = save_path + '/'+ path_tag + '_' + var +  '_' + str(int(time_col)) + '.pdf'
        lab_dict['cbar_title'] = '$v_{N,x,K}/v_{N,x,B}$ [$v_{T}$]'
        lab_dict['save_name'] = save_name
        v_nx_kinetic = np.transpose(v_nx)[cl_index:c_index,:]
        data = v_nx_kinetic/v_nx_classical[cl_index:c_index,:]
        data = np.where(v_nx_classical[cl_index:c_index,:]==0.0,0.0,data)
        lab_dict['vmin'] = -1.0
        lab_dict['vmax'] = 1.0
        data = np.where(np.abs(data)>=lab_dict['vmax'],np.sign(data)*lab_dict['vmax'],data)
    
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.contourf(data,aspect='auto',
                        cmap=colormap,norm=None,extent=lab_dict['lims'],
                        vmin=lab_dict['vmin'],vmax=lab_dict['vmax'])
        plt.colorbar(im,ax=ax,aspect='auto',norm=None,label=lab_dict['cbar_title'])
        print ' path = ', path
        var = 'v_ny_ratio'
        save_name = save_path + '/'+ path_tag + '_' + var +  '_' + str(int(time_col)) + '.pdf'
        lab_dict['cbar_title'] = '$v_{N,y,K}/v_{N,y,B}$ [$v_{T}$]'
        lab_dict['save_name'] = save_name
        v_ny_kinetic = np.transpose(v_ny)[cl_index:c_index,:]
        data = v_ny_kinetic/v_ny_classical[cl_index:c_index,:]
        data = np.where(v_ny_classical[cl_index:c_index,:]==0.0,0.0,data)
        lab_dict['vmin'] = -1.0
        lab_dict['vmax'] = 1.0
        data = np.where(np.abs(data)>=lab_dict['vmax'],np.sign(data)*lab_dict['vmax'],data)
        ##norm = lab_dict['norm']
    
    
    
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        im = ax.contourf(data,aspect='auto',
                        cmap=colormap,norm=None,extent=lab_dict['lims'],
                        vmin=lab_dict['vmin'],vmax=lab_dict['vmax'])
        plt.colorbar(im,ax=ax,aspect='auto',norm=None,label=lab_dict['cbar_title'])
        #cf.produce_2D_fig(x_grid_SI,data,lab_dict,lineout_list=[])
        #print ' path = ', path
        #plt.show()

    nt = 15

    
    #------- plot the ratios now
    if False:
        width = 1.0
        fig = fprl.newfig_generic_3yscale(width,scale_width=1.0,scale_ratio=0.5,clearon=True)
        fig.clf()

        [ax1,ax2,ax3] = [fig.add_subplot(1,3,1),fig.add_subplot(1,3,2),fig.add_subplot(1,3,3)]


        ax_list = [ax1,ax2,ax3]
        map(fprl.fix_ax,[ax1,ax2,ax3])
        #plt.subplots_adjust(hspace=0.1,left=0.2,right=0.95,wspace=0.68)
        fig.subplots_adjust(hspace=0.5,left=0.15,right=0.96,wspace=0.5,bottom=0.2,top=0.9)
        
        t_list = [3,8,12,15]
        t_list = [3,5]

        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        
        c_list = iter(fcmap(np.linspace(0,1,len(t_list))))
        log_on = True
        lab_list = []
        plot_list = []
        for tt in t_list:
            c= next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt
        
            dict_ratio = get_q_ratio(path,time)
            data_x2d_SH,data_y2d_SH,prelab_SH = dict_ratio['SH x'],dict_ratio['SH y'],'\perp'
            data_x2d_RL,data_y2d_RL,prelab_RL = dict_ratio['RL x'],dict_ratio['RL y'],'RL'
            data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time)
            
            
            
            prelab_vN = 'v_N'
            prelab_SH = 'q_{\perp,x}'
            prelab_RL = 'q_{RL,y}'
            prelab_vN = 'v_{N,x}'
            #prelab = 'v_{N'
            data_xSH = get_avgx(data_x2d_SH)
            data_ySH = get_avgx(data_y2d_SH)
            data_xRL = get_avgx(data_x2d_RL)
            data_yRL = get_avgx(data_y2d_RL)
            data_xvN = get_avgx(data_x2d_vN)
            data_yvN = get_avgx(data_y2d_vN)
            
            #p1, =ax1.plot(x_grid_SI[cl_index:c_index],data_x,c=c)

            p1, =ax1.semilogy(x_grid_SI[cl_index:c_index],data_xSH,c=c)
            p2, =ax2.semilogy(x_grid_SI[cl_index:c_index],data_xvN,c=c)  
            p3, =ax3.semilogy(x_grid_SI[cl_index:c_index],data_yRL,c=c)
            #--- legend stuff
            dict_wt = cf.load_dict(path,fpre(path),'wt',time)
            t_l = dict_wt['time']*tstep_factor
            t_lab = '%3i \si{ps}' % t_l
            t_lab = '%3i' % t_l
            
            plot_list.append(p1)
            lab_list.append(t_lab)
        ##prelab_SH,prelab_vN,prelab_RL = 
        ytit = lambda ax, name: ax.set_ylabel('$ %s,_K/ %s,_B$' % (name,name))
        #ytit = lambda ax, name: ax.set_ylabel('$\Theta_{%s,K}/ \Theta_{%s,B}$' % (name,name))
        def clearxaxis(ax):
            ax.set_xlabel('')
            ax.set_xticklabels([])        
            
        map(ytit,[ax1,ax2,ax3],[prelab_SH,prelab_vN,prelab_RL])
        #map(clearxaxis,[ax1,ax2,ax3])
        map(lambda axy: axy.set_xlabel(xlab_rel), [ax1,ax2,ax3])
        map(lambda axy: axy.xaxis.set_major_locator(MaxNLocator(4,prune='both')), [ax1,ax2,ax3])
        #map(lambda axy: axy.set_xlabel(xlab), [ax1,ax2,ax3,ax4,ax5,ax6])

        xlim = ax1.get_xlim()
        map(lambda axy: axy.set_xlim(10.0,xlim[-1]),ax_list)
                
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
        #map(lambda axy: sort_y(axy),[ax1,ax5,ax3])
        # - qSHy, vNy
        #map(lambda axy: axy.set_ylim(0.5e-2,1.5e2),[ax2,ax4,ax6])
        #map(sort_y_ylim,[ax2,ax4,ax6])
        map(lambda axy: axy.set_xlim(-10.0,100.0),[ax1,ax2,ax3])
        map(lambda axy: axy.set_ylim(0.5,1.2),[ax1,ax2,ax3])

        #path = '/Users/'+userid+'/Dropbox/PhD/Poster_Presentations/DDFIW-2018/'
        #ax1.set_ylabel(xtit)
        #ax2.set_ylabel(ytit)
        #ax1.set_xlabel(xlab)
        #ax2.set_xlabel(xlab)
        if False:
            ax5.legend(plot_list,lab_list,
                    bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,
                    title=r't $\si{ps}$')
        if log_on:
            save_name = save_path + 'q_ratio_' +path_tag + '.pdf'
            fig.savefig(save_name)
            print('saving fig as ... ', save_name)
        else:
            save_name = save_path + 'q_ratio_' +path_tag + '.pdf'
            fig.savefig(save_name)
            print('saving fig as ... ', save_name)
    #fprl.savefig_thesis(ratios_savename)
    #plt.show()

    #------- plot the ratios now
    if True:
        #fig = plt.figure()
        #plt.close(fig)
        #fig = fprl.newfig_generic_2yscale(1.0,scale_width=1.0,scale_ratio=1.0,clearon=True)
        #fig = fprl.newfig_generic(1.1)
        #fig2 = fprl.newfig_generic(1.1)
        width = 1.0
        fig = fprl.newfig_generic_3yscale(width,scale_width=1.0,scale_ratio=0.5,clearon=True)
        fig.clf()
        #plt.clf()
        #gs1 = GS.GridSpec(2,3,width_ratios=[1.0,1.0,1.0])
        #pl
        #[ax1,ax2,ax3] = [fig.add_subplot(2,3,1),fig.add_subplot(2,3,6),fig.add_subplot(2,3,2)]
        #[ax4,ax5,ax6]= [fig.add_subplot(2,3,4),fig.add_subplot(2,3,3),fig.add_subplot(2,3,5)]

        [ax1,ax2,ax3] = [fig.add_subplot(1,3,1),fig.add_subplot(1,3,2),fig.add_subplot(1,3,3)]


        ax_list = [ax1,ax2,ax3]
        map(fprl.fix_ax,[ax1,ax2,ax3])
        #plt.subplots_adjust(hspace=0.1,left=0.2,right=0.95,wspace=0.68)
        fig.subplots_adjust(hspace=0.5,left=0.15,right=0.96,wspace=0.5,bottom=0.2,top=0.9)
        
        t_list = [3,8,12,15]
        t_list = [10]

        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        
        c_list = iter(fcmap(np.linspace(0,1,len(t_list))))
        log_on = True
        lab_list = []
        plot_list = []
        for tt in t_list:
            c= next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt
            print ' path again = ', path
            dict_ratio_k,dict_ratio_c = get_q_abs(path,time)
            data_xK_SH,data_yK_SH,prelab_SH = dict_ratio_k['SH x'],dict_ratio_k['SH y'],'\perp'
            data_xK_RL,data_yK_RL,prelab_RL = dict_ratio_k['RL x'],dict_ratio_k['RL y'],'RL'
            data_xK_vN,data_yK_vN,prelab_RL = dict_ratio_k['vN x'],dict_ratio_k['vN y'],'v_N'
            
            data_xC_SH,data_yC_SH,prelab_SH = dict_ratio_c['SH x'],dict_ratio_c['SH y'],'\perp'
            data_xC_RL,data_yC_RL,prelab_RL = dict_ratio_c['RL x'],dict_ratio_c['RL y'],'RL'
            data_xC_vN,data_yC_vN,prelab_RL = dict_ratio_c['vN x'],dict_ratio_c['vN y'],'v_N'

            
            
            
            prelab_vN = 'v_N'
            prelab_SH = 'q_{\perp,x}'
            prelab_RL = 'q_{RL,y}'
            prelab_vN = 'v_{N,x}'
            #prelab = 'v_{N'
            if False:
                data_xSH = get_avgx(data_x2d_SH)
                data_ySH = get_avgx(data_y2d_SH)
                data_xRL = get_avgx(data_x2d_RL)
                data_yRL = get_avgx(data_y2d_RL)
                data_xvN = get_avgx(data_x2d_vN)
                data_yvN = get_avgx(data_y2d_vN)
            
            #p1, =ax1.plot(x_grid_SI[cl_index:c_index],data_x,c=c)
            #custom_dat = lambda data: (get_avgx(data))
            custom_dat = lambda data: ((data[:,10]))

            funcplot  = lambda ax,data,lstyle: ax.plot(x_grid_SI[cl_index:c_index],custom_dat(data),c=c,linestyle=lstyle)
            p1, = funcplot(ax1,data_xK_SH,'-')
            p2, = funcplot(ax2,data_xK_vN,'-')
            p3, = funcplot(ax3,data_yK_RL,'-')

            p1, = funcplot(ax1,data_xC_SH,'--')
            p2, = funcplot(ax2,data_xC_vN,'--')
            p3, = funcplot(ax3,data_yC_RL,'--')

            #--- legend stuff
            
            dict_wt = cf.load_dict(path,fpre(path),'wt',time)
            t_l = dict_wt['time']*tstep_factor
            t_lab = '%3i \si{ps}' % t_l
            print 't_lab = ',t_lab
            t_lab = '%3i' % t_l
            
            plot_list.append(p1)
            lab_list.append(t_lab)
        ##prelab_SH,prelab_vN,prelab_RL = 
        ytit = lambda ax, name: ax.set_ylabel('$ %s,_K/ %s,_B$' % (name,name))
        ytit = lambda ax, name: ax.set_ylabel('$%s$' % (name))
        
        #ytit = lambda ax, name: ax.set_ylabel('$\Theta_{%s,K}/ \Theta_{%s,B}$' % (name,name))
        def clearxaxis(ax):
            ax.set_xlabel('')
            ax.set_xticklabels([])        
            
        map(ytit,[ax1,ax2,ax3],[prelab_SH,prelab_vN,prelab_RL])
        #map(clearxaxis,[ax1,ax2,ax3])
        map(lambda axy: axy.set_xlabel(xlab_rel), [ax1,ax2,ax3])
        map(lambda axy: axy.xaxis.set_major_locator(MaxNLocator(4,prune='both')), [ax1,ax2,ax3])

        map(lambda axy: axy.yaxis.set_major_locator(MaxNLocator(6,prune='both')), [ax1,ax2,ax3])
        #map(lambda axy: axy.set_xlim(-6.0,12.0), [ax1,ax2,ax3])
        
        #map(lambda axy: axy.set_xlabel(xlab), [ax1,ax2,ax3,ax4,ax5,ax6])

        xlim = ax1.get_xlim()
        #map(lambda axy: axy.set_xlim(10.0,xlim[-1]),ax_list)
                
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
        #map(lambda axy: sort_y(axy),[ax1,ax5,ax3])
        # - qSHy, vNy
        #map(lambda axy: axy.set_ylim(0.5e-2,1.5e2),[ax2,ax4,ax6])
        #map(sort_y_ylim,[ax2,ax4,ax6])
        #map(lambda axy: axy.set_xlim(-10.0,100.0),[ax1,ax2,ax3])
        map(lambda axy: axy.set_xlim(-10.0,30.0),[ax1,ax2,ax3])

        #map(lambda axy: axy.set_ylim(0.5,1.2),[ax1,ax2,ax3])

        #path = '/Users/'+userid+'/Dropbox/PhD/Poster_Presentations/DDFIW-2018/'
        #ax1.set_ylabel(xtit)
        #ax2.set_ylabel(ytit)
        #ax1.set_xlabel(xlab)
        #ax2.set_xlabel(xlab)

        print 'path again = ', path
        if log_on:
            save_name = save_path + 'qx_vN_qrl_lineouts' + path + str(tt) + '.pdf'
            fig.savefig(save_name)
            print ' saving fig as: ', save_name
        else:
            save_name = save_path + 'qx_vN_qrl_lineouts' + path + str(tt) + '.pdf'
            fig.savefig(save_name)
            print ' saving fig as: ', save_name
        plt.clf()
        #plt.show()  
    sys.exit()
    
  
    plt.close(fig)
    #sys.exit()
