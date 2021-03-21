'''
    This script is the same as 'plot_ohms_law_components_2D_varZ.py'
    except it calls 
'''
#-----------------------------------------------------------------------

import sys, os, re, getpass, site, numpy as np
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/')
import kinetic_ohmslaw_module_varZ as kbier

import matplotlib.pyplot as plt
import kinetic_ohmslaw_module_varZ as q_mod
import kinetic_ohmslaw_module_varZ as bier
import matplotlib as mpl
from pylab import *
from matplotlib import ticker
#
import chfoil_module as cf
#import figure_prl_twocol as fprl
from chfoil_module import conv_factors_cd5
from chfoil_module import cd5_switches
import figure as fprl

nv = 100
Z = 5.9
ones_nv = np.ones((nv))
vmax = 10.0
#xmin,xmax,nx = 0.0,142.662,40
#ymin,ymax,ny = 0.0,35.0,10
vmin, vmax, nv = 0.0, 10.0, 100

vb_grid = np.arange(vmin, vmax, (vmax - vmin) / (nv + 1))
v_grid = 0.5 * (vb_grid[1:] + vb_grid[:-1])
dv = v_grid[1] - v_grid[0]
v2dv = (v_grid**2) * dv

#v_grid = np.arange(vmin,vmax,(vmax-vmin)/nv)
#---- get vel grid---
#ic = np.arange(nv)
#ib = np.arange(-0.5,nv+0.5)
#vb_grid = cf.interp_data(ic,v_grid,ib)
dvc = vb_grid[1:] - vb_grid[:-1]
Y1_weight = 4.0 * np.pi / 3.0
Y0_weight = 4.0 * np.pi
#---------------------

path_pre = '../PREMAG/2D_RUNS/'

path_speck = path_pre + 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E'

path = path_speck

time = '10'
colormap = 'RdBu_r'
colormap = 'plasma'
x_lim = 500
save_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/RESULTS/' + path.split('/')[-1] + '/'
if not os.path.exists(save_path):
    os.system('mkdir ' + save_path)
path_tag = cf.retrieve_path_tag(path)
#-----------------------------------------------------------------------
log_on = True

SI_on = cd5_switches.SI_on
save_on = cd5_switches.save_on
hlines_on = cd5_switches.hlines_on
grid_on = cd5_switches.grid_on
horizontal_on = cd5_switches.horizontal_on
separate_plots_on = cd5_switches.separate_plots_on

norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + 'norm.log'
[T0, n0, Z0, Bz0] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_path, Z0, Ar=6.51)
cl_index = cd5.cl_index
c_index = cd5.c_index
SI_on = cd5.SI_on
v_te = cd5.v_te
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab
xlab = cd5.xlab_rel

ylab = cd5.ylab
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
yi = cd5.yi

divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
ptag = cf.retrieve_path_tag(path)
fcmap = fprl.plotting_params.lineouts_cmap
#--- saving ----
ratios_savename = 'qSHqRLvN_ratios_' + ptag


#-----------------------------------------------------------------------
def extend_grid_xy_to_vxy(nv, ny, nx, grid_xy):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iv in range(nv):
        grid_vyx[iv, :, :] = grid_xy
    return grid_vyx


def extend_grid_v_to_vxy(nv, ny, nx, grid_v):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iv in range(nv):
        grid_vyx[iv, :, :] = grid_v[iv]
    return grid_vyx


def extend_grid_y_to_vxy(nv, ny, nx, grid_y):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iy in range(ny):
        grid_vyx[:, iy, :] = grid_y[iy]
    return grid_vyx


def extend_grid_x_to_vxy(nv, ny, nx, grid_x):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for ix in range(nx):
        grid_vyx[:, :, ix] = grid_x[ix]
    return grid_vyx


#-----------------------------------------------------------------------
def get_omega(v_grid, ni, Bz):
    nu_ei = (Z**2) * ni * (v_grid**-3)
    omega = Bz * (nu_ei**-1)
    return omega


def get_v2dv(v_grid):
    dv = v_grid[1] - v_grid[0]
    v2dv = (v_grid**2) * dv
    return v2dv


def get_v2dv_vyx(v_grid):
    dv = v_grid[1] - v_grid[0]
    v2dv = (v_grid**2) * dv
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, grid_v)
    return v2dv_vyx


def maxw_dist(v, vte):
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    return f_om


def get_v_mom_m(v_grid, omega, F0, m):
    '''
        int_0^\inf dv V^{m+2}F0/(1 + omega^2*V^6)
    '''

    #ones_nv
    v2dv = get_v2dv(v_grid)
    int = ((v_grid**m) * v2dv * F0) / (ones_nv + (omega**2) * v_grid**6)
    mom = 4.0 * np.pi * np.sum(int)
    return mom


def get_v_mom_m_n(v_grid, omega, F0, m, n):
    '''
        V^m_n
    '''
    mom_m = get_v_mom_m(v_grid, omega, F0, m)
    mom_n = get_v_mom_m(v_grid, omega, F0, n)
    mom_m_n = mom_m / mom_n
    return mom_m_n


def get_delta(v_grid, omega, F0):

    mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    delta = 1.0 + (omega**2) * (mom_8_5**2)
    ##print 'DELTA = ', delta
    return delta


#-----------------------------------------------------------------------
def get_fo_mom_int(m, rho_vyx, omega_vyx, fo):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni
    # omega  = Bz/Z2ni = Bz*rho

    nv, ny, nx = np.shape(fo)
    #print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    ones_vyx = np.ones((nv, ny, nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / (ones_vyx + (omega_vyx**2) * v_grid_vyx**6)
    mom = np.sum(prefactor * int, axis=0)
    return mom


#-----------------------------------------------------------------------


def get_fo_mom_int2(m, rho_vyx, omega_vyx, fo):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni

    nv, ny, nx = np.shape(fo)
    ##print ' int 2 = '
    ##print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    ones_vyx = np.ones((nv, ny, nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / ((ones_vyx + (omega_vyx**2) * (v_grid_vyx**6))**2)
    mom = np.sum(prefactor * int, axis=0)
    return mom


#-----------------------------------------------------------------------


def get_dfodv_int(n, rho_vyx, omega_vyx, fo):
    # omega  = Bz/Z2ni = Bz*rho
    #print ' ---- dfodv --- int -----'
    v_mom_nm3 = get_fo_mom_int(n - 3, rho_vyx, omega_vyx, fo)
    v_mom2_np3 = get_fo_mom_int2(n + 3, rho_vyx, omega_vyx, fo)
    mom_out = -1.0 * (n * v_mom_nm3 - 6.0 * (omega_vyx[0, :, :]**2) * v_mom2_np3
                     )    # unclear whether this minus one should actually be there or not
    return mom_out


def get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- I1 -----'
    m = 5
    mom_x = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x)
    mom_y = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y)
    return mom_x, mom_y


def get_I2_scalar(rho_vyx, omega_vyx, fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' I2_scalar'
    m = 6
    mom = get_dfodv_int(m, rho_vyx, omega_vyx, fo)
    return mom


def get_I4_scalar(rho_vyx, omega_vyx, fo):
    '''
        I4 = get_I4_scalar(rho,omega,fo)
    '''
    #print ' I4_scalar'

    m = 9
    mom = rho_vyx[0, :, :] * get_dfodv_int(m, rho_vyx, omega_vyx, fo)
    return mom


def get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y):
    '''
        I3 =  get_I3_vec(rho,omega,gradfo)
    '''
    #print ' ---- I3 -----'
    m = 8
    mom_x = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x)
    mom_y = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y)
    return mom_x, mom_y


def get_K1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- K1 -----'
    m = 7
    mom_x = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x)
    mom_y = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y)
    return mom_x, mom_y


def get_K2_scalar(rho_vyx, omega_vyx, fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
    m = 8
    mom = 0.5 * get_dfodv_int(m, rho_vyx, omega_vyx, fo)
    return mom


def get_K3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- K3 -----'
    m = 10
    mom_x = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x)
    mom_y = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y)
    return mom_x, mom_y


def get_K4_scalar(rho_vyx, omega_vyx, fo):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
    m = 11
    mom = 0.5 * rho_vyx[0, :, :] * get_dfodv_int(m, rho_vyx, omega_vyx, fo)
    return mom


def get_v_N(grid, rho_vyx, Bz_vyx, fo):

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d(grid, fo)
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    v_nx = (1.0 / Eta) * (I3_x - (I4 / I2) * I1_x)
    v_ny = (1.0 / Eta) * (I3_y - (I4 / I2) * I1_y)

    return v_nx, v_ny


def get_alpha_perp_kinetic(grid, rho_vyx, Bz_vyx, fo):
    '''
        alpha_perp_kin = get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo)
    '''

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d(grid, fo)
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    alpha_perp_kin = -(1.0 / Eta)

    return alpha_perp_kin


def data_var(data, var):
    dict = {}
    dict['data'] = data
    dict['name'] = var
    return dict


def get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo):
    '''
        dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    '''

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d(grid, fo)
    #print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    #print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo)
    #========
    K2 = get_K2_scalar(rho_vyx, omega_vyx, fo)
    K4 = get_K4_scalar(rho_vyx, omega_vyx, fo)
    K1_x, K1_y = get_K1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)
    K3_x, K3_y = get_K3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    v_nx = -(1.0 / Eta) * (I3_x - (I4 / I2) * I1_x)
    v_ny = -(1.0 / Eta) * (I3_y - (I4 / I2) * I1_y)

    #--- Spitzer
    q_SH_x = K1_x - (K2 / Eta) * I1_x - (K2 / Eta) * (I4 / I2) * (Bz_yx**2) * I3_x + (K4 / Eta) * (
        Bz_yx**2) * (I3_x - I1_x * (I4 / I2))
    q_SH_y = K1_y - (K2 / Eta) * I1_y - (K2 / Eta) * (I4 / I2) * (Bz_yx**2) * I3_y + (K4 / Eta) * (
        Bz_yx**2) * (I3_y - I1_y * (I4 / I2))

    #---- RL
    kappa_w_dxT = -(K2 / Eta) * (I3_x - I1_x * (I4 / I2)) + K3_x - (K4 / Eta) * I1_x - (
        Bz_yx**2) * (K4 / Eta) * (I4 / I2) * I3_x
    kappa_w_dyT = -(K2 / Eta) * (I3_y - I1_y * (I4 / I2)) + K3_y - (K4 / Eta) * I1_y - (
        Bz_yx**2) * (K4 / Eta) * (I4 / I2) * I3_y
    q_RL_x = -Bz_yx * kappa_w_dyT
    q_RL_y = Bz_yx * kappa_w_dxT

    #--- Thermo electric perp
    beta_perp = K2 / Eta + (Bz_yx**2) * (K4 * I4) / (Eta * I2)
    q_TE_x = beta_perp * jx
    q_TE_y = beta_perp * jy

    #--- Thermo electric wedge
    beta_wedge = (-(K2 * I4 / (Eta * I2)) + K4 / Eta)
    q_E_x = -Bz_yx * beta_wedge * jy
    q_E_y = Bz_yx * beta_wedge * jx

    dict = {}
    dict['q_SH_x'] = data_var(-1.0 * q_SH_x, 'q_x SH')
    dict['q_SH_y'] = data_var(-1.0 * q_SH_y, 'q_y SH')
    dict['q_RL_x'] = data_var(-1.0 * q_RL_x, 'q_x RL')
    dict['q_RL_y'] = data_var(-1.0 * q_RL_y, 'q_y RL')
    dict['q_TE_x'] = data_var(-1.0 * q_TE_x, 'q_x thermoelectric')
    dict['q_TE_y'] = data_var(-1.0 * q_TE_y, 'q_y thermoelectric')
    dict['q_E_x'] = data_var(-1.0 * q_E_x, 'q_x Ettinghausen')
    dict['q_E_y'] = data_var(-1.0 * q_E_y, 'q_y Ettinghausen')

    return dict


#-----------------------------------------------------------------------


def get_alpha_perp(w, rho, ne, Te, F0, v_grid):
    '''
    w = eBz/m_e (gyrofreq)
    rZZni = rho
    vte**2 = 2.0*Te
    alpha_perp = get_alpha_perp(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho
    vte = (2.0 * Te)**0.5

    omega = w * rZZni

    v_mom_5 = get_v_mom_m(v_grid, omega, F0, 5)
    delta = get_delta(v_grid, omega, F0)
    #print 'shapes delta vmom5 ; ',np.shape(delta), np.shape(v_mom_5)
    alpha = (1.5) * ((vte * vte) / (ne * rZZni)) * ((v_mom_5 * delta)**-1)

    return alpha


#-----------------------------------------------------------------------
def get_alpha_wedge(w, rho, ne, Te, F0, v_grid):
    '''
        alpha_wedge = K *omega*( (1.5 V_8_5/(V_5*Delta)) - 1.0)
        alpha_wedge = get_alpha_wedge(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho
    vte = (2.0 * Te)**0.5

    omega = w * rZZni
    v_mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    v_mom_5 = get_v_mom_m(v_grid, omega, F0, 5)
    delta = get_delta(v_grid, omega, F0)

    alpha = (omega / (rZZni * ne)) * (((1.5) * (vte**2) * (v_mom_8_5 / (v_mom_5 * delta))) - 1.0)
    return alpha


#-----------------------------------------------------------------------
def get_beta_perp(w, rho, ne, Te, F0, v_grid):
    '''
        beta_perp = get_beta_perp(w,rho,ne,Te,F0,v_grid)
    '''
    rZZni = rho
    vte = (2.0 * Te)**0.5
    omega = w * rZZni
    delta = get_delta(v_grid, omega, F0)
    v_mom_7_5 = get_v_mom_m_n(v_grid, omega, F0, 7, 5)
    v_mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    v_mom_10_8 = get_v_mom_m_n(v_grid, omega, F0, 10, 8)

    beta = (v_mom_7_5 + ((omega * v_mom_8_5)**2) * v_mom_10_8) / delta - 2.5
    return beta


#-----------------------------------------------------------------------
def get_beta_wedge(w, rho, ne, Te, F0, v_grid):
    '''
        beta_wedge = get_beta_wedge(w,rho,ne,Te,F0,v_grid)
        w = B-field Normalised
        rho = 1/(Z^2ni)
    '''
    rZZni = rho
    vte = (2.0 * Te)**0.5
    omega = w * rZZni

    v_mom_7_5 = get_v_mom_m_n(v_grid, omega, F0, 7, 5)
    v_mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    v_mom_10_8 = get_v_mom_m_n(v_grid, omega, F0, 10, 8)
    delta = get_delta(v_grid, omega, F0)
    beta = (1.0 / (vte * vte)) * omega * v_mom_8_5 * (v_mom_10_8 - v_mom_7_5) / delta
    return beta


#-----------------------------------------------------------------------
def get_kappa_perp(w, rho, ne, Te, F0, v_grid):
    rZZni = rho
    vte = (2.0 * Te)**0.5

    omega = w * rZZni
    v_mom_7_5 = get_v_mom_m_n(v_grid, omega, F0, 7, 5)
    v_mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    v_mom_10_5 = get_v_mom_m_n(v_grid, omega, F0, 10, 5)
    v_mom_10_8 = get_v_mom_m_n(v_grid, omega, F0, 10, 8)

    v_mom_5 = get_v_mom_m(v_grid, omega, F0, 5)
    v_mom_7 = get_v_mom_m(v_grid, omega, F0, 7)
    v_mom_9 = get_v_mom_m(v_grid, omega, F0, 9)

    delta = get_delta(v_grid, omega, F0)

    prefactor = (1.0 / 3.0) * ((ne * rZZni) / (vte**4))

    kappa = prefactor * (v_mom_9 - ((v_mom_7_5 +
                                     (omega**2) * v_mom_8_5 * v_mom_10_5) * v_mom_7 / delta) +
                         ((omega**2) * v_mom_8_5 * v_mom_10_5 *
                          (v_mom_10_8 - v_mom_7_5) * v_mom_5 / delta))
    return kappa


#-----------------------------------------------------------------------
def get_kappa_wedge(w, rho, ne, Te, F0, v_grid):
    rZZni = rho
    vte = (2.0 * Te)**0.5

    omega = w * rZZni
    v_mom_7_5 = get_v_mom_m_n(v_grid, omega, F0, 7, 5)
    #v_mom_8_5 = get_v_mom_m_n(v_grid,omega,F0,8,5)
    v_mom_10_5 = get_v_mom_m_n(v_grid, omega, F0, 10, 5)
    #v_mom_10_8 = get_v_mom_m_n(v_grid,omega,F0,10,8)

    #v_mom_5 = get_v_mom_m(v_grid,omega,F0,5)
    v_mom_7 = get_v_mom_m(v_grid, omega, F0, 7)
    v_mom_8 = get_v_mom_m(v_grid, omega, F0, 8)
    #v_mom_9 = get_v_mom_m(v_grid,omega,F0,9)
    v_mom_12 = get_v_mom_m(v_grid, omega, F0, 12)

    delta = get_delta(v_grid, omega, F0)

    kappa = (1.0 / 3.0) * ((ne * rZZni) /
                           (vte**4)) * omega * (v_mom_12 - (2.0 * v_mom_10_5 * v_mom_7 / delta) +
                                                (((v_mom_7_5**2) -
                                                  (omega * v_mom_10_5)**2) * v_mom_8 / delta))

    return kappa


#-----------------------------------------------------------------------


def get_v_N_classical(Z2ni, ne, Te, w, dxT, dyT, jx=0.0, jy=0.0):

    #Te = 0.5
    vte = (2.0 * Te)**0.5
    rho = Z2ni**-1

    ##print ' ######################################'
    ##print 'inputs : vte: %3.4f \t vmax: %3.4f \t \n \t w = %3.4f \n\t ne =  %3.4f \t rho %3.4f ' % (vte,vmax, w, ne ,rho)
    ##print ' ######################################'
    ##print ' v_grid ===== ', v_grid
    # ------

    F0 = maxw_dist(v_grid, vte)

    beta_wedge_c = get_beta_wedge(w, rho, ne, Te, F0, v_grid)
    beta_wedge = beta_wedge_c

    v_N_x = -beta_wedge * dxT / w
    v_N_y = -beta_wedge * dyT / w

    return v_N_x, v_N_y


#-----------------------------------------------------------------------


def get_vN_from_path(path, fprefix, time):
    '''
        v_nx,v_ny = get_vN_from_path(path,fprefix,time)
    '''

    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']

    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    grid = dict_fo
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))
    rA = (Z2ni)**-1

    #-------

    dxT, dyT = get_gradT(x_grid, y_grid, Te)
    ##print 'qx_c = ', np.shape(qx_c), 'qy_c = ', np.shape(qy_c), np.shape(Bz_data), np.shape(T_data),np.shape(n_data)

    ny, nx = np.shape(qx_c)    #len(T_data[:,0])-1,len(T_data[0,:])-1
    v_nx = np.zeros((ny, nx))
    v_ny = np.zeros((ny, nx))
    v_nx_hf = np.zeros((ny, nx))
    v_ny_hf = np.zeros((ny, nx))

    for ix in range(1, len(Te[0, :]) - 1):
        for iy in range(1, len(Te[:, 0]) - 1):
            Z2ni = ne[iy, ix]
            v_nx[iy - 1,
                 ix - 1], v_ny[iy - 1,
                               ix - 1] = get_v_N_classical(Z2ni, ne[iy, ix], Te[iy, ix], Bz[iy - 1,
                                                                                            ix - 1],
                                                           dxT[iy - 1, ix - 1], dyT[iy - 1, ix - 1])
            v_nx_hf[iy - 1, ix - 1] = qx_c[iy - 1, ix - 1] / (2.5 * n[iy, ix] * T[iy, ix])
            v_ny_hf[iy - 1, ix - 1] = qy_c[iy - 1, ix - 1] / (2.5 * n[iy, ix] * T[iy, ix])
    ##print 'shapes = ',np.shape(Bz),np.shape(v_nx)
    v_nx = np.where(Bz[1:-1, 1:-1] == 0.0, 0.0, v_nx)
    v_ny = np.where(Bz[1:-1, 1:-1] == 0.0, 0.0, v_ny)
    return v_nx, v_ny, v_nx_hf, v_ny_hf


#-----------------------------------------------------------------------
def get_q_SH(ne, Te, w, dxT, dyT):

    #Te = 0.5
    vte = (2.0 * Te)**0.5
    w = 0.1
    #ne = 1.0
    rho = ne

    #print ' ######################################'
    #print 'inputs : vte: %3.4f \t vmax: %3.4f \t \n \t w = %3.4f \n\t ne =  %3.4f \t rho %3.4f ' % (vte,vmax, w, ne ,rho)
    #print ' ######################################'

    # ------
    F0 = maxw_dist(v_grid, vte)
    kappa_perp_c = get_kappa_perp(w, rho, ne, Te, F0, v_grid)
    ##print 'kp'
    kappa_wedge_c = get_kappa_wedge(w, rho, ne, Te, F0, v_grid)

    kappa_perp = kappa_perp_c * (Te**2.5)
    kappa_wedge = kappa_wedge_c * (Te**2.5)
    q_SH_x = -kappa_perp * dxT
    q_SH_y = -kappa_perp * dyT
    q_RL_x = +kappa_wedge * dyT
    q_RL_y = -kappa_wedge * dxT
    return q_SH_x


def get_gradT(x_grid, y_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
    '''

    if len(np.shape(T_data)) == 1:
        dx = x_grid[2:] - x_grid[:-2]
        dxT = (T_data[2:] - T_data[:-2]) / dx
        return dxT
    else:
        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        dxT = (T_data[:, 2:] - T_data[:, :-2]) / dx
        dyT = (T_data[2:, :] - T_data[:-2, :]) / dy
        return dxT, dyT


def vmom_fo(m, v_grid, fo):
    '''
     mom = vmom_fo(m,v_grid,fo)
     fo can be 1D array or 3D [v,y,x]
    '''
    #dvc = get_dvc(v_grid)
    if len(np.shape(fo)) > 1:
        vmom = np.einsum('ijk,i->ijk', fo, Y0_weight * (v_grid**(m + 2)) * dvc)
    else:
        vmom = Y0_weight * fo * (v_grid**(m + 2)) * dvc
    mom = np.sum(vmom, axis=0)
    return mom


def convert_var_to_str(var):
    for name in globals():
        if eval(name) == var:
            return name


def plot_2D(ax,
            data,
            label,
            middle=None,
            colormap='jet',
            limits=None,
            xlabel=None,
            ylabel=None,
            clabel=None):
    if middle != None:
        try:
            norm = cf.MidPointNorm(middle)
        except ValueError:
            norm = None
    else:
        norm = None

    if limits:
        try:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm, extent=limits)
            plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)
        except ValueError:
            norm = None
            #print 'setting norm to None;'
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm, extent=limits)
            plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)

    else:
        try:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm)
            plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)
        except ValueError:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm)
            plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)

    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    #print ' plotted ---- ', label


#-----------------------------------------------------------------------
def vmom_f1(m, v_grid, f1x_c, f1y_c):
    '''
    
    '''
    #dvc = get_dvc(v_grid)

    if len(np.shape(f1x_c)) > 1:
        momx = np.einsum('ijk,i->ijk', f1x_c, Y1_weight * (v_grid**(m + 3)) * dvc)
        momy = np.einsum('ijk,i->ijk', f1y_c, Y1_weight * (v_grid**(m + 3)) * dvc)
    else:
        momx = Y1_weight * f1x_c * (v_grid**(m + 3)) * dvc
        momy = Y1_weight * f1y_c * (v_grid**(m + 3)) * dvc
    vmomx = np.sum(momx, axis=0)
    vmomy = np.sum(momy, axis=0)
    return vmomx, vmomy


def get_kinetic_heatflow_b(path, time):
    '''
        dict = get_kinetic_heatflow_b(path,time)
        var_list = ['q_SH_x','q_SH_y','q_RL_x','q_RL_y','q_TE_x','q_TE_y','q_E_x','q_E_y']
        for var in range(len(dict.keys())):
            var_name = var_list[var]
            ax = fig.add_subplot(ny_plots,nx_plots,var+1)
            plot_2D(ax,dict[var_name]['data'],dict[var_name]['name'])
    
    '''
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']

    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    grid = dict_fo
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))
    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)
    dict['x_grid'] = dict_fo['x_grid']
    dict['y_grid'] = dict_fo['y_grid']
    dict['v_grid'] = dict_fo['v_grid']
    dict['time'] = dict_fo['time']
    return dict


def convert_name_to_label(name):
    comp = name[-1]
    s = re.search('q_(?P<type>\w+)_', name)
    if s.group('type') == 'SH':
        type = r'\perp'
    else:
        type = s.group('type')

    return r'$q_{' + type + ',' + comp + '}$'


def get_kinetic_heatflow(path, time):
    '''
       dict_x,dict_y = get_kinetic_heatflow(path,time)
    '''
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']

    dict_fxX = cf.load_dict(path, fprefix, 'fxX', time)
    fxX = dict_fxX['mat']

    dict_fxY = cf.load_dict(path, fprefix, 'fxY', time)
    fxY = dict_fxY['mat']

    dict_fyX = cf.load_dict(path, fprefix, 'fyX', time)
    fyX = dict_fyX['mat']

    dict_fyY = cf.load_dict(path, fprefix, 'fyY', time)
    fyY = dict_fyY['mat']

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']

    f1x_c = 0.5 * (fxX[:, :, 1:] + fxX[:, :, :-1])
    f1y_c = 0.5 * (fyY[:, 1:, :] + fyY[:, :-1, :])

    qx_c = 0.5 * (qxX[:, 1:] + qxX[:, :-1])
    qy_c = 0.5 * (qyY[1:, :] + qyY[:-1, :])

    vmomx_0, vmomy_0 = vmom_f1(0, v_grid, f1x_c, f1y_c)
    vmomx_3, vmomy_3 = vmom_f1(3, v_grid, f1x_c, f1y_c)
    vmomx_5, vmomy_5 = vmom_f1(5, v_grid, f1x_c, f1y_c)

    jx = -vmomx_0
    jy = -vmomy_0

    mom0_3 = vmom_fo(3, v_grid, fo)
    mom0_5 = vmom_fo(5, v_grid, fo)
    mom0_7 = vmom_fo(7, v_grid, fo)

    ##print 'np.shape(fxX) = ', np.shape(fxX), np.shape(fyY), np.shape(mom0_3), np.shape(mom0_5)
    ##print ' shapes = ',np.shape(vmomx_3),np.shape(vmomy_3),np.shape(vmomx_5),np.shape(vmomy_5)
    ny, nx = np.shape(vmomy_5)
    ##print ' shape(Z2ni) = ', np.shape(Z2ni)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    ##print 'np.sort(Z2ni) = ', np.sort(Z2ni)[0,0], np.sort(fo)
    ##print ' np.shape(Z2ni) = ', np.shape(Z2ni)
    rA = (Z2ni * 2.0)**-1

    vqx = rA * (vmomx_5 - (4.0 / 3.0) * (mom0_5 / mom0_3) * (vmomx_3))
    vqy = rA * (vmomy_5 - (4.0 / 3.0) * (mom0_5 / mom0_3) * (vmomy_3))

    vqBz_x = -1.0 * vqy * Bz
    vqBz_y = -1.0 * -1.0 * vqx * Bz
    #---------------
    #print np.shape(mom0_7), np.shape(mom0_7)[:], np.shape(mom0_7)[:][0]
    ny, nx = np.shape(mom0_7)

    x_grid = cf.trim_array(x_grid, nx, [1])
    y_grid = cf.trim_array(y_grid, ny, [1])

    ##print ' np.shape(x_Grid) = ', np.shape(x_grid),np.shape(y_grid), np.shape(mom0_7)
    dymom7, dxmom7 = cf.get_grad(y_grid, x_grid, mom0_7)
    dymom5, dxmom5 = cf.get_grad(y_grid, x_grid, mom0_5)
    dymom7, dxmom7 = -1.0 * dymom7, -1.0 * dxmom7
    dymom5, dxmom5 = -1.0 * dymom5, -1.0 * dxmom5

    ny, nx = np.shape(dxmom7)
    eta_star = Z2ni
    eta = Z2ni / (2.0 * mom0_3)    #Z2ni/(2vmom_0)
    #--- trim arrays
    ##print ' qx_c = ', np.shape(qx_c), np.shape(qy_c), np.shape(vqBz_y)
    array_list = [rA, eta, jx, jy, mom0_7, mom0_5, mom0_3, vqBz_x, vqBz_y]
    ##print ' nx = ',nx, 'ny = ', ny, np.shape(dxmom7)
    qx_c = cf.trim_array(np.transpose(qx_c), ny, nx)
    qy_c = cf.trim_array(np.transpose(qy_c), ny, nx)
    ##print ' qx_c = ', np.shape(qx_c), np.shape(qy_c)
    #cf.trim_array(array,ny,nx)
    rA, eta, jx, jy, mom0_7, mom0_5, mom0_3, vqBz_x, vqBz_y = map(
        lambda array: cf.trim_array(array, ny, nx), array_list)

    resist_x = rA * (8.0 / 3.0) * eta * jx * mom0_5
    resist_y = rA * (8.0 / 3.0) * eta * jy * mom0_5

    grad1_x = rA * (1.0 / 3.0) * dxmom7
    grad1_y = rA * (1.0 / 3.0) * dymom7

    grad2_x = -(4.0 / 9.0) * rA * (dxmom5) * (mom0_5 / mom0_3)
    grad2_y = -(4.0 / 9.0) * rA * (dymom5) * (mom0_5 / mom0_3)
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

    return dict_x, dict_y


def get_avgx(mat, ax=1):
    '''
        ax  = [n,..,2,1,0]
        ax = 1 x axis
        ax = 0 y axis
    '''

    avg = np.average(mat, axis=ax)
    return avg


#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
def get_Nernst_ratio(path, time, sep_var=False):
    '''
    data_x,data_y = get_Nernst_ratio(path,time)
    '''
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #------ NERNST -----------------------------------------------------
    v_nx, v_ny = get_v_N(grid, rho_vyx, Bz_vyx, fo)
    v_nx_classical, v_ny_classical, v_nx_hf, v_ny_hf = q_mod.get_vN_from_path(path, fprefix, time)
    if np.any(np.isnan(v_nx_classical)):
        print ' shape v_nx_classical = ', np.shape(v_nx_classical)
        print ' np.where(np.isnan(v_nx_classical))', np.where(np.isnan(v_nx_classical))
        print ' np.where(np.isnan(v_ny_classical))', np.where(np.isnan(v_ny_classical))

        sys.exit()
    v_nx_kinetic = np.transpose(v_nx)[:x_lim, :]
    v_ny_kinetic = np.transpose(v_ny)[:x_lim, :]

    #data_y = v_ny_kinetic/v_ny_classical[:x_lim,:]
    #data_y = np.where(v_ny_classical[:x_lim,:]==0.0,1.0,data_y)
    #vmin,vmax = -1.5,1.5
    #data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)

    if sep_var:
        print '--- Nernst separate ---'
        data_x = get_combined_dict(v_nx_classical[:x_lim, :], v_nx_kinetic)
        data_y = get_combined_dict(v_ny_classical[:x_lim, :], v_ny_kinetic)
    else:
        print ' --- nernst ratios ------'
        data_x = get_ratio_lim(v_nx_classical[:x_lim, :], v_nx_kinetic, vmax=1.5)
        data_y = get_ratio_lim(v_ny_classical[:x_lim, :], v_ny_kinetic, vmax=1.5)
    #data_x = v_nx_kinetic/v_nx_classical[:x_lim,:]
    #data_x = np.where(v_nx_classical[:x_lim,:]==0.0,1.0,data_x)
    #vmin,vmax = -1.0,1.0
    #data_x = np.where(np.abs(data_x)>vmax,np.sign(data_x)*vmax,data_x)
    return data_x, data_y


#-----------------------------------------------------------------------
def get_ratio_lim(q_yc, q_yk, vmax=1.5):
    print ' np.max(q_yc) min = ', np.min(np.abs(q_yc)), np.max(np.abs(q_yc))
    ratio_factor = np.abs(q_yk / q_yc)
    ratio_factor = np.where(np.isnan(q_yc), 0.0, ratio_factor)
    ratio_factor = np.where(np.abs(q_yc) <= (1e-10) * np.max(np.abs(q_yc)), 0.0, ratio_factor)
    #ratio_factor = np.where(np.abs(ratio_factor)>=10.0,np.sign(ratio_factor)*10.0,ratio_factor)
    data_y = ratio_factor
    #data_y = (q_yc - q_yk)/np.abs(q_yc)

    #data_y = np.where(q_yc==0.0,1.0,data_y)
    #data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)
    return data_y


#-----------------------------------------------------------------------
def get_combined_dict(q_yc, q_yk):
    comb_dict = {}
    comb_dict['classical'] = q_yc
    comb_dict['kinetic'] = q_yk
    return comb_dict


#-----------------------------------------------------------------------


def get_q_ratio(path, time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #heat flow ratios

    #--- kinetic
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)
    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[:x_lim]

    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
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
    q_x_VFP = (q_xX[:, 1:, :] + q_xX[:, :-1, :]) * 0.5
    q_y_VFP = (q_yY[:, :, 1:] + q_yY[:, :, :-1]) * 0.5
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
    q_yc = q_RL_y[time_int, :x_lim, :]
    rat_RL_y = get_ratio_lim(q_yc, q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int, :x_lim, :]
    rat_RL_x = get_ratio_lim(q_xc, q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int, :x_lim, :]
    rat_SH_y = get_ratio_lim(q_yc, q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int, :x_lim, :]
    rat_SH_x = get_ratio_lim(q_xc, q_xk)
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int, :x_lim, :]
    q_xc = q_x_B[time_int, :x_lim, :]
    rat_tot_x = get_ratio_lim(q_xc, q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int, :x_lim, :]
    q_yc = q_y_B[time_int, :x_lim, :]
    rat_tot_y = get_ratio_lim(q_yc, q_yk)

    dict_ratio = {}
    dict_ratio['SH x'] = rat_SH_x    #q_SH_x[time_int,:x_lim,:]
    dict_ratio['SH y'] = rat_SH_y    #q_SH_y[time_int,:x_lim,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y

    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y
    return dict_ratio


def extract_q(path, time, x_limit=c_index):
    '''
    dict_kinetic = extract_q(path,time)
    '''
    time_int = int(time)
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #heat flow ratios

    #--- kinetic
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)

    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']

    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[:x_limit]
        dict[name]['label'] = convert_name_to_label(name)
    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
    #--- classical
    '''
    
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
    q_yc = q_RL_y[time_int,:x_lim,:]   
    rat_RL_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int,:x_lim,:]   
    rat_RL_x = get_ratio_lim(q_xc,q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int,:x_lim,:]    
    rat_SH_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int,:x_lim,:]    
    rat_SH_x = get_ratio_lim(q_xc,q_xk)
    #--------
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int,:x_lim,:]    
    rat_SH_y = get_ratio_lim(q_yc,q_yk)

    q_xk = q_E_x
    q_xc = q_SH_x[time_int,:x_lim,:]    
    rat_SH_x = get_ratio_lim(q_xc,q_xk)
    
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int,:x_lim,:]
    q_xc = q_x_B[time_int,:x_lim,:]    
    rat_tot_x = get_ratio_lim(q_xc,q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int,:x_lim,:]
    q_yc = q_y_B[time_int,:x_lim,:]    
    rat_tot_y = get_ratio_lim(q_yc,q_yk)
    
    dict_ratio = {}
    dict_ratio['SH x'] = rat_SH_x#q_SH_x[time_int,:x_lim,:]
    dict_ratio['SH y'] = rat_SH_y#q_SH_y[time_int,:x_lim,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y
    
    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y
    return dict_ratio
    '''
    return dict


def get_q_individ(path, time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    fprefix = fpre(path)
    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    ###print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    ###print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    ##print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #heat flow ratios

    #--- kinetic
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)
    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[:x_lim]

    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
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
    q_x_VFP = (q_xX[:, 1:, :] + q_xX[:, :-1, :]) * 0.5
    q_y_VFP = (q_yY[:, :, 1:] + q_yY[:, :, :-1]) * 0.5
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
    q_yc = q_RL_y[time_int, :x_lim, :]
    rat_RL_y = get_combined_dict(q_yc, q_yk)

    q_xk = q_KRL_x
    q_xc = q_RL_x[time_int, :x_lim, :]
    rat_RL_x = get_combined_dict(q_xc, q_xk)
    #---- SH
    q_yk = q_KSH_y
    q_yc = q_SH_y[time_int, :x_lim, :]
    rat_SH_y = get_combined_dict(q_yc, q_yk)

    q_xk = q_KSH_x
    q_xc = q_SH_x[time_int, :x_lim, :]
    rat_SH_x = get_combined_dict(q_xc, q_xk)
    print '--- now doing the totals ---- x'
    #--- tot x
    q_xk = q_x_VFP[time_int, :x_lim, :]
    q_xc = q_x_B[time_int, :x_lim, :]
    rat_tot_x = get_combined_dict(q_xc, q_xk)
    #---- tot y
    print ' ---- now doing the totals ---y'
    q_yk = q_y_VFP[time_int, :x_lim, :]
    q_yc = q_y_B[time_int, :x_lim, :]
    rat_tot_y = get_combined_dict(q_yc, q_yk)

    dict_ratio = {}
    dict_ratio['SH x'] = rat_SH_x    #q_SH_x[time_int,:x_lim,:]
    dict_ratio['SH y'] = rat_SH_y    #q_SH_y[time_int,:x_lim,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y

    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y
    return dict_ratio


def get_alpha_perp_path(path, time):
    '''
     alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time)
    '''
    #--- load in data ....
    time_int = int(time)
    fprefix = fpre(path)

    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']
    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])
    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne
    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor
    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']

    Z2ni = ne
    nv, ny, nx = np.shape(fo)
    # trim data to correct size...
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    alpha_perp_kinetic = get_alpha_perp_kinetic(grid, rho_vyx, Bz_vyx, fo)

    #---- perform the classical calculation
    T_data = np.transpose(cf.trim_array(Te, nx, ny))
    n_data = np.transpose(cf.trim_array(ne, nx, ny))
    Bz_data = Bz
    #---------------------------------------
    alpha_perp_classical = np.zeros((np.shape(T_data)))
    #------
    #print ' shapes = ', np.shape(T_data),np.shape(n_data), np.shape(Bz_data), np.shape(alpha_perp_kinetic)
    #for ix in range(1,len(T_data[0,:])-1):
    #   for iy in range(1,len(T_data[:,0])-1):
    for ix in range(nx):
        for iy in range(ny):
            Z2ni_loc = n_data[iy, ix]
            Te_loc = T_data[iy, ix]
            ne_loc = n_data[iy, ix]
            w_loc = Bz_data[iy, ix]
            #v_nx[iy-1,ix-1],v_ny[iy-1,ix-1] = get_v_N_classical(Z2ni,n_data[iy,ix],T_data[iy,ix],Bz_data[iy-1,ix-1],dxT[iy-1,ix-1],dyT[iy-1,ix-1])
            #v_nx_hf[iy-1,ix-1] = q_xc[iy-1,ix-1]/(2.5*n_data[iy,ix]*T_data[iy,ix])
            #v_ny_hf[iy-1,ix-1] = q_yc[iy-1,ix-1]/(2.5*n_data[iy,ix]*T_data[iy,ix])

            vte = (2.0 * Te_loc)**0.5
            rho = Z2ni_loc**-1
            F0 = maxw_dist(v_grid, vte)
            alpha_perp_classical[iy, ix] = get_alpha_perp(w_loc, rho, ne_loc, Te_loc, F0, v_grid)
    #----------
    #dict[name]['data'] = np.transpose(dict[name]['data'])[:x_lim]
    return alpha_perp_classical, alpha_perp_kinetic


def extract_ohms(path, time, x_limit=73):
    '''
    dict = extract_ohms(path,time)
    '''
    #data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time)
    #data_x2d_vN_k,data_y2d_vN_k = data_x2d_vN['kinetic'],data_y2d_vN['kinetic']
    #data_x2d_vN_c,data_y2d_vN_c = data_x2d_vN['classical'],data_y2d_vN['classical']
    ##classic_biermann,kinetic_biermann = kbier.get_kinetic_and_classicB(path,path,time)
    #E_P_classical,E_P_kinetic =kbier.get_kinetic_and_classic_E_pressure(path,fpre(path),time,xlim=73)
    #E_x,E_y = kbier.get_hallfield(grid,rho_vyx,Bz_vyx,jx_vyx,jy_vyx,fo)
    #E_TEwedgex,E_TEwedgey = kbier.get_thermoelectricEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)

    #E_P_c_x, E_P_c_y = E_P_classical['E_P_x'],E_P_classical['E_P_y']
    #E_P_k_x, E_P_k_y = E_P_kinetic['E_P_x'],E_P_kinetic['E_P_y']
    fprefix = fpre(path)
    dict_kinetic = kbier.get_kinetic_E(path, fprefix, time, xlim=73)

    ##data_x2d_vN,data_y2d_vN =
    ##alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time):

    dict = {}
    dict['E_betawedge_x'] = {}
    dict['E_betawedge_x']['data'] = dict_kinetic['E_TEwedge_x']
    dict['E_betawedge_x']['label'] = r'$E_{\beta_{\wedge},x}$'

    dict['E_betawedge_y'] = {}
    dict['E_betawedge_y']['data'] = dict_kinetic['E_TEwedge_y']
    dict['E_betawedge_y']['label'] = r'$E_{\beta_{\wedge},y}$'

    dict['E_hall_x'] = {}
    dict['E_hall_x']['data'] = dict_kinetic['E_hall_x']
    dict['E_hall_x']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},x}$'

    dict['E_hall_y'] = {}
    dict['E_hall_y']['data'] = dict_kinetic['E_hall_y']
    dict['E_hall_y']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},y}$'

    dict['E_P_x'] = {}
    dict['E_P_x']['data'] = dict_kinetic['E_P_x']
    dict['E_P_x']['label'] = r'$E_{\nabla P,x}$'

    dict['E_P_y'] = {}
    dict['E_P_y']['data'] = dict_kinetic['E_P_y']
    dict['E_P_y']['label'] = r'$E_{\nabla P,y}$'

    dict['alpha_jx'] = {}
    dict['alpha_jx']['data'] = dict_kinetic['alpha_jx']
    dict['alpha_jx']['label'] = r'$\alpha_{\perp}\mathbf{j}_x$'

    dict['alpha_jy'] = {}
    dict['alpha_jy']['data'] = dict_kinetic['alpha_jy']
    dict['alpha_jy']['label'] = r'$\alpha_{\perp}\mathbf{j}_y$'

    func_l = lambda data: np.transpose(data)    #[:,:]
    for key in dict.keys():
        dict[key]['data'] = func_l(dict[key]['data'])

    return dict


#-----------------------------------------------------------------------
def extract_ohms(path, time, x_limit=73):
    '''
    dict = extract_heatflow(path,time)
        wrapper for get_kinetic_q(path,time)
    '''
    #data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time)
    #data_x2d_vN_k,data_y2d_vN_k = data_x2d_vN['kinetic'],data_y2d_vN['kinetic']
    #data_x2d_vN_c,data_y2d_vN_c = data_x2d_vN['classical'],data_y2d_vN['classical']
    ##classic_biermann,kinetic_biermann = kbier.get_kinetic_and_classicB(path,path,time)
    #E_P_classical,E_P_kinetic =kbier.get_kinetic_and_classic_E_pressure(path,fpre(path),time,xlim=73)
    #E_x,E_y = kbier.get_hallfield(grid,rho_vyx,Bz_vyx,jx_vyx,jy_vyx,fo)
    #E_TEwedgex,E_TEwedgey = kbier.get_thermoelectricEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)

    #E_P_c_x, E_P_c_y = E_P_classical['E_P_x'],E_P_classical['E_P_y']
    #E_P_k_x, E_P_k_y = E_P_kinetic['E_P_x'],E_P_kinetic['E_P_y']
    fprefix = fpre(path)
    dict_kinetic = kbier.get_kinetic_q(path, fprefix, time)

    ##data_x2d_vN,data_y2d_vN =
    ##alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time):

    dict = {}
    dict['q_SH_x'] = {}
    dict['q_SH_x']['data'] = dict_kinetic['q_SH_x']
    dict['q_SH_x']['label'] = r'$q_{\beta_{\wedge},x}$'

    dict['E_betawedge_y'] = {}
    dict['E_betawedge_y']['data'] = dict_kinetic['E_TEwedge_y']
    dict['E_betawedge_y']['label'] = r'$E_{\beta_{\wedge},y}$'

    dict['E_hall_x'] = {}
    dict['E_hall_x']['data'] = dict_kinetic['E_hall_x']
    dict['E_hall_x']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},x}$'

    dict['E_hall_y'] = {}
    dict['E_hall_y']['data'] = dict_kinetic['E_hall_y']
    dict['E_hall_y']['label'] = r'$E_{\mathbf{j}\times \mathbf{B},y}$'

    dict['E_P_x'] = {}
    dict['E_P_x']['data'] = dict_kinetic['E_P_x']
    dict['E_P_x']['label'] = r'$E_{\nabla P,x}$'

    dict['E_P_y'] = {}
    dict['E_P_y']['data'] = dict_kinetic['E_P_y']
    dict['E_P_y']['label'] = r'$E_{\nabla P,y}$'

    dict['alpha_jx'] = {}
    dict['alpha_jx']['data'] = dict_kinetic['alpha_jx']
    dict['alpha_jx']['label'] = r'$\alpha_{\perp}\mathbf{j}_x$'

    dict['alpha_jy'] = {}
    dict['alpha_jy']['data'] = dict_kinetic['alpha_jy']
    dict['alpha_jy']['label'] = r'$\alpha_{\perp}\mathbf{j}_y$'

    func_l = lambda data: np.transpose(data)    #[:,:]
    for key in dict.keys():
        dict[key]['data'] = func_l(dict[key]['data'])

    return dict


def fpre(path):
    return path.split('/')[-1]


def plot_2D(ax,
            data,
            label,
            middle=None,
            colormap='jet',
            limits=None,
            xlabel=None,
            ylabel=None,
            clabel=None,
            fmt='%0.1e'):
    #fmt = ticker.FuncFormatter(fmt)
    ##fmt = ScalarFormatter(useMathText=True)
    fmt = '%i'
    sample = np.average(np.abs(data))
    power = cf.extract_power(sample)
    data *= (10**-power)

    lab_ad = r' $[10^{' + str(power) + '} q_0] $'
    label = label + lab_ad
    print 'label = ', label, (10**-power)
    if middle != None:
        try:
            norm = cf.MidPointNorm(middle)
        except ValueError:
            norm = None
    else:
        norm = None

    if limits:
        try:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm, extent=limits)
            cbar = plt.colorbar(im, ax=ax, aspect='auto', norm=norm, format=fmt, label=clabel)
        except ValueError:
            norm = None
            #print 'setting norm to None;'
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm, extent=limits)
            cbar = plt.colorbar(im, ax=ax, aspect='auto', norm=norm, format=fmt, label=clabel)

    else:
        try:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm)
            cbar = plt.colorbar(im, ax=ax, aspect='auto', norm=norm, format=fmt, label=clabel)
        except ValueError:
            im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm)
            cbar = plt.colorbar(im, ax=ax, aspect='auto', norm=norm, format=fmt, label=clabel)

    #setp(im.get_title(),fontsize=10)
    ##cbar.formatter.set_scientific(True)
    #ax.offsetText
    #cbar.yaxis.set_major_formatter(ScalarFormatter())
    #cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    #cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    #cbar.formatter.set_powerlimits((0,0))
    ax.set_title(label, fontsize=12)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_2D_irregulargrid(ax,
                          xgrid,
                          ygrid,
                          data,
                          label,
                          middle=None,
                          colormap='jet',
                          limits=None,
                          xlabel=None,
                          ylabel=None,
                          clabel=None,
                          fmt='%0.1e'):
    #fmt = ticker.FuncFormatter(fmt)
    ##fmt = ScalarFormatter(useMathText=True)
    fmt = '%i'
    sample = np.average(np.abs(data))
    power = cf.extract_power(sample)
    data *= (10**-power)

    lab_ad = r' $[10^{' + str(power) + '} q_0] $'
    label = label + lab_ad
    print 'label = ', label, (10**-power)
    if middle != None:
        try:
            norm = cf.MidPointNorm(middle)
        except ValueError:
            norm = None
    else:
        norm = None
    print ' shapes = ', np.shape(xgrid), np.shape(ygrid), np.shape(data)

    xs0 = ygrid
    ys0 = xgrid
    XS, YS = np.meshgrid(xs0, ys0)
    xs0 = XS.flatten()
    ys0 = YS.flatten()

    zs0 = data.flatten()
    print ' shapes = ', np.shape(xs0), np.shape(ys0), np.shape(zs0)
    N = 30j
    M = 500j
    #extent = (-np.pi/2,np.pi/2,0,3.5)
    extent = (limits)
    print('extent = ', extent)

    xs, ys = np.mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:M]
    resampled = griddata(xs0, ys0, zs0, xs, ys, interp='linear')
    print(' shape resampled =', np.shape(resampled))
    extent_rev = [extent[0], extent[1], extent[3], extent[2]]
    im = ax.imshow(resampled.T[:, :], aspect='auto', norm=norm, extent=extent_rev)
    cbar = plt.colorbar(im, ax=ax, aspect='auto', norm=norm, format=fmt, label=clabel)
    #ax.set_ylim(-10.0,20.0)
    #plt.plot(xs0, ys0, "r.")
    #plt.plot(xs, ys, "b.")
    #plt.title("imshow for irregularly spaced data using griddata")
    '''
    if limits:
        try:
            im = ax.imshow(resampled.T,aspect='auto',cmap=colormap,norm=norm,extent=limits)
            cbar = plt.colorbar(im,ax=ax,aspect='auto',norm=norm,format= fmt,label=clabel)
        except ValueError:
            norm=None
            #print 'setting norm to None;'
            im = ax.imshow(resampled.T,aspect='auto',cmap=colormap,norm=norm,extent=limits)
            cbar = plt.colorbar(im,ax=ax,aspect='auto',norm=norm,format= fmt,label=clabel)
            
    else:
        try:
            im = ax.imshow(resampled.T,aspect='auto',cmap=colormap,norm=norm)
            cbar = plt.colorbar(im,ax=ax,aspect='auto',norm=norm,format= fmt,label=clabel)
        except ValueError:
            im = ax.imshow(resampled.T,aspect='auto',cmap=colormap,norm=norm)
            cbar = plt.colorbar(im,ax=ax,aspect='auto',norm=norm,format= fmt,label=clabel)
    '''
    #setp(im.get_title(),fontsize=10)
    ##cbar.formatter.set_scientific(True)
    #ax.offsetText
    #cbar.yaxis.set_major_formatter(ScalarFormatter())
    #cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

    #cbar.ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    #cbar.formatter.set_powerlimits((0,0))
    ax.set_title(label, fontsize=12)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


if __name__ == "__main__":
    '''
        dxT = T_data[1:] - T_data[:-1]/(x_grid[1:] - x_grid[:-1])
    '''
    fprefix = fpre(path)

    dict_Bz = cf.load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = cf.load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = cf.load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']
    #print 'np.shape(jxX) = ', np.shape(jxX), np.shape(jyY)
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = cf.load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = cf.load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    #print 'np.shape(qxX) = ', np.shape(qxX), np.shape(qyY)
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])

    dict_wt = cf.load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']

    dict_ne = cf.load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    Z2ni = ne

    dict_te = cf.load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']
    x_grid_SI = x_grid * xstep_factor
    y_grid_SI = y_grid * xstep_factor
    time_col = dict_te['time'] * tstep_factor

    dict_fo = cf.load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    #print ' -------------------- SHAPES --------------------------------'

    grid = dict_fo
    nv, ny, nx = np.shape(fo)
    Z2ni = np.transpose(cf.trim_array(Z2ni, nx, ny))
    Bz = np.transpose(cf.trim_array(Bz, nx, ny))
    jx = np.transpose(cf.trim_array(jx_c, nx, ny))
    jy = np.transpose(cf.trim_array(jy_c, nx, ny))
    # #print ' np.shapae(jx_c) = ', np.shape(jx_c),np.shape(jy_c)
    # #print ' np.shape(Z2ni) = ', np.shape(Z2ni),' np.shape(Bz) = ', np.shape(Bz)

    rA = (Z2ni)**-1

    #===============
    grid = dict_fo
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #------ NERNST -----------------------------------------------------
    v_nx, v_ny = get_v_N(grid, rho_vyx, Bz_vyx, fo)    # kinetic Nernst
    v_nx_classical, v_ny_classical, v_nx_hf, v_ny_hf = q_mod.get_vN_from_path(path, fprefix, time)

    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)

    #print ' dict keys = ', dict.keys()
    if False:
        #fig = plt.figure()
        fig = fprl.newfig_generic_2yscale(width, scale_width=1.0, scale_ratio=1.0, clearon=True)

        ny_plots, nx_plots = 1, 2
        var_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']

        var_list = ['q_SH_x', 'q_SH_y']
        q_SH_x, q_SH_y, q_RL_x, q_RL_y, q_TE_x, q_TE_y, q_E_x, q_E_y = dict['q_SH_x']['data'], dict[
            'q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict['q_TE_x'][
                'data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
        qx_sum = q_SH_x + q_RL_x + q_TE_x + q_E_x
        qy_sum = q_SH_y + q_RL_y + q_TE_y + q_E_y

    limits = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[x_lim], x_grid_SI[0]]
    lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    lab_dict['lims'] = limits    #[y_grid_SI[0],y_grid_SI[-1],x_grid_SI[c_index],x_grid_SI[0]]
    lab_dict['colormap'] = 'RdBu_r'
    lab_dict['xlab'], lab_dict['ylab'] = xlab, ylab
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

    #cf.produce_2D_fig(x_grid_SI,data,lab_dict,lineout_list=[])
    #print ' path = ', path
    #plt.show()

    nt = 15

    #------- plot the ratios now
    if False:
        fig = fprl.newfig_generic(1.0)

        plt.clf()
        #gs1 = GS.GridSpec(2,3,width_ratios=[1.0,1.0,1.0])
        #pl
        [ax1, ax2,
         ax3] = [fig.add_subplot(2, 3, 1),
                 fig.add_subplot(2, 3, 2),
                 fig.add_subplot(2, 3, 3)]
        [ax4, ax5,
         ax6] = [fig.add_subplot(2, 3, 4),
                 fig.add_subplot(2, 3, 5),
                 fig.add_subplot(2, 3, 6)]
        ax_list = [ax1, ax2, ax3, ax4, ax5, ax6]
        map(fprl.fix_ax, [ax1, ax2, ax3, ax4, ax5, ax6])
        #plt.subplots_adjust(hspace=0.1,left=0.2,right=0.95,wspace=0.68)
        plt.subplots_adjust(hspace=0.1, left=0.1, right=0.95, wspace=0.5)
        t_list = [3, 8, 12, 15]
        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))

        c_list = iter(fcmap(np.linspace(0, 1, len(t_list))))
        log_on = True
        lab_list = []
        plot_list = []
        for tt in t_list:
            c = next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt

            dict_ratio = get_q_ratio(path, time)
            data_x2d_SH, data_y2d_SH, prelab_SH = dict_ratio['SH x'], dict_ratio['SH y'], 'q_{\perp'
            data_x2d_RL, data_y2d_RL, prelab_RL = dict_ratio['RL x'], dict_ratio['RL y'], 'q_{RL'
            data_x2d_vN, data_y2d_vN = get_Nernst_ratio(path, time)
            prelab_vN = 'v_{N'
            #prelab = 'v_{N'
            data_xSH = get_avgx(data_x2d_SH)
            data_ySH = get_avgx(data_y2d_SH)
            data_xRL = get_avgx(data_x2d_RL)
            data_yRL = get_avgx(data_y2d_RL)
            data_xvN = get_avgx(data_x2d_vN)
            data_yvN = get_avgx(data_y2d_vN)

            #p1, =ax1.plot(x_grid_SI[:x_lim],data_x,c=c)
            if log_on:
                p1, = ax1.semilogy(x_grid_SI[:x_lim], data_xSH, c=c)
                p2, = ax2.semilogy(x_grid_SI[:x_lim], data_xRL, c=c)
                p3, = ax3.semilogy(x_grid_SI[:x_lim], data_xvN, c=c)

                p4, = ax4.semilogy(x_grid_SI[:x_lim], data_ySH, c=c)
                p5, = ax5.semilogy(x_grid_SI[:x_lim], data_yRL, c=c)
                p6, = ax6.semilogy(x_grid_SI[:x_lim], data_yvN, c=c)

            else:
                p1, = ax1.plot(x_grid_SI[:x_lim], data_xSH, c=c)
                p2, = ax2.plot(x_grid_SI[:x_lim], data_xRL, c=c)
                p3, = ax3.plot(x_grid_SI[:x_lim], data_xvN, c=c)

                p4, = ax4.plot(x_grid_SI[:x_lim], data_ySH, c=c)
                p5, = ax5.plot(x_grid_SI[:x_lim], data_yRL, c=c)
                p6, = ax6.plot(x_grid_SI[:x_lim], data_yvN, c=c)

            #--- legend stuff
            dict_wt = cf.load_dict(path, fprefix, 'wt', time)
            t_l = dict_wt['time'] * tstep_factor
            t_lab = '%3i \si{ps}' % t_l
            t_lab = '%3i' % t_l

            plot_list.append(p1)
            lab_list.append(t_lab)

        xtit = lambda ax, name: ax.set_ylabel('$ %s,x,K}/ %s,x,B}$' % (name, name))
        ytit = lambda ax, name: ax.set_ylabel('$ %s,y,K}/ %s,y,B}$' % (name, name))

        def clearxaxis(ax):
            ax.set_xlabel('')
            ax.set_xticklabels([])

        map(xtit, [ax1, ax2, ax3], [prelab_SH, prelab_RL, prelab_vN])
        map(ytit, [ax4, ax5, ax6], [prelab_SH, prelab_RL, prelab_vN])
        map(clearxaxis, [ax1, ax2, ax3])
        map(lambda axy: axy.set_xlabel(xlab), [ax4, ax5, ax6])
        xlim = ax1.get_xlim()
        map(lambda axy: axy.set_xlim(10.0, xlim[-1]), ax_list)

        # - qSHx, vNx
        def sort_y(axy):
            axy.set_ylim(0.1, 2.0)
            #axy.yaxis.labelpad= 4
            #axy.yaxis.set_ticklabels([0.2,1.0,2.0])

        map(lambda axy: sort_y(axy), [ax1, ax5, ax3])
        # - qSHy, vNy
        map(lambda axy: axy.set_ylim(0.5e-2, 5e2), [ax2, ax4, ax6])

        #ax1.set_ylabel(xtit)
        #ax2.set_ylabel(ytit)
        #ax1.set_xlabel(xlab)
        #ax2.set_xlabel(xlab)
        ax3.legend(plot_list, lab_list, title=r't $\si{ps}$')
        fig.savefig('test.pdf')
        fprl.savefig_thesis(ratios_savename)
        #plt.show()
        sys.exit()

    #------- plot the ratios now ---------------------
    if False:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)
        t_list = [3, 8, 12, 15]
        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        c_list = iter(cm.plasma(np.linspace(0, 1, len(t_list))))

        lab_list = []
        plot_list = []
        for tt in t_list:
            c = next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt

            dict_ratio = get_q_ratio(path, time)
            data_x2d, data_y2d, prelab = dict_ratio['SH x'], dict_ratio['SH y'], 'q_{SH'
            data_x2d_SH, data_y2d_SH, prelab_SH = dict_ratio['SH x'], dict_ratio['SH y'], 'q_{SH'

            data_x2d, data_y2d, prelab = dict_ratio['RL x'], dict_ratio['RL y'], 'q_{RL'
            #data_x2d,data_y2d = get_Nernst_ratio(path,time)
            #prelab = 'v_{N'
            data_xSH = get_avgx(data_x2d_SH)
            data_ySH = get_avgx(data_y2d_SH)

            data_x = get_avgx(data_x2d)
            data_y = get_avgx(data_y2d)
            if log_on:
                p1, = ax1.semilogy(x_grid_SI[:x_lim], data_x, c=c)
                p1, = ax1.semilogy(x_grid_SI[:x_lim], data_xSH, c=c, linestyle='--')

                ax2.semilogy(x_grid_SI[:x_lim], data_y, c=c)
                ax2.semilogy(x_grid_SI[:x_lim], data_xSH, c=c, linestyle='--')

            else:
                p1, = ax1.plot(x_grid_SI[:x_lim], data_x, c=c)
                p1, = ax1.plot(x_grid_SI[:x_lim], data_xSH, c=c, linestyle='--')

                ax2.plot(x_grid_SI[:x_lim], data_y, c=c)
                ax2.plot(x_grid_SI[:x_lim], data_xSH, c=c, linestyle='--')

            #--- legend stuff
            dict_wt = cf.load_dict(path, fprefix, 'wt', time)
            t_l = dict_wt['time'] * tstep_factor
            t_lab = '%3i ps' % t_l
            plot_list.append(p1)
            lab_list.append(t_lab)
        xtit = '$' + prelab + ',x,K}/' + prelab + ',x,B}$'
        ytit = '$' + prelab + ',y,K}/' + prelab + ',y,B}$'

        ax1.set_title(xtit)
        ax2.set_title(ytit)
        ax1.set_xlabel(xlab)
        ax2.set_xlabel(xlab)
        ax2.legend(plot_list, lab_list)
        plt.show()

    #-------- plot resistivitity

    if False:
        fig = plt.figure()
        ax1 = fig.add_subplot(1, 1, 1)
        #ax2 = fig.add_subplot(1,2,2)

        t_list = [3, 12, 15]
        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        c_list = iter(cm.plasma(np.linspace(0, 1, len(t_list))))
        lab_list = []
        plot_list = []

        for tt in t_list:
            c = next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt

            alpha_perp_classical, alpha_perp_kinetic = get_alpha_perp_path(path, time)
            print 'np.shape alpha_perp = ', np.shape(alpha_perp_classical), np.shape(
                alpha_perp_kinetic)
            ###sys.exit()

            data_x_clas = get_avgx(np.transpose(alpha_perp_classical))
            data_x_kin = get_avgx(np.transpose(alpha_perp_kinetic))
            print ' DATA x avg = ', np.shape(data_x_clas), np.shape(data_x_kin)
            #p1, =ax1.plot(x_grid_SI[:x_lim],data_x_clas[:x_lim],c=c,linestyle=':')
            #p2, =ax1.plot(x_grid_SI[:x_lim],data_x_kin[:x_lim],c=c)
            p2, = ax1.plot(x_grid_SI[:x_lim], data_x_kin[:x_lim] / data_x_clas[:x_lim], c=c)

            #--- legend stuff
            dict_wt = cf.load_dict(path, fprefix, 'wt', time)
            t_l = dict_wt['time'] * tstep_factor
            t_lab = '%3i ps' % t_l
            plot_list.append(p2)
            lab_list.append(t_lab)
        #xtit = '$' + prelab+ ',x,K}/' + prelab + ',x,B}$'
        ytit = r'$ \alpha_{\perp}$'

        ax1.set_title(ytit)
        #ax2.set_title(ytit)
        ax1.set_xlabel(xlab)
        #ax2.set_xlabel(xlab)
        ax1.legend(plot_list, lab_list)
        plt.show()
    #----------------------

    #------- plot the ratios now ---------------------
    if False:
        fig = fprl.newfig_generic(1.0)
        ax1 = fig.add_subplot(1, 3, 1)
        ax2 = fig.add_subplot(1, 3, 2)
        ax3 = fig.add_subplot(1, 3, 3)
        t_list = [4, 8, 12, 15]
        map(fprl.fix_ax, [ax1, ax2, ax3])

        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        #map(fprl.fix_ax,[ax1,ax2])

        c_list = iter(cm.plasma(np.linspace(0, 1, len(t_list))))

        lab_list = []
        plot_list = []
        for tt in t_list:
            c = next(c_list)
            print 'tt = ', tt
            time = '%02i' % tt

            dict_ratio = get_q_ratio(path, time)
            data_x2d_RL, data_y2d_RL, prelab_RL = dict_ratio['RL x'], dict_ratio['RL y'], 'q_{RL'
            data_x2d_SH, data_y2d_SH, prelab_SH = dict_ratio['SH x'], dict_ratio['SH y'], 'q_{SH'
            data_x2d_tot, data_y2d_tot, prelab_tot = dict_ratio['tot x'], dict_ratio[
                'tot y'], 'q_{SH'
            data_x2d_vN, data_y2d_vN = get_Nernst_ratio(path, time)
            #prelab = 'v_{N'
            data_x_SH = get_avgx(data_x2d_SH)
            data_y_SH = get_avgx(data_y2d_SH)
            data_x_RL = get_avgx(data_x2d_RL)
            data_y_RL = get_avgx(data_y2d_RL)
            data_x_vN = get_avgx(data_x2d_vN)
            data_y_vN = get_avgx(data_y2d_vN)
            data_x_tot = get_avgx(data_x2d_tot)
            data_y_tot = get_avgx(data_y2d_tot)

            if log_on:
                p1, = ax1.semilogy(x_grid_SI[:x_lim], data_x_tot, c=c)
                #p1, =ax1.semilogy(x_grid_SI[:x_lim],data_xSH,c=c,linestyle='--')
                ax2.semilogy(x_grid_SI[:x_lim], data_y_tot, c=c)
                ax3.semilogy(x_grid_SI[:x_lim], data_x_vN, c=c)

                #ax2.semilogy(x_grid_SI[:x_lim],data_xSH,c=c,linestyle='--')

            else:
                p1, = ax1.plot(x_grid_SI[:x_lim], data_x_tot, c=c)
                ax2.plot(x_grid_SI[:x_lim], data_y_tot, c=c)
                ax3.plot(x_grid_SI[:x_lim], data_x_vN, c=c)
            #--- legend stuff
            dict_wt = cf.load_dict(path, fprefix, 'wt', time)
            t_l = dict_wt['time'] * tstep_factor
            t_lab = '%3i ps' % t_l
            plot_list.append(p1)
            lab_list.append(t_lab)

        #xtit = '$' + prelab+ ',x,K}/' + prelab + ',x,B}$'
        #ytit = '$' + prelab+ ',y,K}/' + prelab + ',y,B}$'
        xtit = 'ratios x'
        ytit = 'ratios y'
        ax1.set_title(xtit)
        ax2.set_title(ytit)
        ax1.set_xlabel(xlab)
        ax2.set_xlabel(xlab)
        ax2.legend(plot_list, lab_list)
        #plt.show()

    #------- plot the individual kinetic + classical heatflows + Nernst ---------------------
    if True:
        scale_ratio = conv_factors_cd5.double_plot_scale_ratio    # reducing this from one will reduce the height relative to the width
        fig = fprl.newfig_generic(1.0)    #
        print ' fig.get_size_inches() = ', fig.get_size_inches()
        fig.subplots_adjust(left=0.08, right=0.95, hspace=0.61,
                            top=0.92)    # adjust bottom so we can fit a legend underneath
        #fig.subplots_adjust(right=0.85) # adjust bottom so we can fit a legend underneath
        fig.subplots_adjust(wspace=0.75)

        #ax1 = fig.add_subplot(1,4,3)
        #ax2 = fig.add_subplot(1,4,4)
        t_list = [12]

        #c_list = iter(cm.jet(np.linspace(0,1,len(t_list))))
        c_q = 'r'
        c_n = 'b'
        lab_list = []
        plot_list = []
        yi = 10
        xi_min = 3

        lims = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[x_lim], x_grid_SI[0]]    #
        lims = [y_grid_SI[0], y_grid_SI[-1], x_grid_SI[0], x_grid_SI[-1]]    #
        #lims = [y_grid_SI[0],y_grid_SI[-1],-10.0,20.0]#
        lims = [y_grid_SI[1], y_grid_SI[-2], -10.0, 20.0]    #

        for tt in [13]:
            #c_q = c
            print 'tt = ', tt
            time = '%02i' % tt
            # ohm's law
            #dict_kinetic = extract_ohms(path,time,x_lim)
            # heat flow
            dict_kinetic = kbier.get_kinetic_q(path, time)
            ax2d = map(lambda ii: fig.add_subplot(2, 4, ii + 1), range(8))
            axli = iter(ax2d)
            print 'keys = ', dict_kinetic.keys()
            # ohm's law
            ##it_list = ['alpha_jx','alpha_jy', 'E_P_x','E_P_y','E_hall_x','E_hall_y','E_betawedge_x','E_betawedge_y']
            # heat flow
            it_list = ['q_SH_x', 'q_RL_x', 'q_TE_x', 'q_E_x', 'q_SH_y', 'q_RL_y', 'q_TE_y', 'q_E_y']

            for key in it_list:
                ax = next(axli)
                data = dict_kinetic[key]['data']
                label = dict_kinetic[key]['label']
                print ' key = ', key, label, np.shape(data)
                if it_list[-1] == 'y':
                    mid = 0.0
                else:
                    mid = None
                #plot_2D_irregulargrid(ax,data,label,middle=mid,colormap='RdBu_r',limits=lims,xlabel=ylab,ylabel=xlab,clabel=' ')
                plot_2D_irregulargrid(ax,
                                      x_grid_SI[1:-1],
                                      y_grid_SI[1:-1],
                                      data,
                                      label,
                                      middle=mid,
                                      colormap='RdBu_r',
                                      limits=lims,
                                      xlabel=ylab,
                                      ylabel=xlab,
                                      clabel=' ')

                ax.tick_params(which='both', direction='in')
                ax.title.set_fontsize(12)

                for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() +
                             ax.get_yticklabels()):
                    item.set_fontsize(10)
            #savepath = 'pics_transport/'
            savename = save_path + 'im_q_' + str(time) + '.pdf'
            plt.savefig(savename)
            plt.show()
            sep_var = True
            #data_x2d_vN,data_y2d_vN = get_Nernst_ratio(path,time,sep_var)
            #data_x2d_vN_k,data_y2d_vN_k = data_x2d_vN['kinetic'],data_y2d_vN['kinetic']
            #data_x2d_vN_c,data_y2d_vN_c = data_x2d_vN['classical'],data_y2d_vN['classical']

            #--- plot kinetic components
            #prelab = 'v_{N'

            #data_x_SH = (data_x2d_SH_k[xi_min:,yi])
            #data_y_SH = (data_y2d_SH_k[xi_min:,yi])
            #data_x_RL = (data_x2d_RL_k[xi_min:,yi])
            #data_y_RL = (data_y2d_RL_k[xi_min:,yi])
            #data_x_vN = (data_x2d_vN_k[xi_min:,yi])
            #data_y_vN = (data_y2d_vN_k[xi_min:,yi])
            #data_x_tot = (data_x2d_tot_k[xi_min:,yi])
            #data_y_tot = (data_y2d_tot_k[xi_min:,yi])

            marker_q = None
            marker_nernst = None
            #0-0000----------------------------
        #fprl.savefig(savename)
        print '----- saving as ... ', savename
        #plt.show()
