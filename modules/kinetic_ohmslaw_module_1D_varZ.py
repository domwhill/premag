'''

    Computes kinetic ohms law components
    
    formerly : 'kinetic_biermann_compare_lineouts.py'
    
    18/1/2019 - Adapted to allow for variable Z profiles
      functions changed but not tested: get_Nernst_ratio (dependency on q_SH_Te_tsteps) - probably obsolete anyway

    For full details of defn of integrals see Table 5.2 & Table 5.3 - D. Hill PhD thesis
    
    l 867
        Bz = np.transpose(cf.trim_array(Bz,nx,ny))<-- this may need reflecting in y direction...?

         dyT = np.transpose(cf.trim_array(dyT,nx,ny))*-1.0<-- check that this minus 1 needs to be here
    
'''

import numpy as np
import matplotlib.pyplot as plt
import sys, re, os
sys.path.extend(["./"])
import modules.chfoil_module as cf

from pylab import *

nv = 100

Y1_weight = 4.0 * np.pi / 3.0
Y0_weight = 4.0 * np.pi

#-------------------------------------->
# get normalisations

tstep_factor = 1.0
#<-------------------------------------


def convert_name_to_label(name):
    comp = name[-1]
    s = re.search('q_(?P<type>\w+)_', name)
    if s.group('type') == 'SH':
        type = r'\perp'
    else:
        type = s.group('type')

    return r'$q_{' + type + ',' + comp + '}$'


def fpre(path):
    return path.split('/')[-1]


def vc_to_vb(vc):
    vb = np.zeros((len(vc) + 1))
    vred = 0.5 * (vc[1:] + vc[:-1])
    vb[1:-1] = vred
    dv0 = vred[1] - vred[0]
    dv1 = vred[-1] - vred[-2]
    vb[0] = vred[0] - dv0
    vb[-1] = vred[-1] + dv1
    return vb


def calc_dv(vc):
    '''
        dv = calc_dv(vc)
    '''
    vb = vc_to_vb(vc)
    return vb[1:] - vb[:-1]


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


def get_v2dv(v_grid):
    dv = calc_dv(v_grid)
    v2dv = (v_grid**2) * dv
    return v2dv


def maxw_dist(v, vte):
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    return f_om


def maxw_dist_impact(v, ne, vte):
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    return f_om


def get_v_mom_m(v_grid, omega, F0, m):
    '''
        int_0^\inf dv V^{m+2}F0/(1 + omega^2*V^6)
    '''

    #ones_nv
    v2dv = get_v2dv(v_grid)
    ones_nv = np.ones((len(v_grid)))
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

    return delta


def get_fo_mom_int(m, rho_vyx, omega_vyx, fo, v_grid):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni
    # omega  = Bz/Z2ni = Bz*rho

    nv, ny, nx = np.shape(fo)

    v2dv = get_v2dv(v_grid)
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)
    ones_vyx = np.ones((nv, ny, nx))

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / (ones_vyx + (omega_vyx**2) * v_grid_vyx**6)
    mom = np.sum(prefactor * int, axis=0)
    return mom


def get_fo_mom_int2(m, rho_vyx, omega_vyx, fo, v_grid):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni

    nv, ny, nx = np.shape(fo)

    v2dv = get_v2dv(v_grid)
    ones_vyx = np.ones((nv, ny, nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / ((ones_vyx + (omega_vyx**2) * (v_grid_vyx**6))**2)
    mom = np.sum(prefactor * int, axis=0)
    return mom


def get_dfodv_int(n, rho_vyx, omega_vyx, fo, v_grid):
    # omega  = Bz/Z2ni = Bz*rho

    v_mom_nm3 = get_fo_mom_int(n - 3, rho_vyx, omega_vyx, fo, v_grid)
    v_mom2_np3 = get_fo_mom_int2(n + 3, rho_vyx, omega_vyx, fo, v_grid)
    mom_out = -1.0 * (n * v_mom_nm3 - 6.0 * (omega_vyx[0, :, :]**2) * v_mom2_np3
                     )    # unclear whether this minus one should actually be there or not
    return mom_out


def get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    m = 5
    mom_x = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho

    m = 6
    mom = get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I4 = get_I4_scalar(rho,omega,fo)
    '''

    m = 9
    mom = rho_vyx[0, :, :] * get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I3 =  get_I3_vec(rho,omega,gradfo)
    '''

    m = 8
    mom_x = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    m = 7
    mom_x = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K2_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    m = 8
    mom = 0.5 * get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_K3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    m = 10
    mom_x = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K4_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho

    m = 11
    mom = 0.5 * rho_vyx[0, :, :] * get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_hallfield(grid, rho_vyx, Bz_vyx, jx_vyx, jy_vyx, fo):
    '''
    # kinetic nernst
    E_x,E_y = get_hallfield(grid,rho_vyx,Bz_vyx,jx_vyx,jy_vyx,fo)
    '''
    v_grid = grid['v_grid']
    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)

    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    jx_yx = jx_vyx[0, :, :]
    jy_yx = jy_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    hallcoeff = (1.0 / Eta) * (I4 / I2)
    # E = -jxomega * hallcoeff
    ##E = omega x v
    #    |  i    j    k |    | -j_y B_z |
    #-1* |  jx   jy   0 |    |  j_x B_z |
    #    |  0    0    Bz|  = |    0     |

    E_x = -Bz_yx * jy_yx * hallcoeff
    E_y = Bz_yx * jx_yx * hallcoeff
    return E_x, E_y


def get_v_N(grid, rho_vyx, Bz_vyx, fo):
    '''
    # kinetic nernst
    v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    '''
    v_grid = grid['v_grid']

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)

    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    v_nx = (1.0 / Eta) * (I3_x - (I4 / I2) * I1_x)
    v_ny = (1.0 / Eta) * (I3_y - (I4 / I2) * I1_y)

    return v_nx, v_ny


def get_thermoelectricEfield_kinetic(grid, rho_vyx, Bz_vyx, fo):
    '''
        E_x,E_y = get_thermoelectricEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)
     '''
    v_nx, v_ny = get_v_N(grid, rho_vyx, Bz_vyx, fo)
    Bz_yx = Bz_vyx[0, :, :]
    ##E = omega x v
    #   i    j    k |    | -V_y B_z |
    #             Bz|    |  v_x B_z |
    #  v_x  v_y   0 |  = |    0     |

    E_x = -Bz_yx * v_ny
    E_y = Bz_yx * v_nx
    return E_x, E_y


def get_pressureEfield_kinetic(grid, rho_vyx, Bz_vyx, fo):
    '''
        ----- should get the kinetic biermann term ----->
       BierE_x,BierE_y =  get_pressureEfield_kinetic(grid,rho_vyx,Bz_vyx,fo)
    '''
    #--
    v_grid = grid['v_grid']

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    BierE_x = (1.0 / Eta) * (I1_x + (omega_yx**2) * (I4 / I2) * I3_x)
    BierE_y = (1.0 / Eta) * (I1_y + (omega_yx**2) * (I4 / I2) * I3_y)
    #---- now take the curl?
    return BierE_x, BierE_y


def get_pressureEfield_classic(grid, rho_vyx, Bz_vyx, ne, Te_vyx):
    '''
        ----- should get the kinetic biermann term ----->
       BierE_x,BierE_y =  get_pressureEfield_classic(grid,rho_vyx,Bz_vyx,ne,Te_vyx)
    '''

    # --- construct fo maxwellian
    v_grid = grid['v_grid']
    nv, ny, nx = np.shape(Bz_vyx)
    v = extend_grid_v_to_vxy(nv, ny, nx, v_grid)
    #Te_vyx = extend_grid_xy_to_vxy(nv,ny,nx,Te)
    vte = (2.0 * Te_vyx)**0.5
    fo_unnorm = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    sum_array = np.sum(fo_unnorm, axis=0)
    sum_vyx = np.zeros((np.shape(Te_vyx)))
    for iv in range(nv):
        sum_vyx[iv, :, :] = sum_array
    fo = ne * ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2) / sum_vyx
    # note IMPACT normalisation of fo different to Epperleins

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    BierE_x = (1.0 / Eta) * (I1_x + (omega_yx**2) * (I4 / I2) * I3_x)
    BierE_y = (1.0 / Eta) * (I1_y + (omega_yx**2) * (I4 / I2) * I3_y)
    return BierE_x, BierE_y


def get_Biermann(grid, rho_vyx, Bz_vyx, fo):
    '''
        ----- should get the kinetic biermann term ----->
       kinetic_biermann =  get_Biermann(grid,rho_vyx,Bz_vyx,fo)
    '''
    #--

    v_grid = grid['v_grid']

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    BierE_x = (1.0 / Eta) * (I1_x + (omega_yx**2) * (I4 / I2) * I3_x)
    BierE_y = (1.0 / Eta) * (I1_y + (omega_yx**2) * (I4 / I2) * I3_y)
    #---- now take the curl?
    # --- get_grad assumes [x,y] indices

    x_grid_yx, xY = np.meshgrid(grid['x_grid'], np.ones(len(grid['y_grid'])))
    yY, y_grid_yx = np.meshgrid(np.ones(len(grid['x_grid'])), grid['y_grid'])

    dxBierE_y = (BierE_y[:, 2:] - BierE_y[:, :-2]) / (x_grid_yx[:, 2:] - x_grid_yx[:, :-2])
    dyBierE_x = (BierE_x[2:, :] - BierE_x[:-2, :]) / (y_grid_yx[2:, :] - y_grid_yx[:-2, :])

    kinetic_biermann_z = -1.0 * (dxBierE_y[1:-1, :] - dyBierE_x[:, 1:-1])

    return kinetic_biermann_z


def get_Biermann_classic(grid, rho_vyx, Bz_vyx, ne, Te_vyx):
    '''
        ----- should get the kinetic biermann term ----->
       kinetic_biermann =  get_Biermann(grid,rho_vyx,Bz_vyx,fo)
    '''

    # --- construct fo maxwellian
    v_grid = grid['v_grid']
    nv, ny, nx = np.shape(Bz_vyx)
    v = extend_grid_v_to_vxy(nv, ny, nx, v_grid)
    vte = (2.0 * Te_vyx)**0.5
    fo_unnorm = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    sum_array = np.sum(fo_unnorm, axis=0)
    sum_vyx = np.zeros((np.shape(Te_vyx)))
    for iv in range(nv):
        sum_vyx[iv, :, :] = sum_array
    fo = ne * ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2) / sum_vyx
    # note IMPACT normalisation of fo different to Epperleins

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    BierE_x = (1.0 / Eta) * (I1_x + (omega_yx**2) * (I4 / I2) * I3_x)
    BierE_y = (1.0 / Eta) * (I1_y + (omega_yx**2) * (I4 / I2) * I3_y)
    #---- now take the curl?
    # --- get_grad assumes [x,y] indices

    x_grid_yx, xY = np.meshgrid(grid['x_grid'], np.ones(len(grid['y_grid'])))
    yY, y_grid_yx = np.meshgrid(np.ones(len(grid['x_grid'])), grid['y_grid'])

    dxBierE_y = (BierE_y[:, 2:] - BierE_y[:, :-2]) / (x_grid_yx[:, 2:] - x_grid_yx[:, :-2])
    dyBierE_x = (BierE_x[2:, :] - BierE_x[:-2, :]) / (y_grid_yx[2:, :] - y_grid_yx[:-2, :])

    kinetic_biermann_z = -1.0 * (dxBierE_y[1:-1, :] - dyBierE_x[:, 1:-1])
    return kinetic_biermann_z


def get_Biermann_gradncrossgradT(grid, ne, Te):
    '''
        ----- should get the kinetic biermann term ----->
       kinetic_biermann =  get_Biermann(grid,rho_vyx,Bz_vyx,fo)
    '''
    x_grid, y_grid = grid['x_grid'], grid['y_grid']

    dxTe, dyTe = cf.get_grad(x_grid, y_grid, Te)
    dxne, dyne = cf.get_grad(x_grid, y_grid, ne)

    print 'shapes biermann cross grad = ', np.shape(dxne), np.shape(dyne), np.shape(ne)
    #sys.exit()

    biermann = (1.0 / ne[1:-1, 1:-1]) * (dyTe * dyne - dyTe * dxne)
    return biermann


def data_var(data, var):
    dict = {}
    dict['data'] = data
    dict['name'] = var
    return dict


def get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo):
    '''
        dict = get_qSH_kinetic(grid,rho_vyx,Bz_vyx,jx,jy,fo)
    '''
    v_grid = grid['v_grid']

    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)

    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)
    #========
    K2 = get_K2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    K4 = get_K4_scalar(rho_vyx, omega_vyx, fo, v_grid)
    K1_x, K1_y = get_K1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    K3_x, K3_y = get_K3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

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


def get_alpha_perp_kinetic(grid, rho_vyx, Bz_vyx, fo):
    '''
        alpha_perp_kin = get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo)
    '''
    v_grid = grid['v_grid']
    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)

    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)

    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    alpha_perp_kin = -(1.0 / Eta)

    return alpha_perp_kin


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

    alpha = (1.5) * ((vte * vte) / (ne * rZZni)) * ((v_mom_5 * delta)**-1)

    return alpha


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


def get_v_N_classical(v_grid, Z2ni, ne, Te, w, dxT, dyT, jx=0.0, jy=0.0):

    #Te = 0.5
    vte = (2.0 * Te)**0.5
    rho = Z2ni**-1

    # ------

    F0 = maxw_dist(v_grid, vte)

    beta_wedge_c = get_beta_wedge(w, rho, ne, Te, F0, v_grid)
    beta_wedge = beta_wedge_c

    v_N_x = -beta_wedge * dxT / w
    v_N_y = -beta_wedge * dyT / w

    return v_N_x, v_N_y


def get_vN_from_path(path, fprefix, time):
    '''
        v_nx,v_ny = get_vN_from_path(path,fprefix,time)
    '''

    dict = cf.load_data_all_1D(path, fpre(path), time)
    ne = dict['ne']
    Te = dict['Te']
    qx = dict['qx']
    qy = dict['qy']
    jx = dict['jx']
    jy = dict['jy']
    dxT = dict['dxT']
    dyT = dict['dyT']
    Bz = dict['Bz']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)
    fo = dict_fo['mat']
    grid = fo
    v_grid = dict_fo['v_grid']
    nv, ny, nx = np.shape(fo)

    v_nx = np.zeros((ny, nx))
    v_ny = np.zeros((ny, nx))
    v_nx_hf = np.zeros((ny, nx))
    v_ny_hf = np.zeros((ny, nx))

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)

    for ix in range(nx):
        for iy in range(ny):
            Z2ni_loc = Z2ni[iy, ix]
            v_nx[iy, ix], v_ny[iy, ix] = get_v_N_classical(v_grid, Z2ni_loc, ne[iy, ix], Te[iy, ix],
                                                           Bz[iy, ix], dxT[iy, ix], dyT[iy, ix])

            v_nx_hf[iy, ix] = qx[iy, ix] / (2.5 * ne[iy, ix] * Te[iy, ix])
            v_ny_hf[iy, ix] = qy[iy, ix] / (2.5 * ne[iy, ix] * Te[iy, ix])

    v_nx = np.where(Bz == 0.0, 0.0, v_nx)
    v_ny = np.where(Bz == 0.0, 0.0, v_ny)
    return v_nx, v_ny, v_nx_hf, v_ny_hf


def vmom_fo(m, v_grid, fo):
    '''
     mom = vmom_fo(m,v_grid,fo)
     fo can be 1D array or 3D [v,y,x]
    '''
    dvc = calc_dv(v_grid)    #get_dvc(v_grid)
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


def vmom_f1(m, v_grid, f1x_c, f1y_c):
    '''
    
    '''
    #dvc = get_dvc(v_grid)
    dvc = calc_dv(v_grid)
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

    dict = cf.load_data_all_1D(path, fpre(path), time)
    ne = dict['ne']
    Te = dict['Te']
    qx = dict['qx']
    qy = dict['qy']
    jx = dict['jx']
    jy = dict['jy']
    dxT = dict['dxT']
    dyT = dict['dyT']
    Bz = dict['Bz']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)
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


def get_avgx(mat, ax=1):
    '''
        ax  = [n,..,2,1,0]
        ax = 1 x axis
        ax = 0 y axis
    '''

    avg = np.average(mat, axis=ax)
    return avg


def get_Nernst_abs(path, time):
    '''
    data_x,data_y = get_Nernst_ratio(path,time)
    '''
    dict = cf.load_data_all_1D(path, fpre(path), time)
    ne = dict['ne']
    Te = dict['Te']
    qx = dict['qx']
    qy = dict['qy']
    jx = dict['jx']
    jy = dict['jy']
    dxT = dict['dxT']
    dyT = dict['dyT']
    Bz = dict['Bz']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)
    fo = dict_fo['mat']

    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #------ NERNST -----------------------------------------------------
    v_nx, v_ny = get_v_N(grid, rho_vyx, Bz_vyx, fo)
    v_nx_classical, v_ny_classical, v_nx_hf, v_ny_hf = get_vN_from_path(path, fpre(path), time)

    v_nx_c = np.transpose(v_nx_classical)
    v_ny_c = np.transpose(v_ny_classical)    #np.transpoe(v_ny_classical)#[:,::-1]
    v_nx_k = np.transpose(v_nx)
    v_ny_k = np.transpose(v_ny)    #np.transpose(v_ny)#[:,::-1]

    return v_nx_k, v_ny_k, v_nx_c, v_ny_c


def get_q_abs(path, time):
    '''
    dict_ratio = get_q_ratio(path,time)
    '''
    time_int = int(time)
    dict = cf.load_data_all_1D(path, fpre(path), time)
    ne = dict['ne']
    Te = dict['Te']
    qx = dict['qx']
    qy = dict['qy']
    jx = dict['jx']
    jy = dict['jy']
    dxT = dict['dxT']
    dyT = dict['dyT']
    Bz = dict['Bz']
    U = dict['U']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)
    fo = dict_fo['mat']
    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)
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
        dict[name]['data'] = np.transpose(dict[name]['data'])[:]

    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
    #--- classical

    q_dir = path + '/'
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

    dict_c = {}
    dict_c['tot x'] = q_x_B[time_int, :, :]    #q_xc
    dict_c['tot y'] = q_y_B[time_int, :, :]

    dict_c['SH x'] = q_SH_x[time_int, :, :]    #q_SH_x[time_int,:x_lim,:]
    dict_c['SH y'] = q_SH_y[time_int, :, :]    #q_SH_y[time_int,:x_lim,:]
    dict_c['RL x'] = q_RL_x[time_int, :, :]
    dict_c['RL y'] = q_RL_y[time_int, :, :]
    dict_c['E x'] = q_E_x[time_int, :, :]
    dict_c['E y'] = q_E_y[time_int, :, :]

    dict_k = {}
    dict_k['tot x'] = q_x_VFP[time_int, :, :]
    dict_k['tot y'] = q_y_VFP[time_int, :, :]

    dict_k['SH x'] = q_KSH_x    #q_SH_x[time_int,:x_lim,:]
    dict_k['SH y'] = q_KSH_y    #q_SH_y[time_int,:x_lim,:]
    dict_k['RL x'] = q_KRL_x
    dict_k['RL y'] = q_KRL_y
    dict_k['E x'] = q_KTE_x + q_KE_x
    dict_k['E y'] = q_KTE_y + q_KE_y

    dict_k['time'] = dict_fo['time']
    dict_c['x_grid'] = x_grid[1:-1]
    dict_c['y_grid'] = y_grid
    dict_k['x_grid'] = x_grid[1:-1]
    dict_k['y_grid'] = y_grid
    dict_k['ne'] = np.transpose(ne)
    dict_k['Te'] = np.transpose(Te)
    dict_k['jx'] = np.transpose(jx)
    dict_k['jy'] = np.transpose(jy)
    dict_k['dxT'] = np.transpose(dxT)
    dict_k['dyT'] = np.transpose(dyT)
    dict_k['Bz'] = np.transpose(Bz)
    dict_k['Z2ni'] = np.transpose(Z2ni)
    dict_k['U'] = np.transpose(U)

    return dict_c, dict_k


def get_kinetic_and_classicB(path, fprefix, time, xlim=-1):
    '''
       classic_biermann,kinetic_biermann = get_kinetic_and_classicB(path,fprefix,time)
    '''
    dict = cf.load_data_all_1D(path, fpre(path), time)
    ne = dict['ne']
    Te = dict['Te']
    qx = dict['qx']
    qy = dict['qy']
    jx = dict['jx']
    jy = dict['jy']
    dxT = dict['dxT']
    dyT = dict['dyT']
    Bz = dict['Bz']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)
    fo = dict_fo['mat']
    time_col = dict_fo['time'] * tstep_factor

    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)

    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    Te_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Te)
    ne_vyx = extend_grid_xy_to_vxy(nv, ny, nx, ne)

    #------ NERNST -----------------------------------------------------
    print '\n\n--------- time ==================== ', time

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    print(' bier init= ', np.shape(rho_vyx), np.shape(Bz_vyx), np.shape(ne_vyx), np.shape(Te_vyx))

    classic_biermann = get_Biermann_classic(grid, rho_vyx, Bz_vyx, ne_vyx, Te_vyx)
    kinetic_biermann = get_Biermann(grid, rho_vyx, Bz_vyx, fo)

    classic_biermann = np.transpose(classic_biermann)
    kinetic_biermann = np.transpose(kinetic_biermann)

    return classic_biermann, kinetic_biermann
