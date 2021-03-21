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
#-----------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import sys, re, os, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import chfoil_module as cf
#import q_SH_Te_tsteps as q_mod
import matplotlib as mpl
from pylab import *

nv = 100
Z = 6.51
#ones_nv = np.ones((nv))
#vmax = 10.0
#xmin,xmax,nx = 0.0,142.662,40
#ymin,ymax,ny = 0.0,35.0,10
#vmin,vmax,nv = 0.0,10.0,100

#vb_grid = np.arange(vmin,vmax,(vmax-vmin)/(nv+1))
#v_grid = 0.5*(vb_grid[1:] + vb_grid[:-1])
#dv = v_grid[1] - v_grid[0]
#v2dv = (v_grid**2)*dv

#v_grid = np.arange(vmin,vmax,(vmax-vmin)/nv)
#---- get vel grid---
#ic = np.arange(nv)
#ib = np.arange(-0.5,nv+0.5)
#vb_grid = cf.interp_data(ic,v_grid,ib)
#dvc = vb_grid[1:] - vb_grid[:-1]
Y1_weight = 4.0 * np.pi / 3.0
Y0_weight = 4.0 * np.pi

#-------------------------------------->
# get normalisations
norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + 'norm.log'
[T0, n0, Z0, Bz0] = np.loadtxt(log_file)
'''
cd5 = cf.conv_factors_custom(norm_path,Z0,Ar=6.51)

print 'Z0 = ', Z0
print 'T0 = ', T0
print 'n0 = ', n0
cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index,c_index = 0,-1
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab_rel
ylab = cd5.ylab
'''
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


#-----------------------------------------------------------------------
def fpre(path):
    return path.split('/')[-1]


#-----------------------------------------------------------------------


def vc_to_vb(vc):
    vb = np.zeros((len(vc) + 1))
    vred = 0.5 * (vc[1:] + vc[:-1])
    vb[1:-1] = vred
    dv0 = vred[1] - vred[0]
    dv1 = vred[-1] - vred[-2]
    vb[0] = vred[0] - dv0
    vb[-1] = vred[-1] + dv1
    return vb


#-----------------------------------------------------------------------


def calc_dv(vc):
    '''
        dv = calc_dv(vc)
    '''
    vb = vc_to_vb(vc)
    return vb[1:] - vb[:-1]


#-----------------------------------------------------------------------
def extend_grid_xy_to_vxy(nv, ny, nx, grid_xy):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iv in range(nv):
        grid_vyx[iv, :, :] = grid_xy
    return grid_vyx


#-----------------------------------------------------------------------


def extend_grid_v_to_vxy(nv, ny, nx, grid_v):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iv in range(nv):
        grid_vyx[iv, :, :] = grid_v[iv]
    return grid_vyx


#-----------------------------------------------------------------------


def extend_grid_y_to_vxy(nv, ny, nx, grid_y):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for iy in range(ny):
        grid_vyx[:, iy, :] = grid_y[iy]
    return grid_vyx


#-----------------------------------------------------------------------


def extend_grid_x_to_vxy(nv, ny, nx, grid_x):
    '''
        
    '''
    grid_vyx = np.zeros((nv, ny, nx))
    for ix in range(nx):
        grid_vyx[:, :, ix] = grid_x[ix]
    return grid_vyx


#-----------------------------------------------------------------------
''' 
OBSOLETE (therefore removed) - becuause this does not work with a dynamically varying Z -
def get_omega(v_grid,ni,Bz):
    nu_ei = (Z**2)*ni*(v_grid**-3)
    omega = Bz*(nu_ei**-1)
    return omega
    
'''
#-----------------------------------------------------------------------


def get_v2dv(v_grid):
    #dv = v_grid[1:] - v_grid[:-1]
    dv = calc_dv(v_grid)
    v2dv = (v_grid**2) * dv
    return v2dv


#-----------------------------------------------------------------------


def get_v2dv_vyx(v_grid):
    dv = v_grid[1] - v_grid[0]
    v2dv = (v_grid**2) * dv
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, grid_v)
    return v2dv_vyx


#-----------------------------------------------------------------------


def maxw_dist(v, vte):
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    return f_om


#-----------------------------------------------------------------------
def maxw_dist_impact(v, ne, vte):
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-(v / vte)**2)
    return f_om


#-----------------------------------------------------------------------


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
    ##print 'DELTA = ', delta
    return delta


#-----------------------------------------------------------------------
def get_fo_mom_int(m, rho_vyx, omega_vyx, fo, v_grid):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni
    # omega  = Bz/Z2ni = Bz*rho

    nv, ny, nx = np.shape(fo)
    #print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    v2dv = get_v2dv(v_grid)
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)
    ones_vyx = np.ones((nv, ny, nx))

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / (ones_vyx + (omega_vyx**2) * v_grid_vyx**6)
    mom = np.sum(prefactor * int, axis=0)
    return mom


#-----------------------------------------------------------------------


def get_fo_mom_int2(m, rho_vyx, omega_vyx, fo, v_grid):
    prefactor = (4.0 * np.pi / 3.0) * rho_vyx    # rho = 1/Z2ni

    nv, ny, nx = np.shape(fo)
    ##print ' int 2 = '
    ##print ' nv = %i ny = %i nx = %i' % (nv,ny,nx)
    v2dv = get_v2dv(v_grid)
    ones_vyx = np.ones((nv, ny, nx))
    v2dv_vyx = extend_grid_v_to_vxy(nv, ny, nx, v2dv)
    v_grid_vyx = extend_grid_v_to_vxy(nv, ny, nx, v_grid)

    int = ((v_grid_vyx**m) * v2dv_vyx * fo) / ((ones_vyx + (omega_vyx**2) * (v_grid_vyx**6))**2)
    mom = np.sum(prefactor * int, axis=0)
    return mom


#-----------------------------------------------------------------------


def get_dfodv_int(n, rho_vyx, omega_vyx, fo, v_grid):
    # omega  = Bz/Z2ni = Bz*rho
    #print ' ---- dfodv --- int -----'
    v_mom_nm3 = get_fo_mom_int(n - 3, rho_vyx, omega_vyx, fo, v_grid)
    v_mom2_np3 = get_fo_mom_int2(n + 3, rho_vyx, omega_vyx, fo, v_grid)
    mom_out = -1.0 * (n * v_mom_nm3 - 6.0 * (omega_vyx[0, :, :]**2) * v_mom2_np3
                     )    # unclear whether this minus one should actually be there or not
    return mom_out


def get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- I1 -----'
    m = 5
    mom_x = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' I2_scalar'
    m = 6
    mom = get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I4 = get_I4_scalar(rho,omega,fo)
    '''
    #print ' I4_scalar'

    m = 9
    mom = rho_vyx[0, :, :] * get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I3 =  get_I3_vec(rho,omega,gradfo)
    '''
    #print ' ---- I3 -----'
    m = 8
    mom_x = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- K1 -----'
    m = 7
    mom_x = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = 0.5 * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K2_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
    m = 8
    mom = 0.5 * get_dfodv_int(m, rho_vyx, omega_vyx, fo, v_grid)
    return mom


def get_K3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid):
    '''
        I1 =  get_I1_vec(rho,omega,gradfo)
    '''

    #print ' ---- K3 -----'
    m = 10
    mom_x = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_x, v_grid)
    mom_y = 0.5 * rho_vyx[0, :, :] * get_fo_mom_int(m, rho_vyx, omega_vyx, gradfo_y, v_grid)
    return mom_x, mom_y


def get_K4_scalar(rho_vyx, omega_vyx, fo, v_grid):
    '''
        I2 = get_I2_scalar(rho,omega,fo)
    '''
    # omega  = Bz/Z2ni = Bz*rho
    #print ' K2_scalar'
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
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
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
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
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


#-----------------------------------------------------------------------
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


#-----------------------------------------------------------------------
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


#-----------------------------------------------------------------------
def get_Biermann(grid, rho_vyx, Bz_vyx, fo):
    '''
        ----- should get the kinetic biermann term ----->
       kinetic_biermann =  get_Biermann(grid,rho_vyx,Bz_vyx,fo)
    '''
    #--
    '''
    nv,ny,nx = np.shape(fo)
    
    v2dv = get_v2dv(v_grid)
    mom2 = np.zeros((ny,nx))
    for iy in range(ny):
        for ix in range(nx):
           mom2[iy,ix] =  get_v_mom_m(v_grid,0.0,fo[:,iy,ix],0)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    im = ax.imshow(mom2)
    plt.colorbar(im,ax=ax)
    plt.show()
    sys.exit()    
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
    BierE_x = (1.0 / Eta) * (I1_x + (omega_yx**2) * (I4 / I2) * I3_x)
    BierE_y = (1.0 / Eta) * (I1_y + (omega_yx**2) * (I4 / I2) * I3_y)
    #---- now take the curl?
    # --- get_grad assumes [x,y] indices
    #print ' np.shape(BierE_x) = ', np.shape(BierE_x)

    x_grid_yx, xY = np.meshgrid(grid['x_grid'], np.ones(len(grid['y_grid'])))
    yY, y_grid_yx = np.meshgrid(np.ones(len(grid['x_grid'])), grid['y_grid'])

    #print '\n\n x_grid = ', x_grid_yx, '\n xY = ',xY, np.shape(x_grid_yx)

    #print ' y_grid = ', y_grid_yx, '\n yY = ', yY, np.shape(y_grid_yx)

    #dyBierE_x, dxBierE_x = cf.get_grad(grid['y_grid'],grid['x_grid'],BierE_x)
    #dyBierE_y, dxBierE_y = cf.get_grad(grid['y_grid'],grid['x_grid'],BierE_y)
    #print ' np.shape(BierE_x) = ',np.shape(I1_x), np.shape(omega_yx), np.shape(grid['y_grid']), np.shape(grid['x_grid']), np.shape(BierE_x),np.shape(BierE_y)

    dxBierE_y = (BierE_y[:, 2:] - BierE_y[:, :-2]) / (x_grid_yx[:, 2:] - x_grid_yx[:, :-2])
    dyBierE_x = (BierE_x[2:, :] - BierE_x[:-2, :]) / (y_grid_yx[2:, :] - y_grid_yx[:-2, :])

    #print ' np.shape(BierE_x) = ', np.shape(BierE_x), np.shape(I1_x), np.shape(omega_yx), np.shape(grid['y_grid']), np.shape(grid['x_grid']), np.shape(dyBierE_x)

    #print ' np.shape(dxBierE_y) = ', np.shape(dxBierE_y), np.shape(dyBierE_x)
    kinetic_biermann_z = -1.0 * (dxBierE_y[1:-1, :] - dyBierE_x[:, 1:-1])
    #print ' shape kinetic_biermann_z = ', np.shape(kinetic_biermann_z)
    #sys.exit()
    return kinetic_biermann_z


#-----------------------------------------------------------------------
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
    print ' np.shape(BierE_x) = ', np.shape(BierE_x), 'len ygrid = ', np.shape(
        grid['y_grid']), 'len x = ', np.shape(grid['x_grid'])

    x_grid_yx, xY = np.meshgrid(grid['x_grid'], np.ones(len(grid['y_grid'])))
    yY, y_grid_yx = np.meshgrid(np.ones(len(grid['x_grid'])), grid['y_grid'])

    print '\n\n x_grid = ', x_grid_yx, '\n xY = ', xY, np.shape(x_grid_yx)
    print ' y_grid = ', y_grid_yx, '\n yY = ', yY, np.shape(y_grid_yx)

    #dyBierE_x, dxBierE_x = cf.get_grad(grid['y_grid'],grid['x_grid'],BierE_x)
    #dyBierE_y, dxBierE_y = cf.get_grad(grid['y_grid'],grid['x_grid'],BierE_y)
    print ' np.shape(BierE_x) = ', np.shape(I1_x), np.shape(omega_yx), np.shape(
        grid['y_grid']), np.shape(grid['x_grid']), np.shape(BierE_x), np.shape(BierE_y), np.shape(
            x_grid_yx), np.shape(y_grid_yx)

    dxBierE_y = (BierE_y[:, 2:] - BierE_y[:, :-2]) / (x_grid_yx[:, 2:] - x_grid_yx[:, :-2])
    dyBierE_x = (BierE_x[2:, :] - BierE_x[:-2, :]) / (y_grid_yx[2:, :] - y_grid_yx[:-2, :])

    kinetic_biermann_z = -1.0 * (dxBierE_y[1:-1, :] - dyBierE_x[:, 1:-1])
    return kinetic_biermann_z


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
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


#-----------------------------------------------------------------------


def data_var(data, var):
    dict = {}
    dict['data'] = data
    dict['name'] = var
    return dict


#-----------------------------------------------------------------------


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


#-----------------------------------------------------------------------


def get_alpha_perp_kinetic(grid, rho_vyx, Bz_vyx, fo):
    '''
        alpha_perp_kin = get_alpha_perp_kinetic(grid,rho_vyx,Bz_vyx,fo)
    '''
    v_grid = grid['v_grid']
    omega_vyx = Bz_vyx * rho_vyx
    gradfo_x, gradfo_y = cf.get_grad_3d_varZ(grid, fo)
    ##print ' np.shape(grad_fo) = ',np.shape(gradfo_x),np.shape(gradfo_y)
    I1_x, I1_y = get_I1_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    I3_x, I3_y = get_I3_vec(rho_vyx, omega_vyx, gradfo_x, gradfo_y, v_grid)
    ##print ' got I1_X = ',I1_x, I1_y,I3_x,I3_y
    I2 = get_I2_scalar(rho_vyx, omega_vyx, fo, v_grid)
    I4 = get_I4_scalar(rho_vyx, omega_vyx, fo, v_grid)

    omega_yx = omega_vyx[0, :, :]
    Bz_yx = Bz_vyx[0, :, :]
    Eta = (I2 + (omega_yx**2) * ((I4**2) / I2))
    alpha_perp_kin = -(1.0 / Eta)

    return alpha_perp_kin


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


def get_v_N_classical(v_grid, Z2ni, ne, Te, w, dxT, dyT, jx=0.0, jy=0.0):

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
    print(' shape before = ', np.shape(ne), np.shape(prof_Z), np.shape(ne), nx, ny)
    print(' shape before = ', np.shape(ne), np.shape(prof_Z), nx, ny)
    #Z2ni = ne*prof_Z#np.transpose(cf.trim_array(ne,nx,ny))*np.transpose(cf.trim_array(prof_Z,nx,ny))
    #Bz = np.transpose(cf.trim_array(Bz,nx,ny))
    #jx = np.transpose(cf.trim_array(jx_c,nx,ny))
    #jy = np.transpose(cf.trim_array(jy_c,nx,ny))
    ''''
    #-------
    dxT,dyT = cf.get_grad(x_grid,y_grid,Te)
    Te = np.transpose(cf.trim_array(Te,nx,ny))
    ne = np.transpose(cf.trim_array(ne,nx,ny))
    Z2ni = np.transpose(cf.trim_array(Z2ni,nx,ny))
    dxT = np.transpose(cf.trim_array(dxT,nx,ny))
    dyT = np.transpose(cf.trim_array(dyT,nx,ny))*-1.0
    qx_c = np.transpose(qx_c)
    qy_c = np.transpose(qy_c)
    
    '''
    #ny,nx = np.shape(qx_c)#len(T_data[:,0])-1,len(T_data[0,:])-1
    v_nx = np.zeros((ny, nx))
    v_ny = np.zeros((ny, nx))
    v_nx_hf = np.zeros((ny, nx))
    v_ny_hf = np.zeros((ny, nx))

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    print('shape after = ', np.shape(Z2ni), np.shape(Te), np.shape(dxT), np.shape(dyT), ny, nx)
    for ix in range(nx):
        for iy in range(ny):
            Z2ni_loc = Z2ni[iy, ix]
            v_nx[iy, ix], v_ny[iy, ix] = get_v_N_classical(v_grid, Z2ni_loc, ne[iy, ix], Te[iy, ix],
                                                           Bz[iy, ix], dxT[iy, ix], dyT[iy, ix])

            v_nx_hf[iy, ix] = qx[iy, ix] / (2.5 * ne[iy, ix] * Te[iy, ix])
            v_ny_hf[iy, ix] = qy[iy, ix] / (2.5 * ne[iy, ix] * Te[iy, ix])
    ##print 'shapes = ',np.shape(Bz),np.shape(v_nx)
    v_nx = np.where(Bz == 0.0, 0.0, v_nx)
    v_ny = np.where(Bz == 0.0, 0.0, v_ny)
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


'''
def get_gradT(x_grid,y_grid,T_data):

    
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
'''

#-----------------------------------------------------------------------


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


#-----------------------------------------------------------------------


def convert_var_to_str(var):
    for name in globals():
        if eval(name) == var:
            return name


#-----------------------------------------------------------------------
'''
def plot_2D(ax,data,label,middle=None,colormap='jet',#limits=None,xlabel=None,ylabel=None,clabel=None):
    if middle != None:
        try:
            norm = cf.MidPointNorm(middle)
        except ValueError:
            norm = None
    else:
        norm = None
    
    if #limits:
        try:
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm,extent=#limits)
            plt.colorbar(im,ax=ax,aspect='auto',norm=norm,label=clabel)
        except ValueError:
            norm=None
            #print 'setting norm to None;'
            im = ax.imshow(data,aspect='auto',cmap=colormap,norm=norm,extent=#limits)
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
'''


#-----------------------------------------------------------------------
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


#-----------------------------------------------------------------------


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


#-----------------------------------------------------------------------

#-----------------------------------------------------------------------


def get_avgx(mat, ax=1):
    '''
        ax  = [n,..,2,1,0]
        ax = 1 x axis
        ax = 0 y axis
    '''

    avg = np.average(mat, axis=ax)
    return avg


#-----------------------------------------------------------------------
def get_Nernst_ratio(path, time):
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
    v_ny_kinetic = np.transpose(v_ny)[:x_lim, :]

    data_y = v_ny_kinetic / v_ny_classical[:x_lim, :]
    data_y = np.where(v_ny_classical[:x_lim, :] == 0.0, 1.0, data_y)
    vmin, vmax = -1.5, 1.5
    data_y = np.where(np.abs(data_y) > vmax, np.sign(data_y) * vmax, data_y)

    data_x = v_nx_kinetic / v_nx_classical[:x_lim, :]
    data_x = np.where(v_nx_classical[:x_lim, :] == 0.0, 1.0, data_x)
    #vmin,vmax = -1.0,1.0
    data_x = np.where(np.abs(data_x) > vmax, np.sign(data_x) * vmax, data_x)
    return data_x, data_y


#-----------------------------------------------------------------------
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
    '''
    data_y = v_ny_kinetic/v_ny_classical[:x_lim,:]
    data_y = np.where(v_ny_classical[:x_lim,:]==0.0,1.0,data_y)
    vmin,vmax = -1.5,1.5
    data_y = np.where(np.abs(data_y)>vmax,np.sign(data_y)*vmax,data_y)

    data_x = v_nx_kinetic/v_nx_classical[:x_lim,:]
    data_x = np.where(v_nx_classical[:x_lim,:]==0.0,1.0,data_x)
    #vmin,vmax = -1.0,1.0
    data_x = np.where(np.abs(data_x)>vmax,np.sign(data_x)*vmax,data_x)
    '''

    return v_nx_k, v_ny_k, v_nx_c, v_ny_c


#-----------------------------------------------------------------------
def get_ratio_lim(q_yc, q_yk, vmax=1.5, thresh=1e-2):
    data_y = q_yk / q_yc
    #data_y = np.where(q_yc==0.0,1.0,data_y)
    #thresh=1e-2
    #data_y = np.where(np.abs(q_yc)<=thresh,1.0,data_y)
    #data_y = np.where(np.abs(q_yk)<=thresh,1.0,data_y)

    data_y = np.where(np.abs(data_y) > vmax, np.sign(data_y) * vmax, data_y)

    return data_y


#-----------------------------------------------------------------------


def get_q_ratio(path, time):
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

    #------ total
    q_xk = q_x_VFP[time_int, :x_lim, :]
    q_xc = q_x_B[time_int, :x_lim, :]
    avg = np.average(np.abs(q_xc))
    rat_tot_x = get_ratio_lim(q_xc, q_xk)
    rat_tot_x = np.where(np.abs(q_xc) <= avg * 1e-4, 0.0, rat_tot_x)

    q_yk = q_y_VFP[time_int, :x_lim, :]
    q_yc = q_y_B[time_int, :x_lim, :]
    rat_tot_y = get_ratio_lim(q_yc, q_yk)
    avg = np.average(np.abs(q_yc))
    rat_tot_y = np.where(np.abs(q_yc) <= avg * 1e-2, 0.0, rat_tot_y)

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

    #------

    dict_ratio = {}
    dict_ratio['tot x'] = rat_tot_x
    dict_ratio['tot y'] = rat_tot_y

    dict_ratio['SH x'] = rat_SH_x    #q_SH_x[time_int,:x_lim,:]
    dict_ratio['SH y'] = rat_SH_y    #q_SH_y[time_int,:x_lim,:]
    dict_ratio['RL x'] = rat_RL_x
    dict_ratio['RL y'] = rat_RL_y

    return dict_ratio


#-----------------------------------------------------------------------


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


#-----------------------------------------------------------------------
def get_q_abs_c(path, time):
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
    dxB = dict['dxB']
    U = dict['U']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    Z2ni = dict['Z2ni']
    prof_Z = dict['Z']
    fo = dict['fo']
    nv = dict['nv']
    ny = dict['ny']
    nx = dict['nx']
    wt = dict['wt']

    rA = (Z2ni)**-1

    dict_fo = cf.load_dict_1D(path, fpre(path), 'fo', time)
    fo = dict_fo['mat']
    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d(grid, fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #--- constructing fo maxw
    ##########
    # extra extensions requried for maxwellian f0
    Te_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Te)
    ne_vyx = extend_grid_xy_to_vxy(nv, ny, nx, ne)

    #heat flow ratios

    #--- kinetic
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)
    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']
    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])[:]

    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']

    #--------------------------------------------------------------------------
    fo_m = get_fo_maxw_3D(dict_fo['v_grid'], ne_vyx, Te_vyx)
    dict_c_in = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo_m)
    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']
    for name in name_list:
        dict_c_in[name]['data'] = np.transpose(dict_c_in[name]['data'])[:]

    q_SH_x, q_SH_y = dict_c_in['q_SH_x']['data'], dict_c_in['q_SH_y']['data']
    q_RL_x, q_RL_y = dict_c_in['q_RL_x']['data'], dict_c_in['q_RL_y']['data']
    q_TE_x, q_TE_y = dict_c_in['q_TE_x']['data'], dict_c_in['q_TE_y']['data']
    q_E_x, q_E_y = dict_c_in['q_E_x']['data'], dict_c_in['q_E_y']['data']
    #--- classical

    q_dir = path + '/'

    q_xX = np.load(q_dir + 'q_xX.txt.npy')
    q_yY = np.load(q_dir + 'q_yY.txt.npy')
    #nt,ny,nx
    q_x_VFP = (q_xX[:, 1:, :] + q_xX[:, :-1, :]) * 0.5
    q_y_VFP = (q_yY[:, :, 1:] + q_yY[:, :, :-1]) * 0.5
    q_x_B = q_RL_x + q_E_x + q_SH_x + q_TE_x
    q_y_B = q_RL_y + q_E_y + q_SH_y + q_TE_y

    #------

    dict_c = {}
    dict_c['tot x'] = q_x_B[:, :]    #q_xc
    dict_c['tot y'] = q_y_B[:, :]

    dict_c['SH x'] = q_SH_x[:, :]    #q_SH_x[time_int,:x_lim,:]
    dict_c['SH y'] = q_SH_y[:, :]    #q_SH_y[time_int,:x_lim,:]
    dict_c['RL x'] = q_RL_x[:, :]
    dict_c['RL y'] = q_RL_y[:, :]
    dict_c['E x'] = q_TE_x + q_E_x[:, :]
    dict_c['E y'] = q_TE_y + q_E_y[:, :]

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
    dict_k['dxB'] = np.transpose(dxB)

    dict_k['dyT'] = np.transpose(dyT)
    dict_k['Bz'] = np.transpose(Bz)
    dict_k['Z2ni'] = np.transpose(Z2ni)
    dict_k['U'] = np.transpose(U)
    dict_k['wt'] = np.transpose(wt)

    return dict_c, dict_k


#------------------------------------------------------------------------------

#-----------------------------------------------------------------------


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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
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

    #----- image plot setup
    #limits = [#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[xlim],#x_grid_SI[0]]
    #lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    #lab_dict['lims'] = #limits#[#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[c_index],#x_grid_SI[0]]
    #lab_dict['colormap'] = 'jet'
    #lab_dict['xlab'],#lab_dict['ylab'] = xlab,ylab
    #lab_dict['cbar_title'] = '   '
    #lab_dict['title'] = '   '
    #lab_dict['middle'] = 0.0

    #------ NERNST -----------------------------------------------------
    print '\n\n--------- time ==================== ', time

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    print(' bier init= ', np.shape(rho_vyx), np.shape(Bz_vyx), np.shape(ne_vyx), np.shape(Te_vyx))

    classic_biermann = get_Biermann_classic(grid, rho_vyx, Bz_vyx, ne_vyx, Te_vyx)
    kinetic_biermann = get_Biermann(grid, rho_vyx, Bz_vyx, fo)
    #print(' bier out = ', np.shape(classic_biermann),np.shape(kinetic_biermann), classic_biermann, kinetic_biermann)

    classic_biermann = np.transpose(classic_biermann)
    kinetic_biermann = np.transpose(kinetic_biermann)

    return classic_biermann, kinetic_biermann


#-----------------------------------------------------------------------


def get_kinetic_and_classic_E_pressure(path, fprefix, time, xlim=73):
    '''
       dict_classical,dict_kinetic =get_kinetic_and_classic_E_pressure(path,fprefix,time,xlim=73)
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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
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

    #----- image plot setup
    #limits = [#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[xlim],#x_grid_SI[0]]
    #lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    #lab_dict['lims'] = #limits#[#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[c_index],#x_grid_SI[0]]
    #lab_dict['colormap'] = 'jet'
    #lab_dict['xlab'],#lab_dict['ylab'] = xlab,ylab
    #lab_dict['cbar_title'] = '   '
    #lab_dict['title'] = '   '
    #lab_dict['middle'] = 0.0

    #------ NERNST -----------------------------------------------------
    print '\n\n--------- time ==================== ', time

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)

    E_Pressure_c_x, E_Pressure_c_y = get_pressureEfield_classic(grid, rho_vyx, Bz_vyx, ne, Te_vyx)
    E_Pressure_k_x, E_Pressure_k_y = get_pressureEfield_kinetic(grid, rho_vyx, Bz_vyx, fo)

    dict_classic = {}
    dict_classic['E_P_x'] = E_Pressure_c_x
    dict_classic['E_P_y'] = E_Pressure_c_y

    dict_kinetic = {}
    dict_kinetic['E_P_x'] = E_Pressure_k_x
    dict_kinetic['E_P_y'] = E_Pressure_k_y

    return dict_classical, dict_kinetic


def get_kinetic_and_classic_E_pressure(path, fprefix, time, xlim=73):
    '''
       dict_classical,dict_kinetic =get_kinetic_and_classic_E_pressure(path,fprefix,time,xlim=73)
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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
    time_col = dict_fo['time'] * tstep_factor

    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)

    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    jx_vyx = extend_grid_xy_to_vxy(nv, ny, nx, jx)
    jy_vyx = extend_grid_xy_to_vxy(nv, ny, nx, jy)

    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    Te_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Te)
    ne_vyx = extend_grid_xy_to_vxy(nv, ny, nx, ne)

    #----- image plot setup
    #limits = [#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[xlim],#x_grid_SI[0]]
    #lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    #lab_dict['lims'] = #limits#[#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[c_index],#x_grid_SI[0]]
    #lab_dict['colormap'] = 'jet'
    #lab_dict['xlab'],#lab_dict['ylab'] = xlab,ylab
    #lab_dict['cbar_title'] = '   '
    #lab_dict['title'] = '   '
    #lab_dict['middle'] = 0.0

    #------ NERNST -----------------------------------------------------
    print '\n\n--------- time ==================== ', time

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)

    E_Pressure_c_x, E_Pressure_c_y = get_pressureEfield_classic(grid, rho_vyx, Bz_vyx, ne, Te_vyx)
    E_Pressure_k_x, E_Pressure_k_y = get_pressureEfield_kinetic(grid, rho_vyx, Bz_vyx, fo)
    E_hall_x, E_hall_y = get_hallfield(grid, rho_vyx, Bz_vyx, jx_vyx, jy_vyx, fo)
    E_TEwedgex, E_TEwedgey = get_thermoelectricEfield_kinetic(grid, rho_vyx, Bz_vyx, fo)

    dict_classic = {}
    dict_classic['E_P_x'] = E_Pressure_c_x
    dict_classic['E_P_y'] = E_Pressure_c_y

    dict_kinetic = {}
    dict_kinetic['E_P_x'] = E_Pressure_k_x
    dict_kinetic['E_P_y'] = E_Pressure_k_y

    return dict_classical, dict_kinetic


#=======================================================================
def get_alpha_perp_path(path, time):
    '''
     alpha_perp_classical,alpha_perp_kinetic = get_alpha_perp_path(path,time)
    '''
    #--- load in data ....
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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
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

    rA = (Z2ni)**-1

    #===============
    #grid = dict['grid']
    #dxfo, dyfo = cf.get_grad_3d_varZ(grid,fo)
    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    alpha_perp_kinetic = get_alpha_perp_kinetic(grid, rho_vyx, Bz_vyx, fo)

    #---- perform the classical calculation
    T_data = Te    #np.transpose(cf.trim_array(Te,nx,ny))
    n_data = ne    #np.transpose(cf.trim_array(ne,nx,ny))
    Bz_data = Bz
    #---------------------------------------
    alpha_perp_classical = np.zeros((np.shape(T_data)))
    #------
    #print ' shapes = ', np.shape(T_data),np.shape(n_data), np.shape(Bz_data), np.shape(alpha_perp_kinetic)
    #for ix in range(1,len(T_data[0,:])-1):
    #   for iy in range(1,len(T_data[:,0])-1):
    for ix in range(nx):
        for iy in range(ny):
            Z2ni_loc = Z2ni[iy, ix]
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


def get_kinetic_E(path, fprefix, time, xlim=73):
    '''
        dict_kinetic =get_kinetic_E(path,fprefix,time,xlim=73)
    '''
    #--- load in data ....
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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
    time_col = dict_fo['time'] * tstep_factor

    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)

    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
    time_col = dict_te['time'] * tstep_factor
    dict_fo = cf.load_dict_1D(path, fprefix, 'fo', time)

    grid = dict['grid']
    nv, ny, nx = np.shape(fo)
    #===============
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    jx_vyx = extend_grid_xy_to_vxy(nv, ny, nx, jx)
    jy_vyx = extend_grid_xy_to_vxy(nv, ny, nx, jy)
    jx = jx_vyx[0, :, :]
    jy = jy_vyx[0, :, :]

    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    Te_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Te)
    ne_vyx = extend_grid_xy_to_vxy(nv, ny, nx, ne)

    #----- image plot setup
    #limits = [#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[xlim],#x_grid_SI[0]]
    #lab_dict = {}
    #lab_dict['figsize'] = (6,6) # x,y inches
    #lab_dict['lims'] = #limits#[#y_grid_SI[0],#y_grid_SI[-1],#x_grid_SI[c_index],#x_grid_SI[0]]
    #lab_dict['colormap'] = 'jet'
    #lab_dict['xlab'],#lab_dict['ylab'] = xlab,ylab
    #lab_dict['cbar_title'] = '   '
    #lab_dict['title'] = '   '
    #lab_dict['middle'] = 0.0

    #------ NERNST -----------------------------------------------------
    print '\n\n--------- time ==================== ', time

    #v_nx,v_ny = get_v_N(grid,rho_vyx,Bz_vyx,fo)
    alpha_perp_classical, alpha_perp_kinetic = get_alpha_perp_path(path, time)

    E_Pressure_c_x, E_Pressure_c_y = get_pressureEfield_classic(grid, rho_vyx, Bz_vyx, ne, Te_vyx)
    E_Pressure_k_x, E_Pressure_k_y = get_pressureEfield_kinetic(grid, rho_vyx, Bz_vyx, fo)
    E_hall_x, E_hall_y = get_hallfield(grid, rho_vyx, Bz_vyx, jx_vyx, jy_vyx, fo)
    E_TEwedge_x, E_TEwedge_y = get_thermoelectricEfield_kinetic(grid, rho_vyx, Bz_vyx, fo)

    dict_kinetic = {}
    dict_kinetic['alpha_perp'] = alpha_perp_kinetic
    dict_kinetic['alpha_jx'] = alpha_perp_kinetic * jx
    dict_kinetic['alpha_jy'] = alpha_perp_kinetic * jy

    dict_kinetic['E_P_x'] = E_Pressure_k_x
    dict_kinetic['E_hall_x'] = E_hall_x
    dict_kinetic['E_TEwedge_x'] = E_TEwedge_x

    dict_kinetic['E_P_y'] = E_Pressure_k_y
    dict_kinetic['E_hall_y'] = E_hall_y
    dict_kinetic['E_TEwedge_y'] = E_TEwedge_y

    return dict_kinetic


#=======================================================================


def get_kinetic_q(path, time):
    '''
        This function gets all the kinetic components of heat flow
    22/2/19 - updated to include irregular Z profile
        formerly: extract_q  changed to---> get_kinetic_q inline with ' get_kinetic_E' (the Ohm's law counterpart)
    
    dict_kinetic = extract_q(path,time)
    
    
    
    '''
    time_int = int(time)
    fprefix = fpre(path)
    #--- load in data ....
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
    #x_grid_SI = x_grid*xstep_factor
    #y_grid_SI = y_grid*xstep_factor
    time_col = dict_fo['time'] * tstep_factor

    rA = (Z2ni)**-1

    #===============
    grid = dict['grid']
    dxfo, dyfo = cf.get_grad_3d_varZ(grid, fo)

    rho = rA
    omega = Bz * rA
    rho_vyx = extend_grid_xy_to_vxy(nv, ny, nx, rA)
    Bz_vyx = extend_grid_xy_to_vxy(nv, ny, nx, Bz)
    omega_vyx = extend_grid_xy_to_vxy(nv, ny, nx, omega)
    #--- kinetic
    dict = get_qSH_kinetic(grid, rho_vyx, Bz_vyx, jx, jy, fo)

    name_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']

    for name in name_list:
        dict[name]['data'] = np.transpose(dict[name]['data'])    #[:x_limit]
        dict[name]['label'] = convert_name_to_label(name)
    q_KSH_x, q_KSH_y, q_KRL_x, q_KRL_y, q_KTE_x, q_KTE_y, q_KE_x, q_KE_y = dict['q_SH_x'][
        'data'], dict['q_SH_y']['data'], dict['q_RL_x']['data'], dict['q_RL_y']['data'], dict[
            'q_TE_x']['data'], dict['q_TE_y']['data'], dict['q_E_x']['data'], dict['q_E_y']['data']
    #--- classical
    return dict
