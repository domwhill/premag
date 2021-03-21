'''
    Calculates Braginksii heat flows for all the different time steps for fil in 'path' 
    Saves as .npy files.
    
    Formerly 'q_SH_Te_tsteps.py'
    
    Redundant packages functions removed. (Nernst + Biermann calculations)
    
    19/02/2019 - updated to include varZ
    
'''
import os, site, getpass, re, sys
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')
import EH_poly_coeff_module as sv
import MODULES.chfoil_module as cf
import matplotlib as mpl

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
#---------------
# n0 = (1.87e22) # [cm^-3]
# T0 = 850 #eV

#os.getcwd()
cwdpath = os.getcwd() + '/' + sys.argv[-1]
#cd5 = cf.conv_factors_custom(cwdpath)
#cd5 = cf.conv_factors_premag()
#cd5 = cf.conv_factors_eos()
#cd5 = cf.conv_factors_custom()
Z0 = cf.get_atomic_Z(os.getcwd() + '/fort.10')
nv, ny, nx = cf.get_dims(os.getcwd() + '/fort.10')

print('Z0 = ', Z0)
print('nv=%i ny=%i nx = %i' % (nv, ny, nx))
#T0 = cd5.T0
#n0 = cd5.n0
#logLambda = cd5.logLambda

SI_on = True
mult = 6
yi = 4
xi = 32
time = '04'
t_list = ['01', '04', '05', '06', '07', '08', '09', '10', '12', '13', '14', '15', '16']

#===================================
ones_nv = np.ones((nv))
#vmax = 10.0
#xmin,xmax,nx = 0.0,289.4736,192
#ymin,ymax,ny = -5.0,45.0,50
#vmin,vmax,nv = 0.0,10.0,100
#v_grid = np.arange(vmin,vmax,(vmax-vmin)/nv)
path = os.getcwd()
fprefix = path.split('/')[-1]
dict_in = cf.load_dict(path, fprefix, 'fo', '01')
v_grid = dict_in['v_grid']


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


dv = calc_dv(v_grid)
#====================================

print 'Z0 = ', Z0
#print 'T0 = ',T0
#print 'n0 = ',n0
#print 'v_te = ', v_te
xstep_factor = 1.0
tstep_factor = 1.0
#------------------------------------------------


def list_to_float(input):
    list = input.split()
    arr = np.zeros((len(list)), dtype=float)
    arr[:] = list[:]
    return arr


#------------------------------------------------


def fpg_get_info(fname):
    '''
        Gets the IMPACT header info
        time,ndims,dim_array = fpg_get_info(fname)
        dim_Aray = [nv,ny,nx] or [ny,nx]
    '''

    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    mat = np.loadtxt(fname, skiprows=10)

    #print mat_info[0]
    time = float(mat_info[1])
    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])

        nx = int(mat_info[5])
        x_grid = np.array(mat_info[6])
        dim_array = [ny, nx]
        #put this in a dictionary?
        #grid_array = [
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])
        dim_array = [nv, ny, nx]
    return time, ndims, dim_array


def construct_fname(path, fprefix, var, time):

    if var == 'fo' or var == 'fxX' or var == 'fxY' or var == 'fyX' or var == 'fyY':
        suffix = '.xyv'
    else:
        suffix = '.xy'

    #fprefix = 'thydro_hi'
    #time = '00'
    fname = path + '/' + fprefix + '_' + var + '_' + time + suffix
    return fname


def get_startline(fname):
    '''
        get the line where the data starts
    '''
    f = open(fname)
    f_list = f.readlines()
    out_line = 9
    for i in range(len(f_list)):
        if re.match(r'\n', f_list[i]):
            ##print 'found match on line: ', i
            out_line = i
            #else:
            #out_line = 9
            #continue
            break

    return out_line


def get_index(ny, nx):
    index = ny + 2 + 1
    return index


def get_time(fname):
    f = open(fname)
    data = f.readlines()
    t = float(data[1]) * tstep_factor
    return t


#------------------------------------------------
def fpg_load(fname):
    '''
        loads the file
        returns hmat,v_grid,y_grid,x_grid = fpg_load(fname)
        now can deal with fo when halo cells are also saved
    '''
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    mat = np.loadtxt(fname, skiprows=10)
    out_l = get_startline(fname)
    ##print mat_info
    time = float(mat_info[1])
    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])
        #print '========'
        ##print 'nx: ',nx,'\nx_grid: \n',x_grid
        #print 'np.shape: xgrid: ', np.shape(x_grid)
        #print '============', 'out_l = ====', out_l
        mat = np.loadtxt(fname, skiprows=8)
        #print np.shape(mat)
        grid = [y_grid, x_grid]
        hmat_final = mat
        v_grid = 0.0
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])

        mat = np.loadtxt(fname, skiprows=out_l)
        dim_array = [nv, ny, nx]
        dim_reverse = [nx, ny, nv]

        mat = mat.reshape(dim_reverse)
        hmat_final = np.transpose(mat)
        grid = [v_grid, y_grid, x_grid]    #nv,ny,nx

        #print '############################'
        #print '   nv, ny, nx: ',nv, ny, nx
        #print '############################'
    return hmat_final, v_grid, y_grid, x_grid


def get_v2dv(v_grid):
    #dv = v_grid[1:] - v_grid[:-1]
    dv = calc_dv(v_grid)
    v2dv = (v_grid**2) * dv
    return v2dv


def maxw_dist(v, vte):
    #print ' vte = ', vte
    # --- hang on a sec here ... v should already be normalised with respect to vte?
    #sys.exit()
    #f_om = ((np.pi*(vte**2))**(-1.5)) * np.exp(-(v/vte)**2)
    f_om = ((np.pi * (vte**2))**(-1.5)) * np.exp(-((v)**2) * (vte**-2))
    #f_om =  np.exp(-((v)**2))
    #f_om =  ((np.pi*(vte**2))**(-1.5)) * np.exp(-((v)**2))

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
def get_beta_perp(w, rZZni, ne, Te, F0, v_grid):
    '''
        beta_perp = get_beta_perp(w,rho,ne,Te,F0,v_grid)
    '''
    #rZZni = rho
    vte = (2.0 * Te)**0.5
    omega = w * rZZni
    delta = get_delta(v_grid, omega, F0)
    v_mom_7_5 = get_v_mom_m_n(v_grid, omega, F0, 7, 5)
    v_mom_8_5 = get_v_mom_m_n(v_grid, omega, F0, 8, 5)
    v_mom_10_8 = get_v_mom_m_n(v_grid, omega, F0, 10, 8)

    #dv = v_grid[1] - v_grid[0]
    test_mom_7 = np.sum(np.exp(-(v_grid / vte)**2) * (v_grid**9) * dv)
    test_mom_5 = np.sum(np.exp(-(v_grid / vte)**2) * (v_grid**7) * dv)

    ###print ' v_mom_7_5 = ', v_mom_7_5,'\n\n tes_mom_7_5 = ', test_mom_7/test_mom_5, test_mom_7, test_mom_5
    beta = (1.0 / (2.0 * Te)) * (v_mom_7_5 + ((omega * v_mom_8_5)**2) * v_mom_10_8) / delta - 2.5
    #---- NOTE 30/11/17 found error in beta_perp, missing factor (1.0/2.0*Te)
    return beta


#-----------------------------------------------------------------------
def get_beta_wedge(w, rho, ne, Te, F0, v_grid):
    '''
        beta_wedge = get_beta_wedge(w,rho,ne,Te,F0,v_grid)
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
def get_q_SH(Z2ni, ne, Te, w, dxT, dyT, jx=0.0, jy=0.0):
    '''
        MAIN ROUTINE - calculates the classical heat flow components
    '''

    #Te = 0.5
    vte = (2.0 * Te)**0.5
    #w = 0.1
    #ne = 1.0
    rho = Z2ni**-1

    #print ' ######################################'
    #print 'inputs : vte: %3.4f \t vmax: %3.4f \t \n \t w = %3.4f \n\t ne =  %3.4f \t rho %3.4f ' % (vte,vmax, w, ne ,rho)
    #print ' ######################################'
    ##print ' v_grid ===== ', v_grid
    # ------

    F0 = maxw_dist(v_grid, vte)

    #alpha_perp = get_alpha_perp(w,rho,ne,Te,F0,v_grid)
    ##print ' ap'
    #alpha_wedge = get_alpha_wedge(w,rho,ne,Te,F0,v_grid)
    ##print ' aw'
    beta_perp_c = get_beta_perp(w, rho, ne, Te, F0, v_grid)
    ##print 'bp'
    beta_wedge_c = get_beta_wedge(w, rho, ne, Te, F0, v_grid)
    ##print 'bw'
    #beta_wedge_c,beta_perp_c = 0.0,0.0

    beta_perp = beta_perp_c * (Te)
    beta_wedge = beta_wedge_c * (Te)

    kappa_perp_c = get_kappa_perp(w, rho, ne, Te, F0, v_grid)
    ##print 'kp'
    kappa_wedge_c = get_kappa_wedge(w, rho, ne, Te, F0, v_grid)

    #kappa_perp = kappa_perp_c*((2.0*Te)**2.5)
    #kappa_wedge = kappa_wedge_c*((2.0*Te)**2.5)

    kappa_perp = kappa_perp_c
    kappa_wedge = kappa_wedge_c

    #print 'kappa_perp: ', kappa_perp, 'kappa_wedge = ', kappa_wedge, 'dxT = ', dxT, ' dyT = ', dyT
    ##print ' v_grid = ', v_grid
    ##print ' F0 = ', F0

    q_SH_x = -kappa_perp * dxT
    q_SH_y = -kappa_perp * dyT
    q_RL_x = +kappa_wedge * dyT
    q_RL_y = -kappa_wedge * dxT

    #q_x = q_SH_x + q_RL_x

    #q_E_x = -(beta_perp*jx - beta_wedge*jy)
    #q_E_y = -(beta_perp*jy + beta_wedge*jx)
    q_E_x = -(beta_perp * jx - beta_wedge * jy)
    q_E_y = -(beta_perp * jy + beta_wedge * jx)

    #print 'q_SH_x: ', q_SH_x,'q_RL_x = ', q_RL_x,'q_E_x = ', q_E_x

    return q_SH_x, q_RL_x, q_E_x, q_SH_y, q_RL_y, q_E_y


#-----------------------------------------------------------------------
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
        dxT = np.zeros((len(dx), len(dy)))
        dyT = np.zeros((len(dx), len(dy)))
        ldx, ldy = len(dx), len(dy)
        #print 'SHAPE DX DY = ', ldx,ldy
        for ii in range(len(dx)):
            #print 'ii=  ', ii
            dyT[ii, :] = (T_data[ii, 2:] - T_data[ii, :-2]) / dy

        for ii in range(len(dy)):
            #print 'ix = ', ii, np.shape((T_data[2:,ii] - T_data[:-2,ii])), np.shape(dxT), np.shape(dx)
            dxT[:, ii] = (T_data[2:, ii] - T_data[:-2, ii]) / dx
        return dxT, dyT


def get_grad_centred(x_grid, y_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
    '''

    if len(np.shape(T_data)) == 1:
        dx = x_grid[1:] - x_grid[:-1]
        dxT = (T_data[1:] - T_data[:-1]) / dx
        dxT = dxT[1:]
        return dxT
    else:
        dx = x_grid[1:] - x_grid[:-1]
        dy = y_grid[1:] - y_grid[:-1]
        dxT = np.zeros((len(dx), len(dy)))
        dyT = np.zeros((len(dx), len(dy)))
        ldx, ldy = len(dx), len(dy)
        #print 'SHAPE DX DY = ', ldx,ldy
        for ii in range(len(dx)):
            #print 'ii=  ', ii
            dyT[ii, :] = (T_data[ii, 1:] - T_data[ii, :-1]) / dy

        for ii in range(len(dy)):
            #print 'ix = ', ii, np.shape((T_data[2:,ii] - T_data[:-2,ii])), np.shape(dxT), np.shape(dx)
            dxT[:, ii] = (T_data[1:, ii] - T_data[:-1, ii]) / dx

        dyT = dyT[1:, :-1]
        dxT = dxT[:-1, 1:]
        return dxT, dyT


def get_q_poly(Z2ni, n, T, B, dxT, dyT, Z=6.51):
    '''
        tau_ei_norm = T_norm**1.5/n_norm
        chi = B*tau_ei_norm
        q_x_SH,q_x_RL,q_y_SH,q_y_RL = get_q_poly(n,T,B,dxT,dyT,Z=6.51)
        
        THIS IS WRONG AND NEEDS TO BE CHANGED --------
    '''
    chi = B * (3.0 * (np.pi**0.5) / 4.0) * ((2.0 * T)**1.5) / (Z2ni)    #B*((T**1.5)/n)

    t_dict, coeff_dict = sv.coeff_poly_fit(chi, Z)

    kappa_perp_c = t_dict['kappa_perp'] * (2**1.5)
    kappa_wedge_c = t_dict['kappa_wedge']
    # Ridgers 2.31
    kappa_coeff = (T**2.5)    #n_e*(T_e_eV*q_e)*t_ei/(m_e)
    alpha_coeff = (T**-1.5)

    kappa_perp = kappa_coeff * kappa_perp_c
    kappa_wedge = kappa_coeff * kappa_wedge_c
    q_x_SH = -kappa_perp * dxT
    q_x_RL = +kappa_wedge * dyT
    q_y_SH = -kappa_perp * dyT
    q_y_RL = -kappa_wedge * dxT
    return q_x_SH, q_x_RL, q_y_SH, q_y_RL


if __name__ == "__main__":
    '''
        dxT = T_data[1:] - T_data[:-1]/(x_grid[1:] - x_grid[:-1])
    '''

    path = os.getcwd()
    var_list = ['Cx', 'Cy', 'n', 'Te', 'ExX', 'EyY', 'jxX', 'jyY', 'Bz', 'qxX', 'qyY']
    var_list = ['n', 'Te', 'Cx', 'ExX', 'Bz', 'qxX']

    tc_list = []
    t_list = []
    list_files = os.listdir(path)
    leg_list = []
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var_list[0] + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')

            t_list.append(t_l)
            t_col = get_time(fname)
            tc_list.append(t_col)
            ##print t_list
            fprefix = s.group('fpre')

    fprefix = path.split('/')[-1]
    ftest = construct_fname(path, fprefix, 'n', '00')
    hmat_final, f_grid, y_grid, x_grid = fpg_load(ftest)
    xstep = (x_grid[1:] - x_grid[:-1]) * xstep_factor
    ystep = (y_grid[1:] - y_grid[:-1]) * xstep_factor

    T_name = construct_fname(path, fprefix, 'Te', time)
    T_data, f_grid, y_grid, x_grid = fpg_load(T_name)
    dxT, dyT = get_gradT(x_grid, y_grid, T_data)

    n0, n1, n2 = len(t_list), len(dxT[:, 0]), len(dxT[0, :])
    z_arr = np.zeros((n0, n1, n2))
    q_SH_x = np.zeros((n0, n1, n2))
    q_SH_y = np.zeros((n0, n1, n2))

    q_RL_x = np.zeros((n0, n1, n2))
    q_RL_y = np.zeros((n0, n1, n2))

    q_SH_xp = np.zeros((n0, n1, n2))
    q_SH_yp = np.zeros((n0, n1, n2))

    q_RL_xp = np.zeros((n0, n1, n2))
    q_RL_yp = np.zeros((n0, n1, n2))

    q_E_x = np.zeros((n0, n1, n2))
    q_E_y = np.zeros((n0, n1, n2))

    qx_name = construct_fname(path, fprefix, 'qxX', t_list[0])
    qy_name = construct_fname(path, fprefix, 'qyY', t_list[0])
    q_xX_0, f_grid, y_grid, x_grid = fpg_load(qx_name)
    q_yY_0, f_grid, y_grid, x_grid = fpg_load(qy_name)
    x1, x2 = np.shape(q_xX_0)
    y1, y2 = np.shape(q_yY_0)

    q_xX = np.zeros((n0, x1, x2))
    q_yY = np.zeros((n0, y1, y2))

    q_x = np.zeros((n0, n1, n2))
    q_y = np.zeros((n0, n1, n2))

    for tt in range(len(t_list)):
        time = t_list[tt]
        T_name = construct_fname(path, fprefix, 'Te', time)
        n_name = construct_fname(path, fprefix, 'n', time)
        qx_name = construct_fname(path, fprefix, 'qxX', time)
        qy_name = construct_fname(path, fprefix, 'qyY', time)
        jx_name = construct_fname(path, fprefix, 'jxX', time)
        jy_name = construct_fname(path, fprefix, 'jyY', time)
        Bz_name = construct_fname(path, fprefix, 'Bz', time)
        Z_name = construct_fname(path, fprefix, 'Z', time)
        Z2ni_name = construct_fname(path, fprefix, 'Z2niX', time)
        t_time = get_time(T_name)
        leg_list.append(t_time)
        q_xX_0, f_grid, y_grid, x_grid = fpg_load(qx_name)
        q_yY_0, f_grid, y_grid, x_grid = fpg_load(qy_name)

        jx, f_grid, y_grid, x_grid = fpg_load(jx_name)
        jy, f_grid, y_grid, x_grid = fpg_load(jy_name)
        Bz_data, f_grid, y_grid, x_grid = fpg_load(Bz_name)
        T_data, f_grid, y_grid, x_grid = fpg_load(T_name)
        Z_data, f_grid, y_grid, x_grid = fpg_load(Z_name)
        n_data, f_grid, y_grid, x_grid = fpg_load(n_name)

        q_xX[tt, :, :] = q_xX_0
        q_yY[tt, :, :] = q_yY_0

        dxT, dyT = get_gradT(x_grid, y_grid, T_data)

        for ix in range(1, len(T_data[0, :]) - 1):
            for iy in range(1, len(T_data[:, 0]) - 1):
                #print 'tt== ', tt, 'ix: = ', ix,'iy = ', iy, np.shape(T_data), np.shape(dxT),np.shape(n_data), np.shape(Bz_data),np.shape(dyT),'jx: ', np.shape(jx),'jy: ',np.shape(jy)
                Z2ni = Z_data[iy, ix] * n_data[iy, ix]
                q_SH_xp, q_RL_xp, q_E_xp, q_SH_yp, q_RL_yp, q_E_yp = get_q_SH(
                    Z2ni, n_data[iy, ix], T_data[iy, ix], Bz_data[iy - 1, ix - 1],
                    dxT[iy - 1, ix - 1], dyT[iy - 1, ix - 1], jx[iy - 1, ix - 1], jy[iy - 1,
                                                                                     ix - 1])
                q_SH_x[tt, iy - 1, ix - 1] = q_SH_xp
                #print ' q_SH_x: ', q_SH_x[tt,iy-1,ix-1]
                q_RL_x[tt, iy - 1, ix - 1] = q_RL_xp
                #print 'RL  q_SH_x: ', q_SH_x[tt,iy-1,ix-1]

                q_E_x[tt, iy - 1, ix - 1] = q_E_xp
                #print ' E q_SH_x: ', q_SH_x[tt,iy-1,ix-1]

                q_SH_y[tt, iy - 1, ix - 1] = q_SH_yp
                #print 'SH y q_SH_x: ', q_SH_x[tt,iy-1,ix-1]

                q_RL_y[tt, iy - 1, ix - 1] = q_RL_yp
                #print ' RL y q_SH_x: ', q_SH_x[tt,iy-1,ix-1]

                q_E_y[tt, iy - 1, ix - 1] = q_E_yp
    #--------- end loop

    np.save('q_xX.txt', q_xX)
    np.save('q_yY.txt', q_yY)
    np.save('q_SH_x.txt', q_SH_x)
    np.save('q_RL_x.txt', q_RL_x)
    np.save('q_SH_y.txt', q_SH_y)
    np.save('q_RL_y.txt', q_RL_y)
    np.save('q_E_x.txt', q_E_x)
    np.save('q_E_y.txt', q_E_y)

    #------------------------------
    q_x = q_SH_x + q_RL_x + q_E_x
    f_lim = 0.15
    ##print q_SH_x
    #print ' np.shape(x_grid): ', np.shape(x_grid), np.shape(q_SH_x),np.shape(q_xX)
    fig1, ax1 = plt.subplots(1, 1)
    #q_F_x = -0.5*n_data*((T_data)**1.5)*f_lim # 2.52 Ridgers, 0.5*ne*me*vt**3
    x_grid = x_grid * xstep_factor
    y_grid = y_grid * xstep_factor
    for tt in range(len(t_list)):
        ax1.plot(x_grid[1:-1], q_SH_x[tt, :, ny / 2], c='r', linestyle='-.')
        ax1.plot(x_grid[1:-1], q_RL_x[tt, :, ny / 2], c='r', linestyle='--')
        ax1.plot(x_grid[1:-1], q_E_x[tt, :, ny / 2], c='r', linestyle=':')
        ax1.plot(x_grid[1:-1], q_x[tt, :, ny / 2], c='r')
        #ax1.plot(x_grid[1:-1],q_SH_y[tt,:,ny/2],c='g')
        #ax1.plot(x_grid[1:-1],q_RL_y[tt,:,ny/2],c='g',linestyle='--')
        ax1.plot(x_grid[1:-1], q_xX[tt, :-1, ny / 2], c='b', linestyle='-')

    ax1.set_xlabel('x pos')

    ax1.set_ylabel('q_SH')
    #ax2 = ax.twinx()
    #ax2.plot(x_grid,Bz_data[:,ny/2],c='b')

    q_x = q_SH_x + q_RL_x + q_E_x
    #q_x = q_SH_x
    #q_xp = q_SH_xp + q_RL_xp
    #q_y = q_SH_y
    q_y = q_SH_y + q_RL_y + q_E_y
    #q_yp = q_SH_yp + q_RL_yp

    fig, ax = plt.subplots(1, 2)
    n = len(t_list)
    color = iter(cm.seismic(np.linspace(0, 1, n)))
    for tt in range(len(t_list)):
        c = next(color)
        ax[0].plot(x_grid[1:-1],
                   q_xX[tt, :-1, ny / 2] / q_x[tt, :, ny / 2],
                   c=c,
                   label=r'\textit{$' + str(leg_list[tt]) + '$}')

    #ax.plot(x_grid[1:-1],q_xX[tt,:-1,ny/2]/q_xp[tt,:,ny/2])
    #ax.plot(x_grid[1:-1],q_F_x[1:-1,ny/2],linestyle=':')
    #leg = plt.legend(title =leg_title)
    #leg.draggable(True)
    #setp(leg.get_title(),fontsize=25)
    ax[0].set_ylabel('$q_x/q_{x,Braginskii}$')
    ax[0].set_xlabel(xlab)

    #fig,ax = plt.subplots(1,1)
    color = iter(cm.seismic(np.linspace(0, 1, n)))
    for tt in range(len(t_list)):
        c = next(color)
        ax[1].plot(y_grid[1:-1],
                   q_yY[tt, xi, :-1] / q_y[tt, xi, :],
                   c=c,
                   label=r'\textit{$' + str(leg_list[tt]) + '$}')
        #ax.plot(y_grid[1:-1],q_yY[tt,xi,:-1]/q_yp[tt,xi,:])

        #ax.plot(x_grid[1:-1],q_F_x[1:-1,ny/2],linestyle=':')

    ax[1].set_ylabel('$q_y/q_{y,Braginskii}$')
    ax[1].set_xlabel('y')
    #plt.show()

    fig, ax = plt.subplots(1, 2)
    for tt in range(len(t_list)):
        ax[0].plot(x_grid[1:-1], q_E_x[tt, :, ny / 2] / q_SH_x[tt, :, ny / 2], c='b', label='E/SH')
        ax[0].plot(x_grid[1:-1],
                   q_RL_x[tt, :, ny / 2] / q_SH_x[tt, :, ny / 2],
                   c='r',
                   label='RL/SH')

    for tt in range(len(t_list)):
        ax[1].plot(y_grid[1:-1], q_E_y[tt, xi, :] / q_SH_y[tt, xi, :], c='b', label='E/SH')
        ax[1].plot(y_grid[1:-1], q_RL_y[tt, xi, :] / q_SH_y[tt, xi, :], c='r', label='RL/SH')

    leg = plt.legend()
    plt.clf()
