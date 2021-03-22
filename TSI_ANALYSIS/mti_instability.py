#---------> MTI INSTABILITY GROWTH RATE CALCULATION
import numpy as np, getpass, site

userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import impact_norms_py3 as ip3
import EH_poly_coeff_module_py3 as ep3
import chfoil_module_py3 as cf
import matplotlib.pyplot as plt
q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
c = 3.0e8

##ne_ref = 5.0e18

path_pre = '../PREMAG/2D_RUNS/'
path = path_pre + 'r5_v40_Z_FEOS_MODNE_5y_matchedf0_in_50T_E'
#path = path_pre + 'r5_v40_vmax20_vnonuni10_400T_10y_Zb_3'
#path = path_pre + 'r5_v40_vmax20_vnonuni10_400T_Zb'
#time = '06'
time = '13'
norm_name = 'r_5v37'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + '/norm.log'
[Te_ref, ne_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
Te_ref *= q_e    # convert to Joules
cd5 = cf.conv_factors_custom(norm_path, Z_ref, Ar=6.51)
cB = 3.0 * np.sqrt(np.pi) / 4.0
lambda_mfp_ref = cd5.lambda_mfp
v_th_ref = cd5.v_te
tau_ref = cd5.tau_ei
Bz_ref = cd5.Bz_ref
#----##


def get_wt_grad(wt,
                Z,
                var_list=[
                    'alpha_perp', 'alpha_wedge', 'beta_perp', 'beta_wedge', 'kappa_perp',
                    'kappa_wedge'
                ]):
    dwt = 0.01
    if wt > dwt:
        wtmin = wt - dwt
        wtmax = wt + dwt
    else:
        wtmin = wt
        wtmax = wt + 2.0 * dwt
    t_dict_p1, c_dict = ep3.coeff_poly_fit(wtmax, Z)
    t_dict_m1, c_dict = ep3.coeff_poly_fit(wtmin, Z)
    dict_out = {}
    for var in var_list:
        dict_out[var] = (t_dict_p1[var] - t_dict_m1[var]) / (2.0 * dwt)
    return dict_out


#----##


def load_all_SI(path, fprefix, time, iy, ix):

    dict_all = cf.load_data_all(path, fprefix, time)

    dict_out = {}
    xgrid = dict_all['x_grid'][1:-1] * lambda_mfp_ref
    ygrid = dict_all['y_grid'][1:-1] * lambda_mfp_ref

    ne = dict_all['ne']
    Te = dict_all['Te']
    Ue = dict_all['U']
    dxT = dict_all['dxT']
    dxn = dict_all['dxn']
    dxBz = dict_all['dxBz']

    wt = dict_all['wt']
    Bz = dict_all['Bz']
    Z2ni = dict_all['Z2ni']
    Z = dict_all['Z']
    v_th = ((2.0 * q_e * Te_ref) / m_e)**0.5

    dict_out['x_grid'] = xgrid
    dict_out['y_grid'] = ygrid
    dict_out['ne'] = ne * (ne_ref * 1e6)    # m^-3
    dict_out['Te'] = Te * (2.0 * Te_ref)
    dict_out['Ue'] = Ue
    dict_out['Bz'] = Bz * Bz_ref
    dict_out['dxT'] = dxT * (2.0 * Te_ref) / (lambda_mfp_ref)
    dict_out['dxn'] = dxn * (ne_ref * 1e6) / (lambda_mfp_ref)
    dict_out['dxBz'] = dxBz * Bz_ref / (lambda_mfp_ref)

    dict_out['wt'] = wt
    dict_out['Z'] = Z * Z_ref
    tau_ei_norm = ((2.0 * Te)**(1.5)) / Z2ni
    tau_ei = (((2.0 * Te)**1.5) / Z2ni) * tau_ref    # seconds

    dict_out['tau_ei'] = tau_ei

    dict_out['lambda_mfp'] = ((2.0 * Te)**0.5) * tau_ei_norm * v_th_ref

    grad_t = tau_ei * dxT
    print(' --: ', np.shape(xgrid), np.shape(ygrid), np.shape(grad_t), np.shape(Te))
    dtau_dxT, dydata = cf.get_grad(ygrid, xgrid, grad_t)
    #dtau_dxT = (grad_t[1:] - grad_t[:-1])/xgrid
    dict_out['grad_tau_T'] = dtau_dxT
    return dict_out


#----##


def get_growth(k, Z_in, lambda_T_in, tau_ei_in, n_e_in, T_e_in, Bz_in, wt_in, dxn_in, dxT_in,
               dxBz_in):
    '''
        gamma_m,gamma_p =  get_growth(k,Z_in,lambda_T_in,tau_ei_in,n_e_in,T_e_in,wt_in,dxn_in,dxT_in,dtau_dxT_in)
    '''

    v_te_in = lambda_T_in / tau_ei_in

    #dxT_in = dxT[iy,ix]
    #dxn_in = dxn[iy,ix]d
    #n_e_in = n_e[iy,ix]
    #Z_in = Z[iy,ix]
    #T_e_in = T_e[iy,ix]
    #wt_in = 0.4#wt[iy,ix]
    #dtau_dxT_in =dtau_dxT[iy,ix]
    #lambda_T_in = lambda_T[iy,ix]
    #tau_ei_in = dict_out['tau_ei'][iy,ix]

    wpe = n_e_in * q_e**2 / (epsilon0 * m_e)
    delta = c / wpe
    Lambda = lambda_T_in / delta

    transport_dict, c_dict = ep3.coeff_poly_fit(wt_in, Z_in)
    t_grad_dict = get_wt_grad(wt_in, Z_in)

    alpha_perp = transport_dict['alpha_perp']
    kappa_perp = transport_dict['kappa_perp']
    beta_wedge = transport_dict['beta_wedge']
    psi_wedge = beta_wedge * 1.0    # Maxwellian transport coefficients

    dkappawedge_dchi = t_grad_dict['kappa_wedge']
    #dbetawedgedchi = t_grad_dict['beta_wedge']
    K = lambda_T_in * k
    K2 = K**2
    ln = n_e_in / dxn_in
    lt = T_e_in / dxT_in
    lb = Bz_in / dxBz_in

    LN = ln / lambda_T_in
    LT = lt / lambda_T_in
    LB = lb / lambda_T_in

    gamma = 5.0 / 3.0
    #vs2 = gamma*(2.0/3.0)*Ue_in/ne_in
    Vs = np.sqrt(gamma * T_e_in * (k_b / m_e)) / v_te_in

    i = complex(0.0, 1.0)
    ##KG2 = dkappawedge_dchi*(cB**2) * (Lambda**2)/(2.0*LT*LN*alpha_para*kappa_para)
    #----> table 8.1 JB thesis
    DT = (cB / 3.0) * kappa_perp
    DR = alpha_perp / (cB * (Lambda**2))
    AN = (cB / (2.0 * wt_in)) * beta_wedge
    AE = ((2.0 * wt_in) / (3.0 * cB * (Lambda**2))) * psi_wedge
    Ck = (cB / 3.0) * wt_in * dkappawedge_dchi
    CG = cB / (2.0 * wt_in)
    sigmaB = 1.0

    SG = (4.0 * Ck * CG) / (LT * LN)
    LambdaG = (1.0 + (LN / LT))
    SP = 4.0 * AN * (sigmaB / LT) * Ck
    SE = 4.0 * AN * AE
    VB = sigmaB * (Ck / LT) * ((LT / LB) - (LT / LN))
    LambdaB = (1.0 - (LB / LN))

    A1 = (i * K * DT - VB) + i * K * DR
    A2 = (i * K * DT - VB) * i * K * DR + 0.25 * SG - 0.25 * SP * i * K + 0.25 * SE * K2
    B1 = (0.6 * (i * K * DT - VB) + i * K * DR)
    B2 = 0.6 * (i * K * DT - VB) * i * K * DR + (3.0 / 20.0) * LambdaG * SG - (
        3.0 / 20.0) * SP * i * K + (3.0 / 20.0) * SE * (K2)

    p = np.zeros((5), dtype=complex)
    p[0] = 1.0
    p[1] = A1
    p[2] = (A2 - Vs**2)
    p[3] = -B1 * Vs**2
    p[4] = -B2 * Vs**2

    #N = (cB/2.0)*dbetawedgedchi*((lambda_T_in**2)/(tau_ei_in*T_e_in))*dtau_dxT_in

    #gamma_p = 0.5*(-1.0*((DT + DR)*(K**2) - N) + np.sqrt(((DT + DR)*(K**2)-N)**2 + 4.0*DT*(K**2)*(DR*(KG2 - K**2) + N)))
    #gamma_m = 0.5*(-1.0*((DT + DR)*(K**2) - N) - np.sqrt(((DT + DR)*(K**2)-N)**2 + 4.0*DT*(K**2)*(DR*(KG2 - K**2) + N)))

    print('COEFFS = ', p[0], p[1], p[2], p[3], p[4])
    roots = np.roots(p)
    print('roots = ', roots)
    #return gamma_m,gamma_p
    return roots


#----##

fprefix = path.split('/')[-1]
#dict_all = cf.load_data_all(path,fprefix,time)

iy, ix = 20, 110

dict_out = load_all_SI(path, fprefix, time, iy, ix)
n_e = dict_out['ne']
lambda_T = dict_out['lambda_mfp']
Z = dict_out['Z']
wt = dict_out['wt']
dtau_dxT = dict_out['grad_tau_T']
T_e = dict_out['Te']
dxn = dict_out['dxn']
dxT = dict_out['dxT']
dxBz = dict_out['dxBz']
Bz = dict_out['Bz']
x_grid = dict_out['x_grid']

fig = plt.figure()

ix_min, ix_max = 70, 220
ax = fig.add_subplot(111)
ax.plot(x_grid, T_e[20, :])
ax.axvline(x_grid[ix_min])
ax.axvline(x_grid[ix_max])
#plt.show()

iy = 20
ng = len(range(ix_min, ix_max))
gamma_m_arr = np.zeros((ng))
gamma_p_arr = np.zeros((ng))
k = (2.0 * np.pi) / (100.0 * lambda_mfp_ref)    #7.0e-6#
root_arr = np.zeros((ng, 4), dtype=complex)
for ix in range(ix_min, ix_max):

    dxT_in = dxT[iy, ix]
    dxn_in = dxn[iy, ix]
    n_e_in = n_e[iy, ix]
    Z_in = Z[iy, ix]
    T_e_in = T_e[iy, ix]
    wt_in = wt[iy, ix]
    Bz_in = Bz[iy, ix]
    dxBz_in = dxBz[iy, ix]

    #  dtau_dxT_in =dtau_dxT[iy,ix]
    lambda_T_in = lambda_T[iy, ix]
    tau_ei_in = dict_out['tau_ei'][iy, ix]
    root_in = get_growth(k, Z_in, lambda_T_in, tau_ei_in, n_e_in, T_e_in, Bz_in, wt_in, dxn_in,
                         dxT_in, dxBz_in)
    idx = np.argsort(np.abs(root_in))
    np.shape(root_in)
    root_arr[ix - ix_min, :] = root_in[idx]
    #gamma_m_arr[ix-ix_min],gamma_p_arr[ix-ix_min] =  get_growth(k,Z_in,lambda_T_in,tau_ei_in,n_e_in,T_e_in,wt_in,dxn_in,dxT_in,dtau_dxT_in)

fig = plt.figure()
ax_list = [fig.add_subplot(141), fig.add_subplot(142), fig.add_subplot(143), fig.add_subplot(144)]

for iax in range(len(ax_list)):
    ax = ax_list[iax]
    ax.plot(x_grid[ix_min:ix_max], root_arr[:, iax].real)
    ax.plot(x_grid[ix_min:ix_max], root_arr[:, iax].imag, linestyle='--')
    ax.axhline(0.0, linestyle=':')
    ax.set_ylim(-1e7, 1e7)

plt.show()
'''
lt = T_e[iy,ix_min:ix_max]/dxT[iy,ix_min:ix_max]
ln = n_e[iy,ix_min:ix_max]/dxn[iy,ix_min:ix_max]
fig = plt.figure()
ax = fig.add_subplot(211)
axg = fig.add_subplot(212)
#axg2 = axg.twinx()
axg2 = ax.twinx()
ax.plot(x_grid[ix_min:ix_max],T_e[20,ix_min:ix_max],c='r')
axg.plot(x_grid[ix_min:ix_max],gamma_m_arr,c='k')
axg.plot(x_grid[ix_min:ix_max],gamma_p_arr,c='k',linestyle='--')
axg.axhline(0.0,linestyle=':',c='b')
axg2.plot(x_grid[ix_min:ix_max],ln*lt)
plt.show()
'''
