'''
	This plot creaetes a single plot of the TSI growth rate + the Temperature amplitude

	Additional debug funtionalitiy is to plot dxT

'''
import sys
import numpy as np
sys.path.extend(["./"])
import matplotlib.pyplot as plt
import matplotlib.gridspec as GS
import matplotlib.ticker as ticker
from pylab import *

import PLOTTERS.gen_Te_lineout as TEL
import gen_dyqy_RL_varZ as dyqy
import gen_dBdt_bier_varZ as gdBdt
import MODULES.figure_prl_twocol as fprl
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
import MODULES.tsi_module as tsi
import MODULES.impact_norms as inorm
from MODULES.plot_utils import run_obj
import EH_poly_coeff_module as ep
# ---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
cB = 3.0 * np.sqrt(np.pi) / 4.0

# -----

# ---> file inputs
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir
norm_dir = paths.norm_dir
log_file = norm_dir + 'norm.log'
[T_ref, n_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir, Z_ref, Ar=6.51)

# ---- user inputs
# aesthetics
xmin, xmax = -10, 60.0    # max/min of x axis on plot
cmap = cm.viridis
style_list = ['--', '--', ':', '-', '--', ':']
mstyle_list = [None, None, None, 'x', '^', 'o']

# ---> TSI parameters
lambda_p_mfp = 5.0
s_list = [2]
bz_in = 400.0    # selected magnetic field [Tesla]
save_name = '%stsi_growth_LT%i_%i_dT_%iT.png' % (save_path, s_list[0], lambda_p_mfp, bz_in)

time_list = ['03', '05', '07']

path_list = []
leg_list = []
get_lt_lab = lambda il: r'$%i$' % il
for ip in range(len(s_list)):
    path_loc = paths.get_path(s_list[ip], bz_in, lambda_p_mfp)
    path_list.append(path_loc)
    leg_list.append(get_lt_lab(s_list[ip]))


# -----
# ==============================================================================
def get_grad(x_grid, data):
    grad = np.zeros((len(data)))
    grad[1:-1] = (data[2:] - data[:-2]) / (x_grid[2:] - x_grid[:-2])
    grad[0] = grad[1]
    grad[-1] = grad[-2]
    return grad


# ==============================================================================
def fpre(path):
    return path.split('/')[-1]


# -----------------------------------------------------------------------
# Compute growth rate
class GetSimData:

    def __init__(self, run_obj, **kwargs):
        path = run_obj.path
        time = kwargs.get('time', run_obj.time_max)
        iy = kwargs.get('iy', 20)
        self.cfg = run_obj.norm

        dict_load = cf.load_data_all(path, fpre(path), time)

        convert_array = lambda key: np.transpose(dict_load[key])[:, iy]
        convert_array_dy = lambda key: np.transpose(dict_load[key])[:, ::-1][:, iy]

        # data in

        self.ne_ref = self.cfg.n0
        self.Te_ref = self.cfg.T0
        self.Bz_ref = (m_e / (q_e * self.cfg.tau_ei)) * 1.0    # assumes that norm Bz is 1T
        self.v_te = self.cfg.v_te
        self.tau_ei = self.cfg.tau_ei
        self.Z_ref = self.cfg.Z
        self.Ar = self.cfg.Ar
        self.lambda_mfp = self.cfg.lambda_mfp
        if False:
            n_factor = self.ne_ref
            Te_factor = (2.0 * self.Te_ref)
            Bz_factor = self.Bz_ref
            C_factor = self.v_te

            P_factor = 2.0 / 3.0
            x_factor = self.cfg.lambda_mfp

        else:
            n_factor = 1.0
            Te_factor = 1.0
            Bz_factor = 1.0
            C_factor = 1.0
            P_factor = 2.0 / 3.0
            x_factor = 1.0
        # Only load in line outs of arrays
        # swap y direction after transpose for y gradients
        self.time = dict_load['time'] * self.tau_ei
        self.dyB = convert_array_dy('dyB')
        self.dyT = convert_array_dy('dyT')

        self.Te1 = cf.get_U_dev(dict_load['Te'])
        self.Bz1 = cf.get_U_dev(dict_load['Bz'])
        self.Z = convert_array('Z') * self.Z_ref
        self.U = convert_array('U') * P_factor
        self.Te = convert_array('Te') * Te_factor
        self.ne = convert_array('ne') * n_factor
        self.Cx = convert_array('Cx') * C_factor

        self.wt = convert_array('wt')
        self.Bz = convert_array('Bz') * Bz_factor
        self.x_grid = cf.trim_array_1D(dict_load['x_grid'], dict_load['nx'],
                                       dict_load['ny']) * x_factor

        self.dxT = convert_array('dxT') * Te_factor / x_factor
        self.dxB = convert_array('dxB') * Bz_factor / x_factor

        self.dxCx = get_grad(self.x_grid, self.Cx)
        self.dxn = get_grad(self.x_grid, self.ne)

        self.dx2T = get_grad(self.x_grid, self.dxT)
        self.dx2B = get_grad(self.x_grid, self.dxB)

    # ==============================================================================
    def load_transport(self, wte):
        # dict_n = inorm.impact_inputs_array(self.ne_ref, self.Te_ref, self.Z_ref, 1.0,
        #                                self.Ar)  # Bz must be in tesla because this is the units of gorgon data
        dict_n = inorm.impact_inputs_array(self.ne, self.Te, self.Z, self.Bz, self.Ar)
        # --------
        self.wpe_over_nu_ei = dict_n['wpe_over_nu_ei']
        self.c_over_vte = dict_n['c_over_vte']
        self.log_lambda = dict_n['log_lambda']
        self.lambda_T = dict_n['lambda_mfp']
        self.tau_T = dict_n['tau_ei']

        delta_c_norm = self.c_over_vte / self.wpe_over_nu_ei
        self.delta = delta_c_norm * self.lambda_T    # m
        dict_n['delta'] = self.delta

        ref_vals = {}
        ref_vals['n_e'] = self.ne_ref
        ref_vals['T_e'] = self.Te_ref
        ref_vals['tau_B'] = self.tau_T * cB

        if type(wte) != np.ndarray:
            wte = np.array([wte])
        self.alpha_para = np.zeros(wte.shape)
        self.alpha_perp = np.zeros(wte.shape)
        self.alpha_wedge = np.zeros(wte.shape)

        self.kappa_para = np.zeros(wte.shape)
        self.dkappawedge_dchi = np.zeros(wte.shape)
        self.dbetawedge_dchi = np.zeros(wte.shape)
        self.dbetaperp_dchi = np.zeros(wte.shape)
        self.dalphawedge_dchi = np.zeros(wte.shape)

        self.beta_wedge = np.zeros(wte.shape)
        self.psi_wedge = np.zeros(wte.shape)
        self.kappa_perp = np.zeros(wte.shape)

        for iw, wt_val in enumerate(wte):
            dict, grad_dict = tsi.get_dimensionful_transport_c(wt_val, self.Z_ref, ref_vals)

            self.alpha_para[iw] = dict['alpha_para']
            self.alpha_perp[iw] = dict['alpha_perp']
            self.alpha_wedge[iw] = dict['alpha_wedge']

            self.kappa_para[iw] = dict['kappa_para']
            self.dkappawedge_dchi[iw] = grad_dict['kappa_wedge']
            self.dbetawedge_dchi[iw] = grad_dict['beta_wedge']
            self.dbetaperp_dchi[iw] = grad_dict['beta_perp']
            self.dalphawedge_dchi[iw] = grad_dict['alpha_wedge']

            self.beta_wedge[iw] = dict['beta_wedge']    # /e?
            self.psi_wedge[iw] = dict['beta_wedge']    # ettinghausen term
            self.kappa_perp[iw] = dict['kappa_perp']


# --------------------
# Plot TSI growthrate
class get_WKB:

    def __init__(self, k, GetSimData):
        self.SimData = GetSimData
        self.k = k * (1.0 / self.SimData.lambda_mfp)

        self.x_grid = self.SimData.x_grid

        self.ix = np.argsort(np.abs(self.SimData.x_grid - 20.0e-6))[0]

        self.T1 = self.SimData.Te1
        self.Bz1 = self.SimData.Bz1

        self.T0 = self.SimData.Te
        self.n0 = self.SimData.ne
        self.Bz0 = self.SimData.Bz

        self.dxT0 = self.SimData.dxT
        self.dxn0 = self.SimData.dxn
        self.dxBz0 = self.SimData.dxB
        self.wt = self.SimData.wt

        GetSimData.load_transport(self.wt)
        self.dkappa_wedge_dchi = self.SimData.dkappawedge_dchi

    def get_lhs_i(self):
        u = (self.T0**2.5) * self.T1
        dT1 = cf.get_grad_dx2_1d(self.x_grid, u)
        return dT1

    def get_lhs_ii(self):
        return -(self.k**2) * self.T1

    def get_rhs_i(self, dT1c, dBz):
        self.r1 = 1j * self.k * ((self.T0**1.5) / self.n0) * self.dkappa_wedge_dchi * (
            self.dxn0 * self.Bz0 - self.n0 * self.dxBz0)
        return self.r1

    def get_rhs_ii(self, dT1c, dBz):
        self.r2 = -self.k * self.dkappa_wedge_dchi * (self.T0**1.5 / self.n0) * self.dxT0 * dBz
        return self.r2

    def get_r(self):
        self.r1 = 1j * self.k * ((self.T0**1.5) / self.n0) * self.dkappa_wedge_dchi * (
            self.dxn0 * self.Bz0 - self.n0 * self.dxBz0)
        # r1 *=0.0
        self.Ln = np.abs(self.dxn0 / self.n0)
        alpha = self.get_alpha()
        self.r2 = alpha * self.k * self.Ln * self.dkappa_wedge_dchi * (
            self.T0**1.5 / self.n0) * self.dxT0    # Note removing 1j assumes that dB \propto -1jdT
        return self.r1, self.r2

    def get_A(self, dT1c=1.0):
        r1, r2 = self.get_r()
        r = r1 + r2

        factor = (self.k**2 + r)**0.25

        A = (dT1c * self.T0[self.ix]**2.5 / (self.T0)) * (factor[self.ix] / factor)
        return A

    def assemble_matrix(self, dT1c, dBz, **kwargs):
        if 'index_lim' in kwargs:
            ix_min, ix_max = kwargs['index_lim'][0], kwargs['index_lim'][-1]

        else:
            ix_min, ix_max = 0, -1
        f_vec = self.get_rhs_i(1.0, dBz)
        f_vec = f_vec[ix_min:ix_max]
        dx = self.x_grid[self.ix + 1] - self.x_grid[self.ix]
        f_vec *= 0.0

        nx = len(f_vec)
        uim1 = (1.0 / (dx**2)) * np.ones((nx - 1))
        ui = -(np.ones(np.shape(f_vec)) * (k**2) + f_vec + (2.0 / (dx**2)))
        uip1 = (1.0 / (dx**2)) * np.ones((nx - 1))

        mat = np.diag(ui, k=0) + np.diag(uip1, k=1) + np.diag(uim1, k=-1)
        # rhs vector
        b = ui[:] * 0.01    #-1j*self.get_rhs_ii( dT1c, dBz)[ix_min:ix_max]

        ##b*=0.0
        # bcs - Dirichlet
        b[-1] -= (1.0 / (dx**2)) * dT1c

        xvec = np.linalg.solve(mat, b)
        return xvec

    #--- numerical sol...

    def int_from_crit(self, dx, data, ix):
        integrand = data * dx
        integral = np.cumsum(integrand[::-1]) - np.sum(integrand[ix:])
        return integral[::-1][:ix]

    def get_alpha(selfs):
        """ get ratio between dB and dT = alpha"""
        return 1.0


lambda_p = 5.0
sim_obj = run_obj(scale_length=1, bz_in=50.0, lambda_p=lambda_p, pert_amp='1p')
sim_data_obj = GetSimData(sim_obj)

k = (2.0 * np.pi / lambda_p)
wkb_obj = get_WKB(k, sim_data_obj)

ixmin = 100
ixmax = 210
dTlineout = np.max(np.abs(wkb_obj.T1), axis=0)
dBlineout = np.max(np.abs(wkb_obj.Bz1), axis=0)

dT1c = dTlineout[ixmax]
dBz = dBlineout[ixmax]
x_vec = wkb_obj.assemble_matrix(dT1c, -1.0 * dBz, index_lim=[ixmin, ixmax])

xgrid = wkb_obj.x_grid[ixmin:ixmax]
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.semilogy(xgrid, np.abs(x_vec))
ax.semilogy(xgrid, np.abs(x_vec.real), c='r', linestyle='--')
ax.semilogy(xgrid, np.abs(x_vec.imag), c='r', linestyle=':')
ax2.plot(xgrid, wkb_obj.T0[ixmin:ixmax], c='k')
plt.show()
