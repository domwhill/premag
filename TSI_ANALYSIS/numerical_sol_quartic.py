'''
	This plot creaetes a single plot of the TSI growth rate + the Temperature amplitude

	Additional debug funtionalitiy is to plot dxT

'''
import sys

sys.path.extend(["./"])
from pylab import *

import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
import MODULES.tsi_module as tsi
import MODULES.impact_norms as inorm
from MODULES.plot_utils import RunInfo

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
cd5 = cf.ConversionFactors(norm_dir, Z_ref, Ar=6.51)

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


def plot_complex(axy, x, **kwargs):
    p0, = axy.plot(x.real, **kwargs)
    axy.plot(x.imag, linestyle='--', **kwargs)
    return p0


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
class SimCoeffs:

    def __init__(self, k, GetSimData):
        self.SimData = GetSimData
        self.k = k * (1.0 / self.SimData.lambda_mfp)

        self.x_grid = self.SimData.x_grid

        self.ix = np.argsort(np.abs(self.SimData.x_grid - 20.0e-6))[0]

        self.T1 = self.SimData.Te1
        self.Bz1 = self.SimData.Bz1

        self.T0 = self.SimData.Te
        self.n0 = self.SimData.ne
        self.Cx0 = self.SimData.Cx

        self.Bz0 = self.SimData.Bz

        self.dxT0 = self.SimData.dxT
        self.dxn0 = self.SimData.dxn
        self.dxBz0 = self.SimData.dxB
        self.dxCx0 = self.SimData.dxCx

        self.dx2T0 = self.SimData.dx2T
        self.dx2Bz0 = self.SimData.dx2B
        self.wt = self.SimData.wt
        self.taub = self.SimData.tau_ei

        GetSimData.load_transport(self.wt)
        self.dkappa_wedge_dchi = self.SimData.dkappawedge_dchi

    def get_quartic_coefficients(self, Kx):
        '''
        g1,g2,g3, a1r,a1i,a2,b1r,b1i, b2r, b2i = self.get_quartic_coefficients( Kx)
        :param Kx:
        :return:
        '''
        n0 = self.n0
        Cx0 = self.Cx0
        Te0 = self.T0
        Bz0 = self.Bz0
        dxTe0 = self.dxT0
        dxBz0 = self.dxBz0
        dxn0 = self.dxn0
        dxCx0 = self.dxCx0
        d2xTe0 = self.dx2T0
        d2xBz0 = self.dx2Bz0

        dc2 = 0.1
        me = 9.11e-31
        k = 1.0
        dkappwedgedchi = self.dkappa_wedge_dchi
        omega1 = 0.0

        ap0 = self.SimData.alpha_perp
        kp0 = self.SimData.kappa_perp
        dbetawedgedchi = self.SimData.dbetawedge_dchi
        taub = self.taub
        a1i = (1 / (2 * me * n0 ** 2)) * (3 * me * n0 ** 3 * omega1 + 2 * dkappwedgedchi * k * taub ** 2 * n0 * \
                                          Te0 ** 4 * dxBz0 - 2 * dkappwedgedchi * k * taub ** 2 * Bz0 * \
                                          Te0 ** 4 * dxn0)

        a1r = (1 / (4 * me)) * (kp0 * taub * (Te0 ** 2) ** (1 / 4.0) * (4 * (k - Kx) * (k + Kx) * Te0 ** 2 - \
                                                                      15 * dxTe0 ** 2 - 10 * Te0 * (2 * Kx * dxTe0 + \
                                                                                                    d2xTe0)))
        a2 = -((1j * dkappwedgedchi * k * taub**2 * Te0**4 * dxTe0) / (me * n0))


        b1r = -((1 / (4 * n0 ** 2 * Te0 ** 4)) * ((Te0 ** 2) ** (1 / 4.0) * \
                                                  (n0 * (2 * dbetawedgedchi * taub * Te0 ** 4 * dxBz0 * \
                                                         (2 * Kx * Te0 + 3 * dxTe0) + 3 * ap0 * dc2 * n0 * \
                                                         (dxBz0 * (2 * Kx * Te0 - 5 * dxTe0) + \
                                                          2 * Te0 * d2xBz0)) + dbetawedgedchi * taub * Bz0 * \
                                                   Te0 ** 3 * (-2 * Te0 * dxn0 * (2 * Kx * Te0 + \
                                                                                  3 * dxTe0) + n0 * (
                                                                           4 * (-k ** 2 + Kx ** 2) * Te0 ** 2 + \
                                                                           3 * dxTe0 ** 2 + 6 * Te0 * (2 * Kx * dxTe0 + \
                                                                                                       d2xTe0))))))

        b1i = (k * dxn0) / n0
        b2r = (1 / (2 * n0 ** 2 * Te0 ** 4)) * (2 * Kx * Cx0 * n0 ** 2 * Te0 ** 4 + \
                                                2 * dbetawedgedchi * taub * Te0 ** 4 * (Te0 ** 2) ** (3.0 / 4.0) * dxn0 * \
                                                dxTe0 + n0 ** 2 * (2 * ap0 * dc2 * (k ** 2 + Kx ** 2) * Te0 * \
                                                                   (Te0 ** 2) ** (3.0 / 4.0) + 2 * Te0 ** 4 * dxCx0 - \
                                                                   3 * ap0 * dc2 * Kx * (Te0 ** 2) ** (3.0 / 4.0) * dxTe0) - \
                                                dbetawedgedchi * taub * n0 * Te0 ** 3 * (Te0 ** 2) ** (3.0 / 4.0) * \
                                                (3 * dxTe0 ** 2 + 2 * Te0 * (Kx * dxTe0 + \
                                                                             d2xTe0)))
        b2i = -omega1
        # g3 Kx^2 + g2 Kx + g1
        g1 = (ap0 * dc2 * k ** 2 * (Te0 ** 2.) ** (3.0 / 4.0)) / Te0 ** 3 + dxCx0 + \
             (dbetawedgedchi * taub * (Te0 ** 2.) ** (3.0 / 4.0) * dxn0 * \
              dxTe0) / n0 ** 2 - \
             (3. * dbetawedgedchi * taub * (Te0 ** 2.) ** (3.0 / 4.0) * dxTe0 ** 2) / \
             (2. * n0 * Te0) - (dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * d2xTe0) / n0
        g2 = Cx0 - (dbetawedgedchi * taub *
                    (Te0**2)**(3.0 / 4.0) * dxTe0) / n0 - (3 * ap0 * dc2 * (Te0**2)**
                                                           (3.0 / 4.0) * dxTe0) / (2. * Te0**4)
        g3 = (ap0 * dc2 * (Te0**2)**(3.0 / 4.0)) / Te0**3
        return g1, g2, g3, a1r, a1i, a2, b1r, b1i, b2r, b2i

    def get_uv_quartic_coefficients(self, Kx):
        '''
        g1,g2,g3, a1r,a1i,a2,b1r,b1i, b2r, b2i = self.get_quartic_coefficients( Kx)
        :param Kx:
        :return:
        '''
        n0 = self.n0
        Cx0 = self.Cx0
        Te0 = self.T0
        Bz0 = self.Bz0
        dxTe0 = self.dxT0
        dxBz0 = self.dxBz0
        dxn0 = self.dxn0
        dxCx0 = self.dxCx0
        d2xTe0 = self.dx2T0
        d2xBz0 = self.dx2Bz0

        dc2 = 0.1
        me = 9.11e-31
        k = 1.0
        dkappwedgedchi = self.dkappa_wedge_dchi
        omega1 = 0.0

        ap0 = self.SimData.alpha_perp
        kp0 = self.SimData.kappa_perp
        dbetawedgedchi = self.SimData.dbetawedge_dchi
        taub = self.taub

        delta1 = kp0
        delta2 = -kp0 * k**2 + 1j * dkappwedgedchi * k * (taub**2) * (
            (Te0**1.5) / n0) * (n0 * dxBz0 - Bz0 * dxn0)

        epsilon1 = -(((1j * dkappwedgedchi * k * taub**2 * Te0**4 * dxTe0) / (me * n0)))


        b1r = -((1 / (4 * n0 ** 2 * Te0 ** 4)) * ((Te0 ** 2) ** (1 / 4.0) * \
                                                  (n0 * (2 * dbetawedgedchi * taub * Te0 ** 4 * dxBz0 * \
                                                         (2 * Kx * Te0 + 3 * dxTe0) + 3 * ap0 * dc2 * n0 * \
                                                         (dxBz0 * (2 * Kx * Te0 - 5 * dxTe0) + \
                                                          2 * Te0 * d2xBz0)) + dbetawedgedchi * taub * Bz0 * \
                                                   Te0 ** 3 * (-2 * Te0 * dxn0 * (2 * Kx * Te0 + \
                                                                                  3 * dxTe0) + n0 * (
                                                                           4 * (-k ** 2 + Kx ** 2) * Te0 ** 2 + \
                                                                           3 * dxTe0 ** 2 + 6 * Te0 * (2 * Kx * dxTe0 + \
                                                                                                       d2xTe0))))))

        eta1 = 1j * ((k * dxn0) / n0) / (Te0**2.5)
        # g3 Kx^2 + g2 Kx + g1
        g1 = (ap0 * dc2 * k ** 2 * (Te0 ** 2) ** (3.0 / 4.0)) / Te0 ** 3 + dxCx0 + \
             (dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * dxn0 * \
              dxTe0) / n0 ** 2 - \
             (3 * dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * dxTe0 ** 2) / \
             (2 * n0 * Te0) - (dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * d2xTe0) / n0
        g2 = Cx0 - (dbetawedgedchi * taub *
                    (Te0**2)**(3.0 / 4.0) * dxTe0) / n0 - (3 * ap0 * dc2 * (Te0**2)**
                                                           (3.0 / 4.0) * dxTe0) / (2 * Te0**4)
        g3 = (ap0 * dc2 * (Te0**2)**(3.0 / 4.0)) / Te0**3

        p0 = (-epsilon1) * eta1 + delta2 * g1
        p1 = delta2 * g2
        p2 = (delta1 * g1 + delta2 * g3)
        p3 = delta1 * g2
        p4 = delta1 * g3

        #poly = np.vstack((p4, p3, p2, p1, p0))
        poly = np.stack((p4, p3, p2, p1, p0), axis=-1)

        return poly

    def get_uv_quartic_coefficients_corrected(self, Kx):
        '''
        g1,g2,g3, a1r,a1i,a2,b1r,b1i, b2r, b2i = self.get_quartic_coefficients( Kx)
        :param Kx:
        :return:
        '''
        n0 = self.n0
        Cx0 = self.Cx0
        Te0 = self.T0
        Bz0 = self.Bz0
        dxTe0 = self.dxT0
        dxBz0 = self.dxBz0
        dxn0 = self.dxn0
        dxCx0 = self.dxCx0
        d2xTe0 = self.dx2T0
        d2xBz0 = self.dx2Bz0

        dc2 = 0.1
        me = 9.11e-31
        k = 1.0
        dkappwedgedchi = self.dkappa_wedge_dchi
        omega1 = 0.0

        ap0 = self.SimData.alpha_perp
        kp0 = self.SimData.kappa_perp
        dbetawedgedchi = self.SimData.dbetawedge_dchi
        taub = self.taub

        delta1 = kp0
        delta2 = -kp0 * k**2 + 1j * dkappwedgedchi * k * (taub**2) * (
            (Te0**1.5) / n0) * (n0 * dxBz0 - Bz0 * dxn0)

        epsilon1 = -(((1j * dkappwedgedchi * k * taub**2 * Te0**4 * dxTe0) / (me * n0)))

        eta1 = (1 / (2. * n0 ** 2 * Te0 ** 6)) * (n0 * (2 * 1j * k * Te0 ** (7. / 2.) * dxn0 + \
               15 * ap0 * dc2 * n0 * dxBz0 * dxTe0 + 2 * dbetawedgedchi * taub * Te0 ** 4 * \
               dxBz0 * dxTe0 - 3 * ap0 * dc2 * n0 * Te0 * d2xBz0) + \
                2. * dbetawedgedchi * taub * Bz0 * Te0 ** 3 * ((-Te0) * dxn0 * dxTe0 + \
                n0 * (-2 * dxTe0 ** 2 + Te0 * (k ** 2 * Te0 + d2xTe0))))

        eta2 = (1 / (2. * n0 ** 2 * Te0 ** 5)) * (-3. * ap0 * dc2 * n0 ** 2 * dxBz0 + \
              2 * dbetawedgedchi * taub * Bz0 * Te0 ** 4 * dxn0 - 2 * dbetawedgedchi * taub * n0 * Te0 ** 3 * \
              (Te0 * dxBz0 - 2 * Bz0 * dxTe0))

        eta3 = -((dbetawedgedchi * taub * Bz0) / (n0 * Te0))
        # g3 Kx^2 + g2 Kx + g1
        g1 = (ap0 * dc2 * k ** 2 * (Te0 ** 2) ** (3.0 / 4.0)) / Te0 ** 3 + dxCx0 + \
             (dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * dxn0 * \
              dxTe0) / n0 ** 2 - \
             (3 * dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * dxTe0 ** 2) / \
             (2 * n0 * Te0) - (dbetawedgedchi * taub * (Te0 ** 2) ** (3.0 / 4.0) * d2xTe0) / n0
        g2 = Cx0 - (dbetawedgedchi * taub *
                    (Te0**2)**(3.0 / 4.0) * dxTe0) / n0 - (3 * ap0 * dc2 * (Te0**2)**
                                                           (3.0 / 4.0) * dxTe0) / (2 * Te0**4)
        g3 = (ap0 * dc2 * (Te0**2)**(3.0 / 4.0)) / Te0**3

        p0 = (-epsilon1) * eta1 + delta2 * g1
        p1 = delta2 * g2 - epsilon1 * eta2
        p2 = (delta1 * g1 + delta2 * g3 - epsilon1 * eta3)
        p3 = delta1 * g2
        p4 = delta1 * g3

        ##poly = np.vstack((p4, p3, p2, p1, p0))
        poly = np.stack((p4, p3, p2, p1, p0), axis=-1)
        return poly

    def get_roots(self, poly):
        roots = np.zeros((len(poly[:, 0]), 4), dtype=complex)
        for ii in range(len(poly[:, 0])):
            roots[ii, :] = np.roots(poly[ii, :])
        return roots


lambda_p = 5.0
sim_obj = RunInfo(scale_length=1, bz_in=50.0, lambda_p=lambda_p, pert_amp='1p')
sim_data_obj = GetSimData(sim_obj)

k = (2.0 * np.pi / lambda_p)
Kx = 1.0
wkb_obj = SimCoeffs(k, sim_data_obj)
#g1, g2, g3, a1r, a1i, a2, b1r, b1i, b2r, b2i = wkb_obj.get_quartic_coefficients(Kx)
poly_approx = wkb_obj.get_uv_quartic_coefficients(Kx)
poly = wkb_obj.get_uv_quartic_coefficients_corrected(Kx)

poly_approx[:, 1] *= 0.0
poly_approx[:, 2] *= 0.0
poly_approx[:, 3] *= 0.0

roots_approx = wkb_obj.get_roots(poly_approx)
roots_true = wkb_obj.get_roots(poly)

sim_data_obj_highBz = GetSimData(sim_obj)
mult = 10000.0
sim_data_obj_highBz.Bz = sim_data_obj_highBz.Bz * mult
sim_data_obj_highBz.wt = sim_data_obj_highBz.wt * mult
wkb_obj_highBz = SimCoeffs(k, sim_data_obj_highBz)
poly_Bz = wkb_obj_highBz.get_uv_quartic_coefficients_corrected(Kx)
roots_highBz = wkb_obj_highBz.get_roots(poly_Bz)

ixmin = 100
ixmax = 210
dTlineout = np.max(np.abs(wkb_obj.T1), axis=0)
dBlineout = np.max(np.abs(wkb_obj.Bz1), axis=0)

fig = plt.figure()
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
p0 = plot_complex(ax, poly[:, 0], c='r')
p1 = plot_complex(ax, poly[:, 1], c='g')
p2 = plot_complex(ax, poly[:, 2], c='b')
p3 = plot_complex(ax, poly[:, 3], c='m')
p4 = plot_complex(ax, poly[:, 4], c='y')

p0 = plot_complex(ax2, poly_Bz[:, 0], c='r')
p1 = plot_complex(ax2, poly_Bz[:, 1], c='g')
p2 = plot_complex(ax2, poly_Bz[:, 2], c='b')
p3 = plot_complex(ax2, poly_Bz[:, 3], c='m')
p4 = plot_complex(ax2, poly_Bz[:, 4], c='y')

ax.set_xlim(100, 207)
ax2.set_ylim(ax.get_ylim())
ax2.set_xlim(100, 207)
p_list = [p0, p1, p2, p3, p4]
leg_list = ['0', '1', '2', '3', '4']
ax.legend(p_list, leg_list)
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
x = np.arange(len(roots_approx[:, 0]))
for ii in range(4):

    ax1.scatter(x, roots_true[:, ii].real, c='b')
    ax1.scatter(x, roots_approx[:, ii].real, c='r')
    ax2.scatter(x, roots_true[:, ii].imag, c='b')
    ax2.scatter(x, roots_approx[:, ii].real, c='r')

    ax1.scatter(x, roots_highBz[:, ii].real, c='g')
    ax2.scatter(x, roots_highBz[:, ii].imag, c='g')

ax1.set_xlim(100, 207)
ax2.set_xlim(100, 207)

plt.show()
