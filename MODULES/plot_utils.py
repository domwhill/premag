'''
27/08/2019
- Rehash of plot_q_lineout.py


'''
import sys
import os
import numpy as np
import re

sys.path.extend(["./"])
import MODULES.figure_prl_twocol as fprl
from pylab import *
import MODULES.chfoil_module as cf
import MODULES.house_keeping as hk
from matplotlib import ticker
import impact_norms as inorm
import tsi_module as tsi

#---> constants...
c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31
#-----
# functions
fpre = lambda path_in: path_in.split('/')[-1]


def clear_yax(ax):
    fprl.clear_yax(ax)
    ax.set_ylabel('')
    pass


def format_yax(ax):
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%i'))
    pass


def b_lab(bz_in):
    if int(bz_in) == -1:
        lab = 'no B'
    else:
        lab = r'$%i\,\si{T}$' % (bz_in)
    return lab


def ensure_list(s):
    # Ref: https://stackoverflow.com/a/56641168/
    return s if isinstance(
        s, list) else list(s) if isinstance(s,
                                            (tuple, set, np.ndarray)) else [] if s is None else [s]


def convert_lists_to_set(a_list, b_list, c_list, d_list):
    var_list = set((a, b, c, d) for a in ensure_list(a_list) for b in ensure_list(np.sort(b_list))
                   for c in ensure_list(c_list) for d in ensure_list(d_list))
    return var_list


def absmaxND(a, axis=None):
    amax = a.max(axis)
    amin = a.min(axis)
    return np.where(-amin > amax, amin, amax)


class RunInfoList:
    """Container for information on a seris

    """

    def __init__(self, var_type, scale_length=1, bz_in=50.0, lambda_p=5, pert_amp='1p', **kwargs):
        # --> defaults
        self.lambda_p = lambda_p
        self.scale_length = scale_length
        self.bz_in = bz_in
        self.lambda_p = lambda_p
        self.pert_amp = pert_amp
        self.var_type = var_type
        self.var_list = convert_lists_to_set(self.scale_length, self.bz_in, self.lambda_p,
                                             self.pert_amp)
        #--->
        self.path_list = {}
        paths = hk.directory_paths()
        run_name = sys.argv[0].split('/')[-1].split('.')[0]
        self.save_tag = '%s%s_' % (paths.save_dir, run_name)
        self.t_max_col = 0

        self.paths = []
        self.tags = []

        self.run_objs = []
        self.paths = []
        self.run_obj_dict = {}
        for ib, var in enumerate(self.var_list):

            run_obj_loc = RunInfo(*var)
            self.run_objs.append(run_obj_loc)
            self.paths.append(run_obj_loc.path)
            tag = run_obj_loc.get_tag()
            # tagged by bz
            self.run_obj_dict['%i' % run_obj_loc.bz] = run_obj_loc
            self.tags.append(tag)
            self.save_tag = self.save_tag + 'LT' + str(scale_length) + '_'
        # find maximum simulation dump time in obj_list
        self.tmax_index = None
        self.get_max_time()
        self.save_tag = self.save_tag + str(kwargs.get("time", self.tmax_index))

    def get_max_time(self):
        """Check all IMPACT runs in paths find the latest time dump in each sim.
        return the highest time step that exists in all simulations.

        :return:
        """
        time_in_list = []
        for path in self.paths:
            t_list, tc_list = cf.get_t_list(path, var='Te')
            time_in_list.append(np.max(list(map(int, t_list))))

        tmax_index = int(np.max(time_in_list))
        self.tmax_index = '%02i' % (tmax_index)
        # time in collision times
        self.tmax_col = tc_list[tmax_index]



class RunInfo:

    def __init__(self, scale_length=1, bz_in=400.0, lambda_p=5, pert_amp='0p', dim='2D'):
        paths = hk.directory_paths()
        self.dim = dim
        self.path = paths.get_path(scale_length, bz_in, lambda_p, pert_amp=pert_amp, dim=self.dim)
        self.bz = bz_in
        self.scale_length = scale_length
        self.lambda_p = lambda_p
        self.pert_amp = pert_amp

        # --->  getting path ->
        t_list, tc_list = cf.get_t_list(self.path, var='Te')
        self.time_max = t_list[-1]
        self.tag = self.get_tag()
        self.save_tag = 'LT%i_%i_%s_%s' % (self.scale_length, self.bz, self.pert_amp, self.time_max)

        norm_dir = paths.norm_dir
        log_file = norm_dir + 'norm.log'
        [T_ref, n_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
        self.norm = cf.ConversionFactors(norm_dir, Z_ref, Ar=6.51)
        self.save_path = paths.save_dir
        self.Te_factor = 2.0 * T_ref * 1e-3

    def get_tag(self, **kwargs):
        if self.bz == '-1':
            tag = 'no B'
        else:
            tag = r'$%i \, [\si{T}]$' % self.bz
        return tag

    def get_idx_list(self):
        if int(self.scale_length) == 2:
            idx_out = [100, 110, 120]
        else:
            idx_out = [100, 120, 140]
        return idx_out


# Compute growth rate
class GetSimData:

    def __init__(self, run_obj, **kwargs):
        path = run_obj.path
        self.time = kwargs.get('time', run_obj.time_max)
        iy = kwargs.get('iy', 20)
        self.dim = run_obj.dim

        self.cfg = run_obj.norm
        if self.dim == '1D':
            dict_load = cf.load_data_all_1D(path, fpre(path), self.time)
            convert_array = lambda key: dict_load[key]
            convert_array_dy = lambda key: dict_load[key]
            self.ny = 1
        else:
            dict_load = cf.load_data_all(path, fpre(path), self.time)
            convert_array = lambda key: dict_load[key]
            convert_array_dy = lambda key: np.transpose(dict_load[key])[:, ::-1][:, iy]
            self.ny = dict_load['ny']
        # data in

        self.ne_ref = self.cfg.n0
        self.Te_ref = self.cfg.T0
        self.Bz_ref = (m_e / (q_e * self.cfg.tau_ei)) * 1.0    # assumes that norm Bz is 1T
        self.v_te = self.cfg.v_te
        self.tau_ei = self.cfg.tau_ei
        self.Z_ref = self.cfg.Z
        self.Ar = self.cfg.Ar
        self.lambda_mfp = self.cfg.lambda_mfp

        # Only load in line outs of arrays
        # swap y direction after transpose for y gradients
        self.time_secs = dict_load['time'] * self.tau_ei
        if self.dim == '1D':
            self.dyB = None
            self.dyT = None
        else:
            self.dyB = convert_array_dy('dyB')
            self.dyT = convert_array_dy('dyT')
        self.dict_load = dict_load

        self.fo = dict_load['fo']

        self.Te1 = cf.get_U_dev(dict_load['Te'].T)
        self.Bz1 = cf.get_U_dev(dict_load['Bz'].T)

        self.Z = convert_array('Z') * self.Z_ref
        self.U = convert_array('U')
        self.Te = convert_array('Te')
        self.ne = convert_array('ne')
        self.Cx = convert_array('Cx')

        self.wt = convert_array('wt')
        self.Bz = convert_array('Bz')
        if self.dim == '1D':
            self.x_grid = cf.trim_array_1D(dict_load['x_grid'], dict_load['nx'], 1)
            self.y_grid = None
        else:
            self.x_grid = cf.trim_array_1D(dict_load['x_grid'], dict_load['nx'], self.ny)
            self.y_grid = cf.trim_array_1D(dict_load['y_grid'], dict_load['ny'], dict_load['nx'])

        self.v_grid = dict_load['v_grid']

        self.dxT = convert_array('dxT')
        self.dxB = convert_array('dxB')

        self.n_factor = 1.0
        self.Te_factor = 1.0
        self.Bz_factor = 1.0
        self.C_factor = 1.0
        self.P_factor = 2.0 / 3.0
        self.x_factor = 1.0

    def conv_to_SI(self):
        self.n_factor = self.ne_ref
        self.Te_factor = (2.0 * self.Te_ref)
        self.Bz_factor = self.Bz_ref
        self.C_factor = self.v_te

        self.P_factor = 2.0 / 3.0
        self.x_factor = self.cfg.lambda_mfp
        self.U *= self.P_factor
        self.Te *= self.Te_factor
        self.ne *= self.n_factor
        self.Cx *= self.C_factor

        self.Bz *= self.Bz_factor
        self.x_grid *= self.x_factor
        self.y_grid *= self.x_factor

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
