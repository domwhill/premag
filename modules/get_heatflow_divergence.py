'''
    plots the divergence of q_RL
'''
import numpy as np, re, os, sys, matplotlib.pyplot as plt

sys.path.extend(["./"])
import modules.chfoil_module as cf
import modules.figure_latex as fprl
import modules.kinetic_ohmslaw_module_varZ as kohb
import modules.house_keeping as hk
from matplotlib import ticker

norm_path = hk.DataDirectoryPaths().norm_dir
cfg = cf.ConversionFactors(norm_path)



def repack_2D(path, time, cfg=cfg):
    ''' Repacks computed transport terms into a dictionary.

        dict_c,dict_k = repack_2D(path,time)
    '''
    dict_qc, dict_qk = kohb.get_q_abs(path, time)
    v_nx_k, v_ny_k, v_nx_c, v_ny_c = kohb.get_Nernst_abs(path, time)

    dict_out_c = {}
    dict_out_k = {}
    var_list = ['SH x', 'SH y', 'RL x', 'RL y', 'E x', 'E y', 'vN x', 'vN y']

    for var in var_list:
        if var[0] != 'v':
            dict_out_k[var] = dict_qk[var]
        elif var == 'vN x':
            dict_out_k[var] = v_nx_k
        elif var == 'vN y':
            dict_out_k[var] = v_ny_k

    for var in var_list:
        if var[0] != 'v':
            dict_out_c[var] = dict_qc[var]
        elif var == 'vN x':
            dict_out_c[var] = v_nx_c
        elif var == 'vN y':
            dict_out_c[var] = v_ny_c

    dict_out_k['U'] = dict_qk['U']
    dict_out_k['Te'] = dict_qk['Te']
    dict_out_k['Bz'] = dict_qk['Bz']
    dict_out_k['x_grid'] = dict_qk['x_grid']
    return dict_out_c, dict_out_k


def load_qdata(path, time='10', cfg=cfg):
    """Load heat flow terms from simulation located at path."""
    kohnew = {}

    tt = int(time)
    string_tt_glb = '%2.2i' % int(time)
    time = string_tt_glb
    print '\n\nkinetic model time = ', time
    dict = kohb.get_kinetic_heatflow_b(path, str(time))
    kohnew[path] = dict

    return kohnew


def get_divqRL(path, time, cfg=cfg, **kwargs):
    '''
        dyqy, qlab = get_divqRL(path,time,cfg=cfg,**kwargs)
    :param path:
    :param time:
    :param cfg:
    :param kwargs:
    :return:
    '''
    switch_kinetic_on = kwargs.pop('switch_kinetic_on', True)

    kohnew = load_qdata(path, time)
    dict_c, dict_k = repack_2D(path, time)

    x_grid_SI = kohnew[path]['x_grid'] * cfg.xstep_factor
    y_grid_SI = kohnew[path]['y_grid'] * cfg.xstep_factor
    # --- do everything in units of 10 eV/ps

    if switch_kinetic_on:
        data_yc = dict_k['RL y']
        kc_str = 'k'
    else:
        data_yc = dict_c['RL y']
        kc_str = 'c'
    data_y = data_yc * 1.0

    dxqy, dyqy = cf.get_grad(x_grid_SI, y_grid_SI, data_y)
    dyqy = np.transpose(dyqy)    # q_SH_y[path][tt,:,:]
    MULT = -1.0

    dyqy *= MULT * cfg.divq_factor * 0.1
    qlab = r'-$\partial_y q_{y,RL,%s}$ %s' % (kc_str, cfg.divq_unit)

    return dyqy, qlab
