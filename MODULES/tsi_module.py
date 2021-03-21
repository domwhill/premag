import numpy as np, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS/MODULES')
import impact_norms as ip
import EH_poly_coeff_module as ep
import matplotlib.pyplot as plt

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
cB = 3.0 * np.sqrt(np.pi) / 4.0
#---------------------------------------------------------------------------


def get_wt_grad(wt,
                Z,
                var_list=[
                    'alpha_perp', 'alpha_wedge', 'beta_perp', 'beta_wedge', 'kappa_perp',
                    'kappa_wedge'
                ]):
    '''
        Status - consider chagning this function to precompute wte gradients then coarsely sample. 
        
    '''

    dwt = 0.01
    if wt > dwt:
        wtmin = wt - dwt
        wtmax = wt + dwt
    else:
        wtmin = wt
        wtmax = wt + 2.0 * dwt
    t_dict_p1, c_dict = ep.coeff_poly_fit(wtmax, Z)
    t_dict_m1, c_dict = ep.coeff_poly_fit(wtmin, Z)
    dict_out = {}
    for var in var_list:
        dict_out[var] = (t_dict_p1[var] - t_dict_m1[var]) / (2.0 * dwt)
    return dict_out


#---------------------------------------------------------------------------
def get_dimensionful_coeff(val_dimensionless, var, ref_vals, dimensionless=True):
    n_e_m3 = ref_vals['n_e'] * 1e6
    T_eJ = q_e * ref_vals['T_e']
    tau_B = ref_vals['tau_B']
    if dimensionless == True:
        return val_dimensionless
    if var == 'kappa':
        val_out = val_dimensionless * (n_e_m3 * T_eJ * tau_B) / m_e
    elif var == 'alpha':
        val_out = val_dimensionless * (m_e * n_e_m3 / tau_B)
    elif var == 'beta':
        val_out = val_dimensionless
    return val_out


#------------------------------------------------------------------------------
def get_dimensionless_transport_c(wte, Z):
    '''
       dict,grad_dict= gget_dimensionless_transport_c(wte,Z)
    '''

    transport_dict, c_dict = ep.coeff_poly_fit(wte, Z)
    t_grad_dict = ep.coeff_poly_fit_diff(wte, Z)
    transport_dict['psi_wedge'] = transport_dict['beta_wedge']    # ettinghausen term

    return transport_dict, t_grad_dict


#------------------------------------------------------------------------------
def get_dimensionful_transport_c(wte, Z, ref_vals):
    '''
       dict,grad_dict= get_dimensionful_tranpsort_c(wte,Z)
    '''

    transport_dict, c_dict = ep.coeff_poly_fit(wte, Z)
    t_grad_dict = ep.coeff_poly_fit_diff(wte, Z)
    var_list = [
        'alpha_para', 'alpha_perp', 'alpha_wedge', 'kappa_para', 'kappa_perp', 'kappa_wedge',
        'beta_perp', 'beta_wedge'
    ]
    dict = {}
    grad_dict = {}
    for var in var_list:
        dict[var] = get_dimensionful_coeff(transport_dict[var], var.split('_')[0], ref_vals)
    for var in t_grad_dict.keys():
        grad_dict[var] = t_grad_dict[var]
    dict['psi_wedge'] = transport_dict['beta_wedge']    # ettinghausen term

    return dict, grad_dict


#-----------------------------------------------------------------------
def get_dimensionful_transport_kinetic(wte, Z, ref_vals):
    '''
       dict,grad_dict= get_dimensionful_tranpsort_c(wte,Z)
    '''

    transport_dict, c_dict = ep.coeff_poly_fit(wte, Z)
    #t_grad_dict = get_wt_grad(wte,Z)
    t_grad_dict = ep.coeff_poly_fit_diff(wte, Z)

    var_list = [
        'alpha_para', 'alpha_perp', 'alpha_wedge', 'kappa_para', 'kappa_perp', 'kappa_wedge',
        'beta_perp', 'beta_wedge'
    ]
    dict = {}
    grad_dict = {}
    for var in var_list:
        dict[var] = get_dimensionful_coeff(transport_dict[var], var.split('_')[0], ref_vals)
    for var in t_grad_dict.keys():
        grad_dict[var] = t_grad_dict[var]
    #grad_dict['kappa_wedge']  = get_dimensionful_coeff(t_grad_dict['kappa_wedge'],'kappa',ref_vals)
    #grad_dict['beta_wedge'] = t_grad_dict['beta_wedge']
    dict['psi_wedge'] = transport_dict['beta_wedge']    # ettinghausen term

    return dict, grad_dict


#-----------------------------------------------------------------------
def compute_growth(k, ne_ref, Te_ref, Z, Bz, Ar, LT, LB):
    '''
        w_p_2,w_m_2, dict_n = compute_growth(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
        
        I think this one gives the best agreement with the table in Bissell's thesis
    '''
    dict_n = ip.impact_inputs(ne_ref, Te_ref, Z, Bz,
                              Ar)    # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei = dict_n['wpe_over_nu_ei']
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']    #*cB # eB/m_e *tau_ei
    wte_in = wte * 1.0
    delta_c_norm = c_over_vte / wpe_over_nu_ei
    delta_c_norm_2 = 1e3 * Z * ((ne_ref / 1e21)**0.5) * log_lambda / (Te_ref**2)
    delta = delta_c_norm * lambda_T    # m
    dict_n['delta'] = delta
    #delta = c/wpe
    #Lambda = lambda_T_in/delta
    Lambda_1 = delta_c_norm**-1
    Lambda_2 = delta_c_norm_2**-1
    #k = (2.0*n[p.pi)/(5.0*lambda_mfp_ref)#7.0e-6#

    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T * cB

    dict, grad_dict = get_dimensionful_transport_c(wte, Z, ref_vals)
    alpha_para = dict['alpha_para']
    alpha_perp = dict['alpha_perp']
    alpha_wedge = dict['alpha_wedge']

    kappa_para = dict['kappa_para']
    dkappawedge_dchi = grad_dict['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge']
    dbetaperp_dchi = grad_dict['beta_perp']
    dalphawedge_dchi = grad_dict['alpha_wedge']

    beta_wedge = dict['beta_wedge']    # /e?
    psi_wedge = dict['beta_wedge']    # ettinghausen term
    kappa_perp = dict['kappa_perp']

    sin_theta = 1.0
    sP = 2.0 * beta_wedge * dkappawedge_dchi * sin_theta * (cB**2) * ((lambda_T**4) /
                                                                      (3.0 * LT * tau_T**2))
    sB = ((cB * wte * lambda_T**2) / (3.0 * LB * tau_T**2)) * dkappawedge_dchi * sin_theta
    sE = ((4.0 * (lambda_T**2) * delta**2) / (3.0 * tau_T**2)) * beta_wedge * psi_wedge
    dR = alpha_perp * (delta**2) / (cB * tau_T)
    dT = (cB * kappa_perp * lambda_T**2) / (3.0 * tau_T)
    pm = 1.0
    '''
    
    a1 = sB*k - (dR + dT)*1j*(k**2)
    b1 = (sB**2)*k**2 
    b2 = (sP + 2.0*sB*(dR - dT)*1j*k**3)
    b3 = - ((dR-dT)**2 + sE)*(k**4)
    coeff_sqrt =  (b1+ b2 + b3)**0.5
    pm = 1.0
    w_p_1 = 0.5*(a1) + pm*0.5*(coeff_sqrt)
    pm = -1.0
    w_m_1 = 0.5*(a1) + pm*0.5*(coeff_sqrt)
    '''
    #----- compute Bissell coefficients
    #--- check
    aN = (cB / (2.0 * wte)) * (lambda_T**2 / (tau_T)) * beta_wedge
    aE = ((2.0 * wte) / (3.0 * cB)) * (delta**2 / (tau_T)) * beta_wedge

    ckappa = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    ckappa = (cB / 3.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * dbetaperp_dchi * (lambda_T**2 / tau_T)

    aC = (wte / cB) * ((delta**2) / (tau_T)) * (1.0 + wte * dbetaperp_dchi)

    cH_hat = (wte / (cB)) * (delta**2 / tau_T) * (alpha_wedge / (wte) - dalphawedge_dchi)

    sigmaB = 1.0
    k2 = k**2
    eT = sigmaB * (ckappa / LB + aC / LB) * sin_theta
    eB = sigmaB * (ckappa / LT + aC / LT) * sin_theta
    fT = sigmaB * (cbeta / LB + (3.0 * cH_hat / (2.0 * LB))) * sin_theta
    fB = sigmaB * (cbeta / LT + 3.0 * cH_hat / (2.0 * LT)) * sin_theta
    vB = (sigmaB / LB) * (ckappa + aC) * sin_theta
    s = 4.0 * (aN * 1j * k2 + fT * k) * (aE * 1j * k2 + eB * k)
    a = (eT + fB) * k - (dT + dR) * 1j * k2
    b = (((dT - dR) * 1j * k2 - (eT - fB) * k)**2 + s)**0.5
    w_p_2 = 0.5 * (a + b)
    w_m_2 = 0.5 * (a - b)
    return w_p_2, w_m_2, dict_n


#-----------------------------------------------------------------------
def compute_growth_array(k, ne_ref, Te_ref, Z, Bz, Ar, LT, LB):
    '''
        w_p_2,w_m_2, dict_n = compute_growth(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
        
        I think this one gives the best agreement with the table in Bissell's thesis
    '''
    dict_n = ip.impact_inputs_array(
        ne_ref, Te_ref, Z, Bz, Ar)    # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei = dict_n['wpe_over_nu_ei']
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']    #*cB # eB/m_e *tau_ei
    wte_in = wte * 1.0
    delta_c_norm = c_over_vte / wpe_over_nu_ei
    delta_c_norm_2 = 1e3 * Z * ((ne_ref / 1e21)**0.5) * log_lambda / (Te_ref**2)
    delta = delta_c_norm * lambda_T    # m
    dict_n['delta'] = delta
    #delta = c/wpe
    #Lambda = lambda_T_in/delta
    Lambda_1 = delta_c_norm**-1
    Lambda_2 = delta_c_norm_2**-1
    #k = (2.0*n[p.pi)/(5.0*lambda_mfp_ref)#7.0e-6#

    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T * cB

    dict, grad_dict = get_dimensionful_transport_c(wte, Z, ref_vals)

    alpha_para = dict['alpha_para']
    alpha_perp = dict['alpha_perp']
    alpha_wedge = dict['alpha_wedge']

    kappa_para = dict['kappa_para']
    dkappawedge_dchi = grad_dict['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge']
    dbetaperp_dchi = grad_dict['beta_perp']
    dalphawedge_dchi = grad_dict['alpha_wedge']

    beta_wedge = dict['beta_wedge']    # /e?
    psi_wedge = dict['beta_wedge']    # ettinghausen term
    kappa_perp = dict['kappa_perp']

    #----- compute Bissell coefficients
    #--- check
    #<---
    sin_theta = 1.0
    sP = 2.0 * beta_wedge * dkappawedge_dchi * sin_theta * (cB**2) * ((lambda_T**4) /
                                                                      (3.0 * LT * tau_T**2))
    sB = ((cB * wte * lambda_T**2) / (3.0 * LB * tau_T**2)) * dkappawedge_dchi * sin_theta
    sE = ((4.0 * (lambda_T**2) * delta**2) / (3.0 * tau_T**2)) * beta_wedge * psi_wedge
    dR = alpha_perp * (delta**2) / (cB * tau_T)
    dT = (cB * kappa_perp * lambda_T**2) / (3.0 * tau_T)
    #--->

    aN = (cB / (2.0 * wte)) * (lambda_T**2 / (tau_T)) * beta_wedge
    aE = ((2.0 * wte) / (3.0 * cB)) * (delta**2 / (tau_T)) * beta_wedge

    ckappa = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    ckappa = (cB / 3.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * dbetaperp_dchi * (lambda_T**2 / tau_T)

    aC = (wte / cB) * ((delta**2) / (tau_T)) * (1.0 + wte * dbetaperp_dchi)

    cH_hat = (wte / (cB)) * (delta**2 / tau_T) * (alpha_wedge / (wte) - dalphawedge_dchi)

    sigmaB = -1
    k2 = k**2
    eT = sigmaB * (ckappa / LB + aC / LB) * sin_theta
    eB = sigmaB * (ckappa / LT + aC / LT) * sin_theta
    fT = sigmaB * (cbeta / LB + (3.0 * cH_hat / (2.0 * LB))) * sin_theta
    fB = sigmaB * (cbeta / LT + 3.0 * cH_hat / (2.0 * LT)) * sin_theta
    vB = (sigmaB / LB) * (ckappa + aC) * sin_theta
    s = 4.0 * (aN * 1j * k2 + fT * k) * (aE * 1j * k2 + eB * k)
    a = (eT + fB) * k - (dT + dR) * 1j * k2
    b = (((dT - dR) * 1j * k2 - (eT - fB) * k)**2 + s)**0.5
    w_p_2 = 0.5 * (a + b)
    w_m_2 = 0.5 * (a - b)
    return w_p_2, w_m_2, dict_n


#-----------------------------------------------------------------------
def compute_growth_array_kinetic(k, ne_ref, Te_ref, Z, Bz, Ar, LT, LB, dict_ratios):
    '''
        w_p_2,w_m_2, dict_n = compute_growth(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
        
        Insert in a dictionary containing the ratios of the kinetic vs. classical transport coeffs. eg.
        dict['kappa_wedge'] = q_{RL,k,x}/q_{RL,c,x}
        
        # 
        
        I think this one gives the best agreement with the table in Bissell's thesis
    '''
    dict_n = ip.impact_inputs_array(
        ne_ref, Te_ref, Z, Bz, Ar)    # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei = dict_n['wpe_over_nu_ei']
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']    #*cB # eB/m_e *tau_ei
    wte_in = wte * 1.0
    delta_c_norm = c_over_vte / wpe_over_nu_ei
    delta_c_norm_2 = 1e3 * Z * ((ne_ref / 1e21)**0.5) * log_lambda / (Te_ref**2)
    delta = delta_c_norm * lambda_T    # m
    dict_n['delta'] = delta
    #delta = c/wpe
    #Lambda = lambda_T_in/delta
    Lambda_1 = delta_c_norm**-1
    Lambda_2 = delta_c_norm_2**-1
    #k = (2.0*n[p.pi)/(5.0*lambda_mfp_ref)#7.0e-6#

    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T * cB

    dict, grad_dict = get_dimensionful_transport_c(wte, Z, ref_vals)

    alpha_para = dict['alpha_para'] * dict_ratio['alpha_para']
    alpha_perp = dict['alpha_perp'] * dict_ratio['alpha_perp']
    alpha_wedge = dict['alpha_wedge'] * dict_ratio['alpha_wedge']

    kappa_para = dict['kappa_para'] * dict_ratio['kappa_para']
    dkappawedge_dchi = grad_dict['kappa_wedge'] * dict_ratio['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge'] * dict_ratio['beta_wedge']
    dbetaperp_dchi = grad_dict['beta_perp'] * dict_ratio['beta_perp']
    dalphawedge_dchi = grad_dict['alpha_wedge'] * dict_ratio['alpha_wedge']

    beta_wedge = dict['beta_wedge'] * dict_ratio['beta_wedge']    # /e?
    psi_wedge = dict['beta_wedge'] * dict_ratio['beta_wedge']    # ettinghausen term
    kappa_perp = dict['kappa_perp'] * dict_ratio['kappa_perp']

    #----- compute Bissell coefficients
    #--- check
    #<---
    sin_theta = 1.0
    sP = 2.0 * beta_wedge * dkappawedge_dchi * sin_theta * (cB**2) * ((lambda_T**4) /
                                                                      (3.0 * LT * tau_T**2))
    sB = ((cB * wte * lambda_T**2) / (3.0 * LB * tau_T**2)) * dkappawedge_dchi * sin_theta
    sE = ((4.0 * (lambda_T**2) * delta**2) / (3.0 * tau_T**2)) * beta_wedge * psi_wedge
    dR = alpha_perp * (delta**2) / (cB * tau_T)
    dT = (cB * kappa_perp * lambda_T**2) / (3.0 * tau_T)
    #--->

    aN = (cB / (2.0 * wte)) * (lambda_T**2 / (tau_T)) * beta_wedge
    aE = ((2.0 * wte) / (3.0 * cB)) * (delta**2 / (tau_T)) * beta_wedge

    ckappa = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    ckappa = (cB / 3.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * dbetaperp_dchi * (lambda_T**2 / tau_T)

    aC = (wte / cB) * ((delta**2) / (tau_T)) * (1.0 + wte * dbetaperp_dchi)

    cH_hat = (wte / (cB)) * (delta**2 / tau_T) * (alpha_wedge / (wte) - dalphawedge_dchi)

    sigmaB = 1
    k2 = k**2
    eT = sigmaB * (ckappa / LB + aC / LB) * sin_theta
    eB = sigmaB * (ckappa / LT + aC / LT) * sin_theta
    fT = sigmaB * (cbeta / LB + (3.0 * cH_hat / (2.0 * LB))) * sin_theta
    fB = sigmaB * (cbeta / LT + 3.0 * cH_hat / (2.0 * LT)) * sin_theta
    vB = (sigmaB / LB) * (ckappa + aC) * sin_theta
    s = 4.0 * (aN * 1j * k2 + fT * k) * (aE * 1j * k2 + eB * k)
    a = (eT + fB) * k - (dT + dR) * 1j * k2
    b = (((dT - dR) * 1j * k2 - (eT - fB) * k)**2 + s)**0.5
    w_p_2 = 0.5 * (a + b)
    w_m_2 = 0.5 * (a - b)
    return w_p_2, w_m_2, dict_n


#-----------------------------------------------------------------------
def compute_growth_k(k, ne_ref, Te_ref, Z, Bz, Ar, LT, LB, cfg):
    '''
        w_p_2,w_m_2 = compute_growth(k,ne_ref,Te_ref,Z,Bz,6.5)
        
        I think this one gives the best agreement with the table in Bissell's thesis
    '''
    dict_n = ip.impact_inputs(ne_ref, Te_ref, Z, Bz,
                              Ar)    # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei = dict_n['wpe_over_nu_ei']
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']    #*cB # eB/m_e *tau_ei
    wte_in = wte * 1.0
    delta_c_norm = c_over_vte / wpe_over_nu_ei
    delta_c_norm_2 = 1e3 * Z * ((ne_ref / 1e21)**0.5) * log_lambda / (Te_ref**2)
    delta = delta_c_norm * lambda_T    # m
    dict_n['delta'] = delta
    #delta = c/wpe
    #Lambda = lambda_T_in/delta
    Lambda_1 = delta_c_norm**-1
    Lambda_2 = delta_c_norm_2**-1
    #k = (2.0*n[p.pi)/(5.0*lambda_mfp_ref)#7.0e-6#

    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T * cB

    #dict,grad_dict= get_dimensionful_transport_c(wte,Z,ref_vals)

    alpha_para = dict['alpha_para']
    alpha_perp = dict['alpha_perp']
    alpha_wedge = dict['alpha_wedge']

    kappa_para = dict['kappa_para']
    dkappawedge_dchi = grad_dict['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge']
    dbetaperp_dchi = grad_dict['beta_perp']
    dalphawedge_dchi = grad_dict['alpha_wedge']

    beta_wedge = dict['beta_wedge']    # /e?
    psi_wedge = dict['beta_wedge']    # ettinghausen term
    kappa_perp = dict['kappa_perp']

    sin_theta = 1.0
    sP = 2.0 * beta_wedge * dkappawedge_dchi * sin_theta * (cB**2) * ((lambda_T**4) /
                                                                      (3.0 * LT * tau_T**2))
    sB = ((cB * wte * lambda_T**2) / (3.0 * LB * tau_T**2)) * dkappawedge_dchi * sin_theta
    sE = ((4.0 * (lambda_T**2) * delta**2) / (3.0 * tau_T**2)) * beta_wedge * psi_wedge
    dR = alpha_perp * (delta**2) / (cB * tau_T)
    dT = (cB * kappa_perp * lambda_T**2) / (3.0 * tau_T)
    pm = 1.0

    a1 = sB * k - (dR + dT) * 1j * (k**2)
    b1 = (sB**2) * k**2
    b2 = (sP + 2.0 * sB * (dR - dT) * 1j * k**3)
    b3 = -((dR - dT)**2 + sE) * (k**4)
    coeff_sqrt = (b1 + b2 + b3)**0.5
    pm = 1.0
    w_p_1 = 0.5 * (a1) + pm * 0.5 * (coeff_sqrt)
    pm = -1.0
    w_m_1 = 0.5 * (a1) + pm * 0.5 * (coeff_sqrt)

    #----- compute Bissell coefficients
    #--- check
    aN = (cB / (2.0 * wte)) * (lambda_T**2 / (tau_T)) * beta_wedge
    aE = ((2.0 * wte) / (3.0 * cB)) * (delta**2 / (tau_T)) * beta_wedge

    ckappa = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    ckappa = (cB / 3.0) * wte * dkappawedge_dchi * (lambda_T**2 / tau_T)
    cbeta = (cB / 2.0) * dbetaperp_dchi * (lambda_T**2 / tau_T)

    aC = (wte / cB) * ((delta**2) / (tau_T)) * (1.0 + wte * dbetaperp_dchi)

    cH_hat = (wte / (cB)) * (delta**2 / tau_T) * (alpha_wedge / (wte) - dalphawedge_dchi)

    sigmaB = 1
    k2 = k**2
    eT = sigmaB * (ckappa / LB + aC / LB) * sin_theta
    eB = sigmaB * (ckappa / LT + aC / LT) * sin_theta
    fT = sigmaB * (cbeta / LB + (3.0 * cH_hat / (2.0 * LB))) * sin_theta
    fB = sigmaB * (cbeta / LT + 3.0 * cH_hat / (2.0 * LT)) * sin_theta
    vB = (sigmaB / LB) * (ckappa + aC) * sin_theta
    s = 4.0 * (aN * 1j * k2 + fT * k) * (aE * 1j * k2 + eB * k)
    a = (eT + fB) * k - (dT + dR) * 1j * k2
    b = (((dT - dR) * 1j * k2 - (eT - fB) * k)**2 + s)**0.5
    w_p_2 = 0.5 * (a + b)
    w_m_2 = 0.5 * (a - b)
    return w_p_1, w_m_1, w_p_2, w_m_2, dict_n


#-----------------------------------------------------------------------
