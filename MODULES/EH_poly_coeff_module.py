import numpy as np
import os
from scipy import interpolate as SI
import getpass
q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27

userid = getpass.getuser()


class MyDict(dict):
    pass


def set_vars():
    dict = {}
    Z = 10.0
    Ar = 2.0 * Z
    dict['Z'] = Z
    dict['Ar'] = Ar
    dict['mi'] = Ar * (1.66e-27)
    dict['flux_lim'] = 1.0
    dict['kappa'] = 1e-33    # kappa/mi^-3.5
    dict['Mach_crit'] = 0.95
    dict['Y'] = 0.0    # energy flow from the ablation profile into the solid
    dict['v_crit'] = 1.0
    dict['rho_crit'] = 1.0

    #T0_over_mi = 1.0e13 # cm^2/sec^2
    #rho0 = 0.3 #g/cm^3
    #v0 = 2.2e5#cm/s
    #g = 0.0#3.0e15 # cm/sec^2
    #T0 = T0_over_mi

    #xmax = (15e-6)*100.0 # 15 mu m
    #x0 = 0.0
    #xmin = -0.000866666666667 #1e-5#-(15.0e-6)*100.0
    #xmax = (15e-6)*100.0 # 15 mu m

    return dict


#  Takes plasma reference values and converts to useful IMPACT values.
def impact_inputs(ne, Te, Z, Bz, Ar):
    # Reference density [cm**-3], temperature [eV], ionisation and B-field [T]
    #ne = 1.0e+20#10.0e19#0.5*9.08e21
    #Te = 30.0
    #Z = 2.0 #
    #Bz = 3.75 #3.75
    #Ar = 40.0

    # Convert ne to 10**21 cm**-3
    ne = ne / 1.0e21
    ni = ne / Z

    # Calculate Coulomb np.logarithm.
    lambda_bar = (2.76e-10) / np.sqrt(Te)
    b0 = (1.44e-9) * (Z / Te)

    if (lambda_bar > b0):
        log_lambda = np.log(Te) - np.log(ne) / 2.0 - 0.16
    else:
        log_lambda = (3.0 / 2.0) * np.log(Te) - np.log(Z) - np.log(ne) / 2.0 - 1.8

    # Reference values
    #   v0 : thermal velocity [m/s],      t0 : collision time [s].
    #   nu0 : collision frequency [s**-1], l0 : mean free path [m])
    v0 = (0.5931e6) * (Te**0.5)
    t0 = (2.588e-16) * (Te**(3.0 / 2.0)) / (Z * Z * ni * log_lambda)
    nu0 = 1.0 / (t0)
    l0 = (1.535e-10) * (Te * Te) / (Z * Z * ni * log_lambda)

    # IMPACT inputs
    wpe_by_nu_ei = 0.4618 * (Te**(1.5)) / (Z * np.sqrt(ne) * log_lambda)
    # for tau_B instead of tau_ei as in IMPACTA - note that ln_lam will be
    # different too
    #wpe_by_nu_ei = 0.4618 * (Te**(3/2)) / (Z * sqrt(ne) * log_lambda) * (3*sqrt(pi)/4)
    c_by_v0 = 1.0 / ((1.9784e-3) * np.sqrt(Te))
    prof_Bz_ave = (1.756e11) * Bz * t0

    # Display
    print '\nINPUT QUANTITIES'
    print 'Density \t', ne * 1e21, '[cm**3]'

    print 'Temperature \t', Te, '[eV]'
    print 'Ionisation \t', Z, '[ ]'
    print 'Bz \t\t', Bz, '[T]'

    print '\n IMPACT VARIABLES'
    print 'log_lambda \t', log_lambda

    print 'wpe_by_nu_ei \t', wpe_by_nu_ei

    print 'c / v0 \t\t', c_by_v0
    print 'prof_Bz_ave \t', prof_Bz_ave

    #clear arry
    print '\n\nPLASMA REFERENCE VARIABLES'
    print 'Reference thermal velocity \t %1.5e [m/s]' % (v0)
    print 'Reference collision time \t %1.5e [s]' % (t0)
    print 'Reference collision frequency \t%1.5e [s**-1]' % (nu0)
    print 'Reference mfp \t\t \t %1.5e [m]\n' % (l0)
    dict = MyDict()
    dict['vte'] = v0
    dict['tau_ei'] = t0
    dict['nu_ei'] = nu0
    dict['lambda_mfp'] = l0
    dict['c_over_vte'] = c_by_v0
    dict['log_lambda'] = log_lambda
    dict['wpe_over_nu_ei'] = wpe_by_nu_ei
    dict['Bz_norm'] = prof_Bz_ave
    # --- get transport coeffs
    kappa_c = 1.0    # dimensionless transport coeff
    ne_m3 = (ne * 1.0e21) / (1.0e6)
    T_Joules = q_e * Te
    kappa = (ne_m3 * T_Joules * t0 / m_e) * kappa_c    # m^-1 s^-1
    # if T in Kelvin...
    # kappa*gradT = (1/(msec))*(J/K) * (K/m) = kappa*k_b *grad (T/Kelvin)
    # W/m^2= [kappa*gradT] = (1/(m*sec))*J/m = J/(m^2 sec)

    return dict


def load_EHcoeffs():
    dirname, filename = os.path.split(os.path.abspath(__file__))
    fname = '%s/EH_coeffs.txt' % dirname
    data = np.loadtxt(fname, dtype=str)
    ylen = len(data[:, 0])
    xlen = len(data[0, :])
    f_data = np.zeros((xlen - 1), dtype=float)
    dict = {}

    for nn in range(ylen):
        f_data = np.zeros((xlen - 1), dtype=float)
        f_data[:] = data[nn, 1:]
        dict[data[nn, 0]] = f_data[:]

    return dict


def interp_data(x_in, y_in, x_data_smooth):
    #'p'
    #y_out = interp_data(x_in,y_in,x_out)
    #'p'
    f = SI.PchipInterpolator(x_in, y_in, extrapolate=True)
    y_data_smooth = f(x_data_smooth)    #linear

    return y_data_smooth


def get_coeff(coeff_name, Z):
    dict = load_EHcoeffs()
    Z_arr = dict['Z']
    Z_arr[-1] = 1000.0
    coeff_val = interp_data(Z_arr, dict[coeff_name], Z)
    return coeff_val


def get_coeff_dict(Z):
    dict = load_EHcoeffs()
    Z_arr = dict['Z']
    Z_arr[-1] = 1000.0
    out_dict = {}
    for name in dict.keys():
        coeff_val = interp_data(Z_arr, dict[name], Z)
        out_dict[name] = coeff_val
    return out_dict


def coeff_poly_fit(chi, Z):
    '''
        t_dict,coeff_dict = coeff_poly_fit(chi,Z)
    '''
    alpha_0 = get_coeff("alpha_0", Z)
    alpha_0_p = get_coeff("alpha_0'", Z)
    alpha_1_p = get_coeff("alpha_0'", Z)
    a_0_p = get_coeff("a_0'", Z)
    a_1_p = get_coeff("a_1'", Z)
    alpha_0_wedge = get_coeff("alpha_0_wedge", Z)
    alpha_0_pp = get_coeff("alpha_0''", Z)
    alpha_1_pp = get_coeff("alpha_1''", Z)
    a_0_pp = get_coeff("a_0''", Z)
    a_1_pp = get_coeff("a_1''", Z)
    a_2_pp = get_coeff("a_2''", Z)
    beta_0 = get_coeff("beta_0", Z)
    beta_0_p = get_coeff("beta_0'", Z)
    beta_1_p = get_coeff("beta_1'", Z)
    b_0_p = get_coeff("b_0'", Z)
    b_1_p = get_coeff("b_1'", Z)
    b_2_p = get_coeff("b_2'", Z)
    beta_0_wedge = get_coeff("beta_0_wedge", Z)
    beta_0_pp = get_coeff("beta_0''", Z)
    beta_1_pp = get_coeff("beta_1''", Z)
    b_0_pp = get_coeff("b_0''", Z)
    b_1_pp = get_coeff("b_1''", Z)
    b_2_pp = get_coeff("b_2''", Z)
    gamma_0 = get_coeff("gamma_0", Z)
    gamma_0_p = get_coeff("gamma_0'", Z)
    gamma_1_p = get_coeff("gamma_1'", Z)
    c_0_p = get_coeff("c_0'", Z)
    c_1_p = get_coeff("c_1'", Z)
    c_2_p = get_coeff("c_2'", Z)
    gamma_0_wedge = get_coeff("gamma_0_wedge", Z)
    gamma_0_pp = get_coeff("gamma_0''", Z)
    gamma_1_pp = get_coeff("gamma_1''", Z)
    c_0_pp = get_coeff("c_0''", Z)
    c_1_pp = get_coeff("c_1''", Z)
    c_2_pp = get_coeff("c_2''", Z)

    alpha_para = 1.0 - (alpha_0_p / a_0_p)
    alpha_perp = 1.0 - (alpha_1_p * chi + alpha_0_p) / (chi**2 + a_1_p * chi + a_0_p)
    alpha_wedge = chi * (alpha_1_pp * chi + alpha_0_pp) / ((
        (chi**3) + (a_2_pp * (chi**2)) + a_1_pp * chi + a_0_pp)**(8.0 / 9.0))
    beta_para = (beta_0_p / (b_0_p**(8.0 / 9.0)))
    beta_perp = (beta_1_p * chi + beta_0_p) / (((chi**3) +
                                                (b_2_p *
                                                 (chi**2)) + b_1_p * chi + b_0_p)**(8.0 / 9.0))
    beta_wedge = chi * (beta_1_pp * chi + beta_0_pp) / (chi**3 + b_2_pp *
                                                        (chi**2) + b_1_pp * chi + b_0_pp)

    kappa_para = gamma_0
    kappa_perp = (gamma_1_p * chi + gamma_0_p) / ((chi**3) + (c_2_p *
                                                              (chi**2)) + c_1_p * chi + c_0_p)
    kappa_wedge = chi * ((gamma_1_pp * chi) + gamma_0_pp) / ((chi**3) + c_2_pp *
                                                             (chi**2) + c_1_pp * chi + c_0_pp)

    t_dict = {}
    t_dict['alpha_para'] = alpha_para
    t_dict['alpha_perp'] = alpha_perp
    t_dict['alpha_wedge'] = alpha_wedge
    t_dict['beta_para'] = beta_para
    t_dict['beta_perp'] = beta_perp
    t_dict['beta_wedge'] = beta_wedge
    t_dict['kappa_para'] = kappa_para
    t_dict['kappa_perp'] = kappa_perp
    t_dict['kappa_wedge'] = kappa_wedge

    coeff_dict = get_coeff_dict(Z)

    return t_dict, coeff_dict


def coeff_poly_fit_diff(chi, Z):
    '''
        Finds the differential forms of the EH transport coefficients.
        
        Uses the analytic forms ( derived in Mathematica).
        t_dict = coeff_poly_fit_diff(chi,Z)
    '''
    alpha0 = get_coeff("alpha_0", Z)
    alpha0p = get_coeff("alpha_0'", Z)
    alpha1p = get_coeff("alpha_0'", Z)
    a0p = get_coeff("a_0'", Z)
    a1p = get_coeff("a_1'", Z)
    alpha0wedge = get_coeff("alpha_0_wedge", Z)
    alpha0pp = get_coeff("alpha_0''", Z)
    alpha1pp = get_coeff("alpha_1''", Z)
    a0pp = get_coeff("a_0''", Z)
    a1pp = get_coeff("a_1''", Z)
    a2pp = get_coeff("a_2''", Z)
    beta0 = get_coeff("beta_0", Z)
    beta0p = get_coeff("beta_0'", Z)
    beta1p = get_coeff("beta_1'", Z)
    b0p = get_coeff("b_0'", Z)
    b1p = get_coeff("b_1'", Z)
    b2p = get_coeff("b_2'", Z)
    beta0wedge = get_coeff("beta_0_wedge", Z)
    beta0pp = get_coeff("beta_0''", Z)
    beta1pp = get_coeff("beta_1''", Z)
    b0pp = get_coeff("b_0''", Z)
    b1pp = get_coeff("b_1''", Z)
    b2pp = get_coeff("b_2''", Z)
    gamma0 = get_coeff("gamma_0", Z)
    gamma0p = get_coeff("gamma_0'", Z)
    gamma1p = get_coeff("gamma_1'", Z)
    c0p = get_coeff("c_0'", Z)
    c1p = get_coeff("c_1'", Z)
    c2p = get_coeff("c_2'", Z)
    gamma0wedge = get_coeff("gamma_0_wedge", Z)
    gamma0pp = get_coeff("gamma_0''", Z)
    gamma1pp = get_coeff("gamma_1''", Z)
    c0pp = get_coeff("c_0''", Z)
    c1pp = get_coeff("c_1''", Z)
    c2pp = get_coeff("c_2''", Z)

    dalphaperp = (a1p * alpha0p - a0p * alpha1p + chi *
                  (2 * alpha0p + alpha1p * chi)) / (a0p + chi * (a1p + chi))**2

    dalphawedge =  -((0.8888888888888888*chi*(alpha0pp + alpha1pp*chi)*(a1pp + chi*(2*a2pp + 3*chi)))/(a0pp + chi*(a1pp + chi*(a2pp + chi)))**
     1.8888888888888888) + (alpha1pp*chi)/ \
    (a0pp + chi*(a1pp + chi*(a2pp + chi)))**0.8888888888888888 + \
    (alpha0pp + alpha1pp*chi)/(a0pp + chi*(a1pp + chi*(a2pp + chi)))**0.8888888888888888

    dbetaperp = -((0.8888888888888888*(beta0p + beta1p*chi)*(b1p + chi*(2*b2p + 3*chi)))/\
    (b0p + chi*(b1p + chi*(b2p + chi)))**1.8888888888888888) + \
    beta1p/(b0p + chi*(b1p + chi*(b2p + chi)))**0.8888888888888888

    dbetawedge = (b0pp*(beta0pp + 2*beta1pp*chi) - chi**2*(b2pp*beta0pp - b1pp*beta1pp + \
     2*beta0pp*chi + beta1pp*chi**2))/(b0pp + chi*(b1pp + chi*(b2pp + chi)))**2

    dkappaperp = ((-c1p)*gamma0p + c0p*gamma1p - chi*((2*c2p + 3*chi)*gamma0p + \
     chi*(c2p + 2*chi)*gamma1p))/(c0p + chi*(c1p + chi*(c2p + chi)))**2

    dkappawedge = ((c0pp - chi**2*(c2pp + 2*chi))*gamma0pp + chi*(2*c0pp + c1pp*chi - chi**3)*gamma1pp) \
    /(c0pp + chi*(c1pp + chi*(c2pp + chi)))**2

    t_dict = {}
    t_dict['alpha_para'] = 0.0
    t_dict['alpha_perp'] = dalphaperp
    t_dict['alpha_wedge'] = dalphawedge

    t_dict['beta_para'] = 0.0
    t_dict['beta_perp'] = dbetaperp
    t_dict['beta_wedge'] = dbetawedge

    t_dict['kappa_para'] = 0.0
    t_dict['kappa_perp'] = dkappaperp
    t_dict['kappa_wedge'] = dkappawedge

    return t_dict


if __name__ == "__main__":

    dict = load_EHcoeffs()
    print dict.keys()
    print '\n\n dict: \n: ', dict
    print '\n\n Z: ', dict['Z']

    print '\n\n alpha_0: ', dict['alpha_0']
    print '\n\n alpha_0_prime: ', dict["alpha_0'"]
    Z = 13.0
    chi = 0.1    #omega*tau
    t_dict, coeff_dict = coeff_poly_fit(chi, Z)
    print t_dict
    #alpha_0'= coeff'oly_fit(Z)
