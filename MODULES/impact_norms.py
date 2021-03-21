#----
import numpy as np

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12


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
    t0 = (2.588e-16) * (Te**(1.5)) / (Z * Z * ni * log_lambda)
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

    dict = {}
    dict['vte'] = v0
    dict['tau_ei'] = t0
    dict['nu_ei'] = nu0
    dict['lambda_mfp'] = l0
    dict['c_over_vte'] = c_by_v0
    dict['log_lambda'] = log_lambda
    dict['wpe_over_nu_ei'] = wpe_by_nu_ei
    dict['Bz_norm'] = prof_Bz_ave
    # --- get transport coeffs
    # if T in Kelvin...
    # kappa*gradT = (1/(msec))*(J/K) * (K/m) = kappa*k_b *grad (T/Kelvin)
    # W/m^2= [kappa*gradT] = (1/(m*sec))*J/m = J/(m^2 sec)

    return dict


#=======================================================================
def impact_inputs_array(ne, Te, Z, Bz, Ar):
    # Reference density [cm**-3], temperature [eV], ionisation and B-field [T]

    # Convert ne to 10**21 cm**-3
    ne = ne / 1.0e21
    ni = ne / Z

    # Calculate Coulomb np.logarithm.
    lambda_bar = (2.76e-10) / np.sqrt(Te)
    b0 = (1.44e-9) * (Z / Te)

    log_lambda = (3.0 / 2.0) * np.log(Te) - np.log(Z) - np.log(ne) / 2.0 - 1.8
    log_lambda[lambda_bar > b0] = np.log(
        Te[lambda_bar > b0]) - np.log(ne[lambda_bar > b0]) / 2.0 - 0.16 * np.ones(
            (np.shape(Te[lambda_bar > b0])))

    # Reference values
    #   v0 : thermal velocity [m/s],      t0 : collision time [s].
    #   nu0 : collision frequency [s**-1], l0 : mean free path [m])
    v0 = (0.5931e6) * (Te**0.5)
    t0 = (2.588e-16) * (Te**(1.5)) / (Z * Z * ni * log_lambda)
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

    dict = {}
    dict['vte'] = v0
    dict['tau_ei'] = t0
    dict['nu_ei'] = nu0
    dict['lambda_mfp'] = l0
    dict['c_over_vte'] = c_by_v0
    dict['log_lambda'] = log_lambda
    dict['wpe_over_nu_ei'] = wpe_by_nu_ei
    dict['Bz_norm'] = prof_Bz_ave
    # --- get transport coeffs
    # if T in Kelvin...
    # kappa*gradT = (1/(msec))*(J/K) * (K/m) = kappa*k_b *grad (T/Kelvin)
    # W/m^2= [kappa*gradT] = (1/(m*sec))*J/m = J/(m^2 sec)

    return dict
