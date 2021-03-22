import numpy as np
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
    #print 'dict: ', dict.keys()
    return dict


#  Takes plasma reference values and converts to useful IMPACT values.


def load_EHcoeffs():
    fname = '/Users/dominichill/Dropbox/PhD/Code/AblationMaths/EH_coeffs.txt'
    #fname = '~/Dropbox/PhD/Code/AblationMaths/EH_coeffs.txt'
    #fname = '/Users/DH814/Dropbox/PhD/Code/AblationMaths/EH_coeffs.txt'
    fname = '/Users/' + userid + '/Dropbox/PhD/Code/AblationMaths/EH_coeffs_py3.txt'
    # fname = '/Users/dominichill/Dropbox/PhD/Code/AblationMaths/EH_coeffs_py3.txt'
    f = open(fname, 'r')
    data = f.readlines()
    #------> modulate -->
    nd = len(data[0].split())

    #data_out = np.zeros((len(data)-1,nd-1),dtype=float)
    i = 0
    #lab_list = []
    dict = {}
    for il in range(0, len(data)):
        var_in = data[il].split()[0]
        data_in = data[il].split()
        #print(data_in)
        #print(len(data_in[1:]))
        #lab_list.append(var_in)
        dict[var_in] = data_in[1:]
        #data_out[il-1,:] = data_in[1:]
    #<-----------------
    '''
    print(' data = ', data)
    ylen = len(data[:,0])
    xlen = len(data[0,:])
    f_data = np.zeros((xlen-1),dtype=float)
    dict = {}

    for nn in range(ylen):
        f_data = np.zeros((xlen-1),dtype=float)
        print('data = ', data[nn,1:])
        f_data[:] = data[nn,1:]
        dict[data[nn,0]] = f_data[:]
        #print data[nn,0],'f_data: ', f_data
        #print 'dict[data[nn,0]: ', dict[data[nn,0]]
    '''
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
    alpha_0_p = get_coeff("alpha_0p", Z)
    alpha_1_p = get_coeff("alpha_0p", Z)
    a_0_p = get_coeff("a_0p", Z)
    a_1_p = get_coeff("a_1p", Z)
    alpha_0_wedge = get_coeff("alpha_0_wedge", Z)
    alpha_0_pp = get_coeff("alpha_0pp", Z)
    alpha_1_pp = get_coeff("alpha_1pp", Z)
    a_0_pp = get_coeff("a_0pp", Z)
    a_1_pp = get_coeff("a_1pp", Z)
    a_2_pp = get_coeff("a_2pp", Z)
    beta_0 = get_coeff("beta_0", Z)
    beta_0_p = get_coeff("beta_0p", Z)
    beta_1_p = get_coeff("beta_1p", Z)
    b_0_p = get_coeff("b_0p", Z)
    b_1_p = get_coeff("b_1p", Z)
    b_2_p = get_coeff("b_2p", Z)
    beta_0_wedge = get_coeff("beta_0_wedge", Z)
    beta_0_pp = get_coeff("beta_0pp", Z)
    beta_1_pp = get_coeff("beta_1pp", Z)
    b_0_pp = get_coeff("b_0pp", Z)
    b_1_pp = get_coeff("b_1pp", Z)
    b_2_pp = get_coeff("b_2pp", Z)
    gamma_0 = get_coeff("gamma_0", Z)
    gamma_0_p = get_coeff("gamma_0p", Z)
    gamma_1_p = get_coeff("gamma_1p", Z)
    c_0_p = get_coeff("c_0p", Z)
    c_1_p = get_coeff("c_1p", Z)
    c_2_p = get_coeff("c_2p", Z)
    gamma_0_wedge = get_coeff("gamma_0_wedge", Z)
    gamma_0_pp = get_coeff("gamma_0pp", Z)
    gamma_1_pp = get_coeff("gamma_1pp", Z)
    c_0_pp = get_coeff("c_0pp", Z)
    c_1_pp = get_coeff("c_1pp", Z)
    c_2_pp = get_coeff("c_2pp", Z)

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


if __name__ == "__main__":

    dict = load_EHcoeffs()
    '''
        print dict.keys()
        print '\n\n dict: \n: ', dict
        print '\n\n Z: ', dict['Z']
        
        print '\n\n alpha_0: ', dict['alpha_0']
        print '\n\n alpha_0_prime: ', dict["alpha_0'"]
        Z = 13.0
        chi = 0.1#omega*tau
        t_dict,coeff_dict = coeff_poly_fit(chi,Z)
        print t_dict
        #alpha_0'= coeff'oly_fit(Z)
        #print '\n\nalpha_0' at Z=13= ',alpha_0'
        '''
