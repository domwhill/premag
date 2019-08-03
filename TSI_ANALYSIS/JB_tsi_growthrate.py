#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 16:28:44 2019

@author: dominichill
"""
import numpy as np, getpass, site
userid = getpass.getuser()
site.addsitedir('/Users/'+ userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
import impact_norms_py3 as ip3
import EH_poly_coeff_module_py3 as ep3
import matplotlib.pyplot as plt


c = 3e8
q_e = 1.602e-19
k_B = 1.38e-23
m_e = 9.11e-31

def get_wt_grad(wt,Z,var_list=['alpha_perp','alpha_wedge','beta_perp','beta_wedge','kappa_perp','kappa_wedge']):
    dwt = 0.01
    if wt>dwt:
        wtmin = wt-dwt
        wtmax = wt+dwt
    else:
        wtmin = wt
        wtmax = wt + 2.0*dwt
    t_dict_p1,c_dict = ep3.coeff_poly_fit(wtmax,Z)
    t_dict_m1,c_dict = ep3.coeff_poly_fit(wtmin,Z)
    dict_out = {}
    for var in var_list:
        dict_out[var] = (t_dict_p1[var] - t_dict_m1[var])/(2.0*dwt)
    return dict_out
#---------------------------------------------------------------------------
def get_dimensionful_coeff(val_dimensionless,var,ref_vals,dimensionless=True):
    n_e_m3 = ref_vals['n_e']*1e6
    T_eJ = q_e*ref_vals['T_e']
    tau_B = ref_vals['tau_B']
    if dimensionless == True:
         return val_dimensionless       
    if var == 'kappa':
        val_out = val_dimensionless*(n_e_m3*T_eJ*tau_B)/m_e
    elif var == 'alpha':
        val_out = val_dimensionless*(m_e*n_e_m3/tau_B)
    elif var == 'beta':
        val_out =val_dimensionless
    return val_out
#------------------------------------------------------------------------------
def get_dimensionful_transport_c(wte,Z,ref_vals):
    '''
       dict,grad_dict= get_dimensionful_tranpsort_c(wte,Z)
    '''
     
    transport_dict,c_dict =  ep3.coeff_poly_fit(wte,Z)
    t_grad_dict = get_wt_grad(wte,Z)   
    var_list = ['alpha_para','alpha_perp','alpha_wedge','kappa_para',
                'kappa_perp','kappa_wedge','beta_perp','beta_wedge']
    dict = {}
    grad_dict = {}
    for var in var_list:
        dict[var] = get_dimensionful_coeff(transport_dict[var],var.split('_')[0],ref_vals)
    for var in t_grad_dict.keys():
        grad_dict[var]= t_grad_dict[var]
    #grad_dict['kappa_wedge']  = get_dimensionful_coeff(t_grad_dict['kappa_wedge'],'kappa',ref_vals)
    #grad_dict['beta_wedge'] = t_grad_dict['beta_wedge']     
    dict['psi_wedge'] = transport_dict['beta_wedge'] # ettinghausen term

    return dict,grad_dict

#------------------------------------------------------------------------------
def get_JB_params(Z,ne_cm3,Te_eV,Bz):
    '''
        
    '''
    
    log_lambda = 6.9 - np.log(Z/10.0) + 1.5*np.log(Te_eV/1000.0) - 0.5*np.log(ne_cm3/1e21)
    # ps
    tau_T_ps = (1.0/6.0)*((log_lambda/5.0)**-1)*((Z/10.0)**-1)*((ne_cm3/1e21)**-1)*((Te_eV/1000.0)**1.5)
    # um
    lambda_T_um = 3.0*(5.0/log_lambda)*(10.0/Z)*(1e21/ne_cm3)*((Te_eV/1000.0)**2)
    
    v_th = lambda_T_um*(1e-6)/(tau_T_ps*1e-12)
    omega_g = q_e*Bz/m_e
    rL_um = (v_th/omega_g)*1e6
    chi = cB*(lambda_T_um)*(rL_um**-1)
    
    chi2 = 0.25*(tau_T_ps)*Bz
    Lambda = ((5.48e-2)**-1)*((ne_cm3/1e21)**-0.5)*(10.0/Z)*(5.0/log_lambda)*((Te_eV/1000.0)**2)
    dict = {}
    dict['wte'] = chi*1.0
    dict['wte2'] = chi2*1.0
    
    dict['lambda_T'] = lambda_T_um*1e-6 
    dict['tau_T'] = tau_T_ps*1e-12
    dict['log_lambda'] = log_lambda
    dict['Lambda'] = Lambda
    
    #
    return dict
#------------------------------------
def compute_growth(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB):
    '''
        w_p_2,w_m_2 = compute_growth(k,ne_ref,Te_ref,Z,Bz,6.5)
    '''
    dict_n = ip3.impact_inputs(ne_ref,Te_ref,Z,Bz,Ar) # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei =  dict_n['wpe_over_nu_ei'] 
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']#*cB # eB/m_e *tau_ei
    wte_in = wte*1.0
    delta_c_norm = c_over_vte/wpe_over_nu_ei
    delta_c_norm_2 = 1e3*Z*((ne_ref/1e21)**0.5)*log_lambda/(Te_ref**2)
    delta= delta_c_norm*lambda_T # m
    dict_n['delta'] = delta
    #delta = c/wpe
    #Lambda = lambda_T_in/delta
    Lambda_1 = delta_c_norm**-1
    Lambda_2 = delta_c_norm_2**-1
    #k = (2.0*n[p.pi)/(5.0*lambda_mfp_ref)#7.0e-6#
    
    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T*cB
    
    dict,grad_dict= get_dimensionful_transport_c(wte,Z,ref_vals)
    
    alpha_para = dict['alpha_para']
    alpha_perp = dict['alpha_perp']
    alpha_wedge = dict['alpha_wedge']
    
    
    kappa_para = dict['kappa_para']
    dkappawedge_dchi =grad_dict['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge'] 
    dbetaperp_dchi = grad_dict['beta_perp'] 
    dalphawedge_dchi = grad_dict['alpha_wedge'] 
    
    
    
    beta_wedge = dict['beta_wedge'] # /e?
    psi_wedge = dict['beta_wedge'] # ettinghausen term
    kappa_perp = dict['kappa_perp']
    
    sin_theta = 1.0
    sP = 2.0*beta_wedge*dkappawedge_dchi*sin_theta*(cB**2)*((lambda_T**4)/(3.0*LT*tau_T**2))
    sB = ((cB*wte*lambda_T**2)/(3.0*LB*tau_T**2))*dkappawedge_dchi*sin_theta
    sE = ((4.0*(lambda_T**2)*delta**2)/(3.0*tau_T**2))*beta_wedge*psi_wedge
    dR = alpha_perp*(delta**2)/(cB*tau_T)
    dT = (cB*kappa_perp*lambda_T**2)/ (3.0*tau_T)
    pm = 1.0
    
    a1 = sB*k - (dR + dT)*1j*(k**2)
    b1 = (sB**2)*k**2 
    b2 = (sP + 2.0*sB*(dR - dT)*1j*k**3)
    b3 = - ((dR-dT)**2 + sE)*(k**4)
    coeff_sqrt =  (b1+ b2 + b3)**0.5
    pm = 1.0
    w_p_1 = 0.5*(a1) + pm*0.5*(coeff_sqrt)
    pm = -1.0
    w_m_1 = 0.5*(a1) + pm*0.5*(coeff_sqrt)
    
    #----- compute Bissell coefficients
    #--- check
    aN = (cB/(2.0*wte))*(lambda_T**2/(tau_T))*beta_wedge
    aE = ((2.0*wte)/(3.0*cB))*(delta**2/(tau_T))*beta_wedge
    
    ckappa = (cB/2.0)*wte*dkappawedge_dchi*(lambda_T**2/tau_T)
    cbeta = (cB/2.0)*wte*dkappawedge_dchi*(lambda_T**2/tau_T)
    ckappa = (cB/3.0)*wte*dkappawedge_dchi*(lambda_T**2/tau_T)
    cbeta = (cB/2.0)*dbetaperp_dchi*(lambda_T**2/tau_T)

    aC = (wte/cB)*((delta**2)/(tau_T))*(1.0 + wte*dbetaperp_dchi)
    
    cH_hat = (wte/(cB))*(delta**2/tau_T)*(alpha_wedge/(wte) - dalphawedge_dchi)
    
    sigmaB = 1
    k2 = k**2
    eT = sigmaB*(ckappa/LB + aC/LB)*sin_theta
    eB = sigmaB*(ckappa/LT + aC/LT)*sin_theta
    fT = sigmaB*(cbeta/LB + (3.0*cH_hat/(2.0*LB)))*sin_theta
    fB = sigmaB*(cbeta/LT + 3.0*cH_hat/(2.0*LT))*sin_theta
    vB = (sigmaB/LB)*(ckappa + aC)*sin_theta
    s = 4.0*(aN*1j*k2 + fT*k)*(aE*1j*k2 + eB*k)
    a = (eT+fB)*k- (dT+ dR)*1j*k2
    b = (((dT- dR)*1j*k2-(eT-fB)*k)**2 +s)**0.5
    w_p_2 = 0.5*(a + b)
    w_m_2 = 0.5*(a - b)
    return w_p_1,w_m_1,w_p_2,w_m_2, dict_n
#------------------------------------------------------------
#------------------------------------
def compute_growth_7_32(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB):
    '''
        OMEGA_p,OMEGA_m, dict_n= compute_growth_7_32(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
    '''
    dict_n = ip3.impact_inputs(ne_ref,Te_ref,Z,Bz,Ar) # Bz must be in tesla because this is the units of gorgon data
    #--------
    wpe_over_nu_ei =  dict_n['wpe_over_nu_ei'] 
    c_over_vte = dict_n['c_over_vte']
    log_lambda = dict_n['log_lambda']
    lambda_T = dict_n['lambda_mfp']
    v_T = dict_n['vte']
    tau_T = dict_n['tau_ei']
    wte = dict_n['Bz_norm']*cB # eB/m_e *tau_ei
    delta_c_norm = c_over_vte/wpe_over_nu_ei
    delta_c_norm_2 = 1e3*Z*((ne_ref/1e21)**0.5)*log_lambda/(Te_ref**2)
    delta = delta_c_norm*lambda_T # m
    dict_n['delta'] = delta
    #delta = c/wpe
    # Lambda = lambda_T/delta_c
    Lambda = delta_c_norm**-1
    
    ref_vals = {}
    ref_vals['n_e'] = ne_ref
    ref_vals['T_e'] = Te_ref
    ref_vals['tau_B'] = tau_T*cB
    
    dict,grad_dict= get_dimensionful_transport_c(wte,Z,ref_vals)
    
    alpha_para = dict['alpha_para']
    alpha_perp = dict['alpha_perp']
    alpha_wedge = dict['alpha_wedge']
    
    
    kappa_para = dict['kappa_para']
    kappa_perp = dict['kappa_perp']
    kappa_wedge = dict['kappa_wedge']
    
    beta_perp = dict['beta_perp']
    beta_wedge = dict['beta_wedge'] # /e?
    psi_wedge = dict['beta_wedge'] # ettinghausen term
    psi_perp = dict['beta_perp'] -1.0
    
    dkappawedge_dchi =grad_dict['kappa_wedge']
    dbetawedge_dchi = grad_dict['beta_wedge'] 
    dpsiperp_dchi = grad_dict['beta_perp'] 

    dbetaperp_dchi = grad_dict['beta_perp'] 
    dalphawedge_dchi = grad_dict['alpha_wedge'] 
    
    
    sigmaB = 1
    sin_theta = 1.0

    #--- compute source terms ----> copied form previous
    
    sP = 2.0*beta_wedge*dkappawedge_dchi*sin_theta*(cB**2)*((lambda_T**4)/(3.0*LT*tau_T**2))
    sB = ((cB*wte*lambda_T**2)/(3.0*LB*tau_T**2))*dkappawedge_dchi*sin_theta
    sE = ((4.0*(lambda_T**2)*delta**2)/(3.0*tau_T**2))*beta_wedge*psi_wedge
    #---->

    
    #----- compute Bissell coefficients
    #--- check
    AN2 = (cB/(2.0*wte))*beta_wedge
    AE2 = ((2.0*wte)/(3.0*cB))*(beta_wedge/(Lambda**2))
    
    ckappa = (cB/3.0)*wte*dkappawedge_dchi*(lambda_T**2/tau_T)
    cbeta = (cB/2.0)*dbetaperp_dchi*(lambda_T**2/tau_T)
    aC = (wte/cB)*((delta**2)/(tau_T))*(1.0 + wte*dbetaperp_dchi)
    vB = (sigmaB/LB)*(ckappa + aC)*sin_theta
    
    
    cH_hat = (wte/(cB))*(delta**2/tau_T)*(alpha_wedge/(wte) - dalphawedge_dchi)
    
    # Eq. 7.36 p. 155
    omega_g = q_e*Bz/m_e
    Omega = omega_g*tau_T
    K = k*lambda_T
    K2 = K**2
    LT /= lambda_T
    LB /= lambda_T
    #----->
    DT = (cB*kappa_perp)/ (3.0)
    DR = alpha_perp/(cB*(Lambda**2))
    AN = cB*beta_wedge/(2.0*wte)
    AE = 2.0*wte*beta_wedge/(3.0*cB*Lambda**2)
    Ckappa = (cB*wte/3.0)*dkappawedge_dchi
    #------>
    AC = (2.0*wte/(3.0*cB*Lambda**2))*(psi_perp + 2.5 + 1.5*wte*dpsiperp_dchi - beta_perp)
    AC2 = (wte/(cB*Lambda**2))*(1.0 + wte*dpsiperp_dchi)#(2.0*wte/(3.0*cB*Lambda**2))*(psi_perp + 2.5 + 1.5*wte*dpsiperp_dchi - beta_perp)
    
    VB2 = vB/v_T
    SP2 = sP*tau_T**2/(lambda_T**3)
    SP = (4.0*AN*sigmaB/LT)*(Ckappa+AC)*sin_theta
    VB = SP*LT/(4.0*AN*LB)
    SE = sE*tau_T**2/(lambda_T**4)
    SE2 = 4.0*AN*AE
    # --- I got here...
    A1 =  VB*K-(DT+DR)*1j*K2
    B1 = ((DT-DR)*1j*K2 - VB*K)**2 
    B2 = SP*1j*(K**3) - SE*(K**4)
    B = np.sqrt(B1 + B2)
    
    OMEGA_p = 0.5*(A1 + B)
    OMEGA_m = 0.5*(A1 - B)
    #-- here
    print('AN = ', AN, AN2,'AE = ',AE,AE2)
    print('AC = ', AC, AC2)
    print('VB = ', VB, VB2)
    print('SP = ', SP, SP2)
    print('SE = ', SE, SE2)
    
    return OMEGA_p/tau_T,OMEGA_m/tau_T, dict_n
#----------------------------------------------------------------------------
#dict_all = cf.load_data_all(path,fprefix,time)
c = 3e8
cB = 3.0*np.sqrt(np.pi)/4.0
# User inputs - JB params PRL 2010
##Te_ref,ne_ref,Z,Bz_ref = norms
##k = 0.0-1.0 um
Te_ref = 20.0 # ev
ne_ref = 1.6e19 # cm^{-3}
LN = 100.0e-6 #n_e_in/dxn_in
LT = 100.0e-6#T_e_in/dxT_in
LB = -1.0*LT#10.0e-6 #n_e_in/dxn_in

# Table 7.2 - p. 159 J.B. thesis
if True:
    Bz = 6.0
    Z = 7.0
    Ar = 14.0
    Te_ref = 0.398*1000.0 # eV
    ne_ref = 0.0149e21 # cm^{-3}
    LT = -160.0e-6
    LB = 146.0e-6
else:
    Bz = 12.0
    Z = 7.0
    Ar = 14.0
    Te_ref = 0.436*1000.0 # eV
    ne_ref = 0.0142e21 # cm^{-3}
    LT = -152.0e-6
    LB = 152.0e-6
    LN = 235.0e-6    
dict_JB = get_JB_params(Z,ne_ref,Te_ref,Bz)

# 1/LT,B << |k| << 1/\lambda_T

k = np.linspace(0.0,0.9,1000)*(1e6) # 0-1 um, 100
omega_p,omega_m,w_p_2,w_m_2,dict_n = compute_growth(k,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()

s_factor = LT*LN
ax.plot(k*1e-6,(omega_p.imag)*1e-9,c='b')
#ax.plot(k*1e-6,(omega_m.imag),c='r')
ax2.plot(k*1e-6,w_p_2.imag*1e-9,c='b',linestyle='--')

ax.set_ylabel('$\omega$ [$ns^{-1}$]' )
ax2.set_ylabel('$\omega$ [$ns^{-1}$]' )
ax.set_title('Magnetothermal Growth rate - recreation of Fig. 7.6, p.149, JB Thesis')
#ax2.plot(k*1e-6,w_m_2.imag,c='r',linestyle='--')
ymin,ymax = ax2.get_ylim()
ax2.set_ylim(0.0,ymax)
ax.set_xlabel(r'k [$\mu m^{-1}$]')
#plt.show()

#----------------------------------------------
#----
##Froula Te 0.4,4T, 0.015e21, 7.0 LT = 150.0e-6, wte = 2.5, Lambda = 20.0
#Froula
Te_ref = 0.4*1000.0
Bz = 4.0
ne_ref = 0.015e21
Z = 7.0
LT = -150.0e-6
LB = -1.0*LT
wte = 2.5
k_f = 2.0*np.pi/(40.0e-6)
omega_pf,omega_mf,w_p_2_f,w_m_2_f,dict_nf = compute_growth(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
OMEGA_P,OMEGA_M,dict_n = compute_growth_7_32(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
omega_exp = 20.0
print('1) FROULA GAMMA_P[9] = %3.4f ns^{-1}  gamma_m  = %3.4f ns^{-1}--> we expect %3.4f ns^{-1}' % (w_p_2_f.imag*1e-9,w_m_2_f.imag*1e-9,omega_exp))
print('2) FROULA GAMMA_P [9] = %3.4f ns^{-1} = %3.4f ns^{-1} --> we expect %3.4f ns^{-1}' % (OMEGA_P.imag*1e-9,OMEGA_M.imag*1e-9,omega_exp))

#Froula [10]
Te_ref = 0.2*1000.0
Bz = 3.0
ne_ref = 0.00075e21
Z = 2.0
LT = -200.0e-6
LB = -1.0*LT
wte = 40.0
k_f = 2.0*np.pi/(25.0e-6)
oemga_exp = 10.0 # ns^{-1}

omega_pf,omega_mf,w_p_2_f,w_m_2_f,dict_nf = compute_growth(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
OMEGA_P,OMEGA_M,dict_n = compute_growth_7_32(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
print('1) FROULA GAMMA_P[10] = %3.4f ns^{-1}  gamma_m  = %3.4f ns^{-1}--> we expect %3.4f ns^{-1}' % (w_p_2_f.imag*1e-9,w_m_2_f.imag*1e-9,omega_exp))
print('2) FROULA GAMMA_P [10] = %3.4f ns^{-1} = %3.4f ns^{-1} --> we expect %3.4f ns^{-1}' % (OMEGA_P.imag*1e-9,OMEGA_M.imag*1e-9,omega_exp))

Te_ref = 0.4*1000.0
Bz = 10.0
ne_ref = 0.0035e21
Z = 3.5
LT = -100.0e-6
LB = -1.0*LT
wte = 40.0
k_f = 2.0*np.pi/(10.0e-6)
omega_exp = 25.0 # ns^{-1}
exp_string = 'Li et al. [106]'
omega_pf,omega_mf,w_p_2_f,w_m_2_f,dict_nf = compute_growth(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
OMEGA_P,OMEGA_M,dict_n = compute_growth_7_32(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)

print('Li et al [106] = ', w_p_2_f.imag*1e-9,'ns^{-1}',omega_pf.imag*1e-9, '---> we wexpect --> 25 ns^{-1}')
print('1) %s GAMMA_P= %3.4f ns^{-1}  gamma_m  = %3.4f ns^{-1}--> we expect %3.4f ns^{-1}' % (exp_string,w_p_2_f.imag*1e-9,w_m_2_f.imag*1e-9,omega_exp))
print('2) %s GAMMA_P= %3.4f ns^{-1} = %3.4f ns^{-1} --> we expect %3.4f ns^{-1}' % (exp_string,OMEGA_P.imag*1e-9,OMEGA_M.imag*1e-9,omega_exp))
#---->
Te_ref = 2.5*1000.0
Bz = 40.0
ne_ref = 0.25e21
Z = 2.0
LT = -300.0e-6
LB = -1.0*LT
wte = 60.0
k_f = 2.0*np.pi/(5.0e-6)
omega_exp = 50.0 # ns^{-1}
exp_string = 'ICF. [31]'
omega_pf,omega_mf,w_p_2_f,w_m_2_f,dict_nf = compute_growth(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
OMEGA_P,OMEGA_M,dict_n = compute_growth_7_32(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)

print('1) %s GAMMA_P= %3.4f ns^{-1}  gamma_m  = %3.4f ns^{-1}--> we expect %3.4f ns^{-1}' % (exp_string,w_p_2_f.imag*1e-9,w_m_2_f.imag*1e-9,omega_exp))
print('2) %s GAMMA_P= %3.4f ns^{-1} = %3.4f ns^{-1} --> we expect %3.4f ns^{-1}' % (exp_string,OMEGA_P.imag*1e-9,OMEGA_M.imag*1e-9,omega_exp))
#---->---->
#----> MIF
Te_ref = 1.0*1000.0
Bz = 1000.0
ne_ref = 200.0e21
Z = 1.0
LT = -5.0e-6
LB = 2.0*LT
wte = 2.0
k_f = 2.0*np.pi/(5.0e-6)
omega_exp = 15.0 # ns^{-1}
exp_string = 'MIF. [8,7]'
omega_pf,omega_mf,w_p_2_f,w_m_2_f,dict_nf = compute_growth(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)
OMEGA_P,OMEGA_M,dict_n = compute_growth_7_32(k_f,ne_ref,Te_ref,Z,Bz,Ar,LT,LB)

print('1) %s GAMMA_P= %3.4f ns^{-1}  gamma_m  = %3.4f ns^{-1}--> we expect %3.4f ns^{-1}' % (exp_string,w_p_2_f.imag*1e-9,w_m_2_f.imag*1e-9,omega_exp))
print('2) %s GAMMA_P= %3.4f ns^{-1} = %3.4f ns^{-1} --> we expect %3.4f ns^{-1}' % (exp_string,OMEGA_P.imag*1e-9,OMEGA_M.imag*1e-9,omega_exp))
print('#---->---->\n\n')


wpe_over_nu_ei =  dict_n['wpe_over_nu_ei'] 
c_over_vte = dict_n['c_over_vte']
log_lambda = dict_n['log_lambda']
lambda_T = dict_n['lambda_mfp']
v_T = dict_n['vte']
tau_T = dict_n['tau_ei']
delta = dict_n['delta']
wte = dict_n['Bz_norm']
Lambda = lambda_T/delta

print(' lambda_T = ', lambda_T, '<-- check --> ', dict_JB['lambda_T'])
print(' tau_T = ', tau_T, '<-- check --> ', dict_JB['tau_T'])
print(' logLambda = ', log_lambda, '<--- check ---> ', dict_JB['log_lambda'])
print(' delta = ', delta, '<--- check ---> ', dict_JB['Lambda']*dict_JB['lambda_T'])

print('JB table 8')
Te_JB = 0.392
ne_JB = 0.0149
wte_JB = 2.4
Lambda_JB = 20.5
lambda_T_um_JB = 28.0
tau_T_ps_JB = 2.4
parse_str = (Te_JB,ne_JB,wte_JB, Lambda_JB,lambda_T_um_JB,tau_T_ps_JB)

parse_str2 = (Te_ref*1e-3,ne_ref/1e21,dict_JB['wte'], dict_JB['Lambda'],dict_JB['lambda_T']*1e6,dict_JB['tau_T']*1e12)
parse_str3 = (Te_ref*1e-3,ne_ref/1e21,wte,Lambda,lambda_T*1e6,tau_T*1e12)

print(' Te_keV = %3.4f \t ne/10^21 cm^{-3} = %3.4f \t  wte = %3.4f \t  Lambda = %3.4f  \t lambda_T = %3.4f um  \t tau_T = %4.3f ps' % parse_str)
print('--- calc from JB appendix ====')
print(' Te_keV = %3.4f \t ne/10^21 cm^{-3} = %3.4f \t  wte = %3.4f \t  Lambda = %3.4f  \t lambda_T = %3.4f um  \t tau_T = %4.3f ps' % parse_str2)
print('--- var in===')
print(' Te_keV = %3.4f \t ne/10^21 cm^{-3} = %3.4f \t  wte = %3.4f \t  Lambda = %3.4f  \t lambda_T = %3.4f um  \t tau_T = %4.3f ps' % parse_str3)
#print('Lambda_1 ={} <-should be same as--> {} '.format(Lambda_1,Lambda_2))
plt.show()
'''
STATUS : - lambda_T etc. don't exactly match JB values, however we get similar curves...
- check Z is correct
- check lambda_T etc. against JB appendix 2.
'''
