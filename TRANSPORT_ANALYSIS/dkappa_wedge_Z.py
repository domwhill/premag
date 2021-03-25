#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
var dkappa_wedge/dchi
@author: DH814
"""
#---------> MTI INSTABILITY GROWTH RATE CALCULATION
import numpy as np, getpass, site

userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')
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
time = '06'
norm_name = 'r_5v37'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + '/norm.log'
[Te_ref, ne_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
Te_ref *= q_e    # convert to Joules
cd5 = cf.ConversionFactors(norm_path, Z_ref, Ar=6.51)
cB = 3.0 * np.sqrt(np.pi) / 4.0
lambda_mfp_ref = cd5.lambda_mfp
v_th_ref = cd5.v_te
tau_ref = cd5.tau_ei
Bz_ref = cd5.Bz_ref
#----##
# fix Bz
Bz = 0.1

nz = 50
nwt = 20
Z_list = np.linspace(0.1, 5.0, nz)
Bz_list = 10**(np.linspace(-2.0, 2.0, nwt))
kwedge_arr = np.zeros((nwt, nz))

for ibz in range(len(Bz_list)):
    for iz in range(len(Z_list)):
        Z_in = Z_list[iz]
        Bz_in = Bz_list[ibz]
        transport_dict, c_dict = ep3.coeff_poly_fit(Bz_in / Z_in, Z_in)
        kwedge = transport_dict['kappa_wedge']

        kwedge_arr[ibz, iz] = kwedge

fig = plt.figure()
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax.imshow(kwedge_arr, extent=[Z_list[0], Z_list[-1], Bz_list[0], Bz_list[-1]], aspect='auto')
ax.set_xlabel('Z')
ax.set_ylabel(r'$B_0$')

p_list = []
leg_list = []
for iz in [0, 1, 4, 8, 10, 15]:
    p1, = ax2.loglog(Bz_list, kwedge_arr[:, iz])
    p_list.append(p1)
    leg_list.append('%1.1f' % Z_list[iz])
ax2.legend(p_list, leg_list)
plt.show()
