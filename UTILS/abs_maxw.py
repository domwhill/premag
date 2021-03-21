#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 17:11:00 2019
Converts laser heating amplitude to a J eV measure
@author: DH814
"""

ne_ref = 5.0e18
Z, Ar = 1.0, 1.0    # Hydrogen
Te_ref = 20.0    # eV
Bz_ref = 1.0    # Tesla
tau_ei = 8.66511e-13
Io_W_cm2 = 1.5527e15
d_a = 1e-3
#heat_amp_xyt = ((Io_W_cm2*1d4)/(d_a))*(2.0d0/(3.0d0*ne(iy,ix)))
heat_amp_xyt = ((Io_W_cm2 * 1e4) / (d_a)) * (2.0 / (3.0 * ne_ref * 1e6))    # J/s

heat_amp_si = heat_amp_xyt * tau_ei / (2.0 * Te_ref)

print('heat_amp_xyt = ', heat_amp_xyt, ' J/s heat_amp_xyt norm si = ', heat_amp_si * 1e12, 'tau/ps')
print('heat amp =', heat_amp_xyt * 1e12 / (q_e), ' eV/ps ')
