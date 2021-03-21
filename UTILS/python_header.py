import re, os, sys, site, getpass, numpy as np
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/SIM_DATA/ANALYSIS')
site.addsitedir(path_pre)
import MODULES.chfoil_module as cf

SI_on = cd5_switches.SI_on
save_on = cd5_switches.save_on
hlines_on = cd5_switches.hlines_on
grid_on = cd5_switches.grid_on
horizontal_on = cd5_switches.horizontal_on
separate_plots_on = cd5_switches.separate_plots_on
'''
cwdpath = os.getcwd() +'/'+ path_pre
print cwdpath
#cd5 = cf.conv_factors_custom(cwdpath)
cd5 = cf.conv_factors_premag()
Z0 = cd5.Z
T0 = cd5.T0
n0 = cd5.n0
'''
norm_name = 'p400nFL_5v37/'
norm_path = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/p400nFL_5v37'
log_file = '/Users/' + userid + '/Dropbox/York/Pre-magnetised/gorgon_import-11.10.17/' + norm_name + 'norm.log'
[T0, n0, Z0, Bz0] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_path, Z0, Ar=6.51)

print 'Z0 = ', Z0
print 'T0 = ', T0
print 'n0 = ', n0
cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index, c_index = 0, 500
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab
xlab_rel = cd5.xlab_rel
ylab = cd5.ylab
t_unit = cd5.t_unit
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
lh_style = cd5.lh_style
dashes = cd5.dashes

ylab = cd5.ylab
leg_title = cd5.leg_title
color_lineout = cd5.color_lineout
lineout_list = cd5.lineout_list
yi = cd5.yi

divq_factor = cd5.divq_factor
divq_unit = cd5.divq_unit
#ptag = cf.retrieve_path_tag(path)
ptag = path_tag
fcmap = fprl.plotting_params.lineouts_cmap
#--- saving ----
ratios_savename = 'qSHqRLvN_ratios_' + ptag + 'DDFIW'
