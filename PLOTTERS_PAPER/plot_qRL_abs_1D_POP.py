"""Plot 1D lineouts of transport terms for kinetic and classical predictions.

Usage example:

    python PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py <transport_component> <impact variable to plot> <save_directory>
eg.
 $  python2 PLOTTERS_PAPER/plot_qRL_abs_1D_POP.py 'RL y' 'Bz' 'OUTPUT/'

Transport components that can be plotted:
    'SH x':  diffusive heat flow x component
    'SH y':  diffusive heat flow y component
    'RL x':  Righi-Leduc heat flow x component
    'RL y':  Righi-Leduc heat flow y component
    'E x':   Ettinghausen heat flow x component
    'E y':   Ettinghausen heat flow y component
    'vN x':  Nernst velocity x component
    'vN y':  Nernst velocity y component


Also plots for reference 1D profiles of IMPACT variables
{Te: electron temperature,
ne: electron number density,
Bz: magnetic field (z/out of plane) component
}

"""

import sys, re, os, getpass, site, pdb
userid = getpass.getuser()
from pylab import *
from matplotlib.legend_handler import HandlerBase
import argparse
sys.path.extend(["./"])
import MODULES.chfoil_module as cf
import MODULES.figure_prl_twocol as fprl
import MODULES.house_keeping as hk
import MODULES.kinetic_ohmslaw_module_1D_varZ as q_mod
q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12

#---> file inputs
paths = hk.directory_paths()
src_dir = paths.src_dir
data_dir = paths.data_dir_2D
save_path = paths.save_dir
norm_dir = paths.norm_dir
log_file = norm_dir + 'norm.log'
[Te_ref, n_ref, Z_ref, Bz_ref] = np.loadtxt(log_file)
cd5 = cf.conv_factors_custom(norm_dir, Z_ref, Ar=6.51)
#---> USER inputs
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "transport_component",
    description=
    "Transport term to plot RL x = x component of Righi-Leduc, SH y = y component of Spitzer-HArm/diffusive heat flow,"
    "vN y = y-component of Nernst")
parser.add_argument("thermodynamic_variable",
                    description="string for impact variable: Te, ne, Bz, Cx (ion velocity)")
parser.add_argument("output_directory",
                    default="",
                    description="output directory path (relative to root)")
args = parser.parse_args()
var1 = args.transport_component    # eg. 'RL y' for y component of Righi-Leduc heat flow, '
var2 = args.thermodynamic_variable
save_path = args.output_directory
print('Input argument= ', args)
if save_path[-1] != '/':
    save_path = save_path + '/'
print("var1 = " + var1 + " var2 = " + var2 + " output_directory = " + save_path)

save_name = '%sqlineout_%s_abs_%s_1d.png' % (save_path, re.sub(' ', '', var1), var2)
#----------------

cl_index = int(cd5.cl_index)
c_index = int(cd5.c_index)
cl_index, c_index = 0, -1
SI_on = cd5.SI_on
tau_ei = cd5.tau_ei
nu_ei = cd5.nu_ei
lambda_mfp = cd5.lambda_mfp
v_te = cd5.v_te
xstep_factor = cd5.xstep_factor
tstep_factor = cd5.tstep_factor
xlab = cd5.xlab_rel
ylab = cd5.ylab
Bz_ref = (m_e / (q_e * tau_ei))

#<--------------
# inputs
#------------->
scale_len = 1
dim = 1
lambda_p_mfp = 5    # doesn't matter since looking for 2d runs
path_0T = paths.get_path(scale_len, '0', lambda_p_mfp, dim=dim)
path_50T = paths.get_path(scale_len, '50', lambda_p_mfp, dim=dim)
path_400T = paths.get_path(scale_len, '400', lambda_p_mfp, dim=dim)
print('path 0 = #', path_0T, '\n', path_50T, '\n', path_400T)
Nernst_on = True
time = '15'
iy = 1
xmin, xmax = -5.0, 20.0
thresh_N_ratio = 5e-2

var_list = ['SH x', 'SH y', 'RL x', 'RL y', 'E x', 'E y', 'vN x', 'vN y']
lab_list = {
    'SH x': r'$q_{\perp,x}$',
    'SH y': r'$q_{\perp,y}$',
    'RL x': r'$q_{RL,x}$',
    'RL y': r'$q_{RL,y}$',
    'E x': r'$q_{E,x}$',
    'E y': r'$q_{E,y}$',
    'vN x': r'$v_{N,x}$',
    'vN y': r'$v_{N,y}$'
}

# Boolean to determine whether the 400T kinetic data should set axis ylims
# or the 50T classical
if var1[:2] == 'RL' or var1[0] == 'E':
    var_lim_400 = True
else:
    var_lim_400 = False


#<-------- norm_labels
def onesigfig(val):
    return float('%0.1g' % val)


def unit_item(val, unit):
    red_val = '%1.1e' % onesigfig(val)
    if red_val[0] != '1':
        lab = r'$\,[\SI{%1.1e}{%s}]$' % (val, unit)
    else:
        pow = int(red_val.split('e')[-1])
        if pow != 0:
            lab = r'$\,[10^{%i}\,\si{%s}]$' % (pow, unit)
        else:
            lab = r'$\,[\si{%s}]$' % (unit)
    return lab


#--------------------------------------------->
class norm_labels(object):

    def __init__(self, var_in, var_data=[]):
        if var_in == 'Bz':
            self.lab = r'$B_z$[\si{T}]'
            self.factor = Bz_ref

        elif var_in == 'Te':
            multT = 2.0
            self.lab = r'$T_e$[\SI{%1.1f}{\kilo\electronvolt}]' % (onesigfig(
                (multT**-1) * 2.0 * Te_ref * 1e-3))
            self.factor = 1.0
        elif var_in in var_list:
            if var_in[0] == 'v':
                multv = 1.0
                self.factor = multv * (cd5.v_te * (1e6) * (1e-12))    #16.0#cd5.v_te*(1e6)*(1e-12)
                v_unit = unit_item(1.0 / multv, unit=r'\micro\meter/ps')
                self.lab = lab_list[var1] + v_unit

            else:

                norm_factor = (n_ref * 1e6) * m_e * (v_te**3)
                if len(var_data) != 0:
                    var_lim = self.get_lim(norm_factor * var_data)
                    pow = int(('%1.1e' % (var_lim)).split('e')[-1])

                else:
                    pow = 18
                multQ = 10**-pow
                self.factor = norm_factor * multQ
                # kg s^-3 => kg m^-1 s^-2
                q_unit = r'$\,[\SI{%0.1g}{W/m^2}]$' % (multQ**-1)
                q_unit = unit_item(multQ**-1, 'W/m^2')
                self.lab = lab_list[var1] + q_unit

    def get_lim(self, var_data):
        vmin, vmax = np.min(var_data), np.max(var_data)
        var_lim = (np.max(np.array([np.abs(vmin), np.abs(vmax)])))
        return var_lim


class AnyObjectHandler(HandlerBase):

    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.7 * height, 0.7 * height],
                        linestyle=orig_handle[1],
                        color='k')
        l2 = plt.Line2D([x0, y0 + width], [0.3 * height, 0.3 * height], color=orig_handle[0])
        return [l1, l2]


#--- functions


def repack(path, time):
    '''
        dict_c,dict_k = repack(path,time)
    '''
    dict_qc, dict_qk = q_mod.get_q_abs(path,
                                       time)    #get_q_ratio(path,time)#get_q_individ(path,time)
    v_nx_k, v_ny_k, v_nx_c, v_ny_c = q_mod.get_Nernst_abs(path, time)

    dict_out_c = {}
    dict_out_k = {}
    for var in var_list:
        if var[0] != 'v':
            dict_out_k[var] = dict_qk[var][:, 0]
        elif var == 'vN x':
            dict_out_k[var] = v_nx_k[:, iy]
        elif var == 'vN y':
            dict_out_k[var] = v_ny_k[:, iy]

    for var in var_list:
        if var[0] != 'v':
            dict_out_c[var] = dict_qc[var][:, 0]
        elif var == 'vN x':
            dict_out_c[var] = v_nx_c[:, iy]
        elif var == 'vN y':
            dict_out_c[var] = v_ny_c[:, iy]

    dict_out_k['U'] = dict_qk['U'][:, iy]
    dict_out_k['Te'] = dict_qk['Te'][:, iy]
    dict_out_k['Bz'] = dict_qk['Bz'][:, iy]
    dict_out_k['x_grid'] = dict_qk['x_grid'] * xstep_factor

    return dict_out_c, dict_out_k


def set_ylim_max_2data(ax_in, grid, data1, data2, y_mult=[1.0, 1.0], xlim=[-5.0, 20.0]):
    xmin, xmax = xlim[0], xlim[1]

    ymin1 = np.min(data1[(grid < xmax) * (grid >= xmin)])
    ymax1 = np.max(data1[(grid < xmax) * (grid >= xmin)])

    ymin2 = np.min(data2[(grid < xmax) * (grid >= xmin)])
    ymax2 = np.max(data2[(grid < xmax) * (grid >= xmin)])
    ymin = min(ymin1, ymin2)
    ymax = max(ymax1, ymax2)

    dy = 0.05 * np.abs(ymax - ymin)
    if ymin != ymax:

        ax_in.set_ylim(ymin - dy * y_mult[0], ymax + dy * y_mult[1])
    else:
        print('error ymin = ymax data all the same!')
    print('ymin = %4.4e ymax = %4.4e ' % (ymin, ymax))

    return


if __name__ == "__main__":
    #----> data loading
    #dict_0T_c,dict_0T_k = repack(path_0T,time)
    dict_50T_c, dict_50T_k = repack(path_50T, time)
    dict_400T_c, dict_400T_k = repack(path_400T, time)

    x_grid = dict_50T_k['x_grid']

    #<--------------- 0T loading---
    #----> plotting
    fig = fprl.newfig_generic_twinx(1.0)
    ax = fig.add_subplot(111)
    ax2 = ax.twinx()
    ax.tick_params(which='both', direction='in')
    ax2.tick_params(which='both', direction='in')

    data_k = dict_50T_k[var1]

    if var_lim_400:
        v1_labs = norm_labels(var1, var_data=dict_400T_k[var1])
    else:
        v1_labs = norm_labels(var1, var_data=dict_50T_c[var1])
    v2_labs = norm_labels(var2)

    vmin = np.min(data_k[np.abs(x_grid - xmin) <= xmax]) * 0.8
    vmax = np.max(data_k[np.abs(x_grid - xmin) <= xmax]) * 1.2

    #p0, = ax.plot(x_grid,rat_0T,c='g')

    p50, = ax.plot(x_grid, dict_50T_k[var1] * v1_labs.factor, c='r')
    p50, = ax.plot(x_grid, dict_50T_c[var1] * v1_labs.factor, c='r', linestyle='--')

    p400, = ax.plot(x_grid, dict_400T_k[var1] * v1_labs.factor, c='b')
    p400, = ax.plot(x_grid, dict_400T_c[var1] * v1_labs.factor, c='b', linestyle='--')

    pT50, = ax2.plot(x_grid, dict_50T_k[var2] * v2_labs.factor, c='k', linestyle='--')
    pT400, = ax2.plot(x_grid, dict_400T_k[var2] * v2_labs.factor, c='k', linestyle=':')

    leg_list = ['$\SI{50}{T}$', '$\SI{400}{T}$']
    plt.legend([("r", "--"), ("b", ":")],
               leg_list,
               handler_map={tuple: AnyObjectHandler()},
               loc='lower right',
               frameon=False)

    ax.set_xlabel(xlab)
    ax.set_ylabel(v1_labs.lab)
    ax2.set_ylabel(v2_labs.lab)
    #ax.set_ylim(-0.15*v1_labs.factor,0.01*v1_labs.factor)#(-0.1,0.01) - Nernst
    if var_lim_400:
        set_ylim_max_2data(ax,
                           x_grid,
                           dict_400T_c[var1] * v1_labs.factor,
                           dict_400T_k[var1] * v1_labs.factor,
                           y_mult=[1.0, 1.0],
                           xlim=[xmin, xmax])

    else:
        set_ylim_max_2data(ax,
                           x_grid,
                           dict_50T_c[var1] * v1_labs.factor,
                           dict_50T_k[var1] * v1_labs.factor,
                           y_mult=[1.0, 1.0],
                           xlim=[xmin, xmax])
    #ax2.set_ylim(0.0,1000.0)

    ax.set_xlim(xmin, xmax)
    plt.savefig(save_name, dpi=600)
    print('saving as: ', save_name)
    print(' copy and paste: open -a preview ' + save_name)
    os.system('open -a preview ' + save_name)
