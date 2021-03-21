'''
    Module file ....
    
    
    3/8/2019 - ported from 
'''
import sys, os, getpass, re
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import matplotlib.ticker as ticker
import matplotlib as mpl
from scipy import interpolate as SI
import pdb
import house_keeping as hk
import impact_norms as INN
#-----------------------------------------------------------------------
userid = getpass.getuser()
paths = hk.directory_paths()


class MidPointNorm(Normalize):

    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self, vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0)    # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax), mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat > 0] /= abs(vmax - midpoint)
            resdat[resdat < 0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

#-----------------------------------------------------------------------

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if mpl.cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val - 0.5)
            val[val > 0] *= abs(vmax - midpoint)
            val[val < 0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0:
                return val * abs(vmin - midpoint) + midpoint
            else:
                return val * abs(vmax - midpoint) + midpoint


#-----------------------------------------------------------------------
class MyDict(dict):
    pass


#-----------------------------------------------------------------------


class DraggableLegend:

    def __init__(self, legend):
        self.legend = legend
        self.gotLegend = False
        legend.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
        legend.figure.canvas.mpl_connect('pick_event', self.on_pick)
        legend.figure.canvas.mpl_connect('button_release_event', self.on_release)
        legend.set_picker(self.my_legend_picker)

    def on_motion(self, evt):
        if self.gotLegend:
            dx = evt.x - self.mouse_x
            dy = evt.y - self.mouse_y
            loc_in_canvas = self.legend_x + dx, self.legend_y + dy
            loc_in_norm_axes = self.legend.parent.transAxes.inverted().transform_point(
                loc_in_canvas)
            self.legend._loc = tuple(loc_in_norm_axes)
            self.legend.figure.canvas.draw()

    def my_legend_picker(self, legend, evt):
        return self.legend.legendPatch.contains(evt)

    def on_pick(self, evt):
        if evt.artist == self.legend:
            bbox = self.legend.get_window_extent()
            self.mouse_x = evt.mouseevent.x
            self.mouse_y = evt.mouseevent.y
            self.legend_x = bbox.xmin
            self.legend_y = bbox.ymin
            self.gotLegend = 1

    def on_release(self, event):
        if self.gotLegend:
            self.gotLegend = False


class constants(object):
    q_e = 1.602e-19
    m_e = 9.11e-31
    m_p = 1.67e-27
    k_b = 1.38e-23
    epsilon0 = 8.854e-12


class cd5_switches(object):
    SI_on = True
    save_on = True
    hlines_on = True
    grid_on = True
    horizontal_on = True
    separate_plots_on = True
    marker_on = False


#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
class conv_factors_custom(object):
    # --- arguments given to the class are passed straight to init...
    #
    def __init__(self, norm_path=paths.norm_dir, Z=5.954127788543701172, Ar=6.51):
        self.norm_name = self.search_path_class('norm.txt', norm_path)
        print 'self.norm_name = ', self.norm_name, norm_path
        if norm_path[-1] != '/':
            norm_path = norm_path + '/'
        if os.path.exists(norm_path + 'indices.txt'):
            icl, icu = np.loadtxt(norm_path + 'indices.txt')
            self.cl_index = int(icl)    #80#
            self.c_index = int(icu)    #207#
        else:
            self.cl_index = 0
            self.c_index = 500
        print ' norm_name =', self.norm_name, norm_path
        print 'indices = ', self.cl_index, self.c_index
        #Z,Ar = 3.5,6.51
        # self.var are instance variables specific to this instance of the class

        self.SI_on = cd5_switches.SI_on
        self.Z = Z
        self.Ar = Ar
        [self.T0, self.n0, self.logLambda] = np.loadtxt(norm_path + self.norm_name)
        dict_norm = INN.impact_inputs(self.n0, self.T0, Z, 1.0, Ar)
        self.v_te = dict_norm['vte']    #     m/s  #1.93e7 # ms^-1
        self.tau_ei = dict_norm['tau_ei']    # s   7.44e-14 # s
        self.nu_ei = 1.0 / self.tau_ei    # s^-1
        self.lambda_mfp = dict_norm['lambda_mfp']    # m
        self.Bz_ref = dict_norm['Bz_norm']
        self.Bz0 = dict_norm['Bz_norm']

        self.yi = 10
        self.lineout_list = [20, 40, 50, 73]
        #self.lineout_list = [20,47,73]
        #self.lineout_list = [110,144,190]
        self.lineout_list = map(int, np.linspace(self.cl_index, self.c_index, 3))
        self.lineout_list = [90, 120, 140]

        interval = self.c_index - self.cl_index
        #self.lineout_list = map(int,[self.cl_index+0.1*interval,self.cl_index+0.3*interval,self.c_index*0.7])
        print 'lineout list = ', self.lineout_list

        self.color_lineout = cm.bone(np.linspace(0.0, 1, len(self.lineout_list) + 1))
        self.color_lineout = self.color_lineout[:-1]
        self.color_lineout = ['b', 'g', 'r']

        self.lh_style = 'None'
        self.dashes = (4, 2)
        self.double_plot_scale_ratio = 0.8

        if self.SI_on:
            self.xstep_factor = self.lambda_mfp * 1e6
            self.tstep_factor = self.tau_ei * 1e12
            self.xlab = r'$x$ [ $\si{\micro\meter}$ ]'
            self.xlab_rel = r'$x-x_{abl}$ [ $\si{\micro\meter}$ ]'

            self.ylab = r'$y$ [ $\si{\micro\meter}$ ]'
            self.leg_title = r'time [ $\si{\pico\second}$ ]'
            self.divq_unit = r'[$\si{\eV}/\si{\pico\second}$]'

            self.divq_factor = (self.T0 / (self.tau_ei * 1e12))
            self.t_unit = r'[ $\si{\pico\second}$ ]'
            self.tlab = r't ' + self.t_unit

        else:
            self.xstep_factor = 1.0
            self.tstep_factor = 1.0
            self.xlab = r'\textit{x-x_{abl} [ $\lambda_{mfp}$ ]}'
            self.ylab = r'\textit{x [ $\lambda_{mfp}$ ]}'
            self.leg_title = r'\text{$time$ [t$_{col}$ ]}'
            self.divq_factor = 1.0
            self.divq_unit = r' [$T_n/ t_{col}$ ]'
            self.t_unit = r'[ $t_{col}$ ]'
            self.tlab = r't ' + self.t_unit
        #----> switches --->
        self.SI_on = True
        self.save_on = True
        self.hlines_on = True
        self.grid_on = True
        self.horizontal_on = True
        self.separate_plots_on = True
        self.marker_on = False
        self.MSIZE = 4
        #<----

        return

    #-------------------------------------------------------------------
    def search_path_class(self, regexp, path):
        '''
            out_name = search_path(regexp,path)
        '''
        out_name = None
        for pp in os.listdir(path):
            rr = re.search(regexp, pp)
            if rr:
                print 'path = ', pp
                out_name = pp
                break
        if out_name is None:
            raise ValueError("File not found for regex in directory")

        return out_name


#-------------------------------------------------------------------====
def fpre(path):
    return path.split('/')[-1]


#-----------------------------------------------------------------------


def list_to_float(input):
    list = input.split()
    arr = np.zeros((len(list)), dtype=float)
    arr[:] = list[:]
    return arr


#-----------------------------------------------------------------------


def construct_fname(path, fprefix, var, time):

    if var == 'fo' or var == 'fxX' or var == 'fxY' or var == 'fyX' or var == 'fyY':
        suffix = '.xyv'
    else:
        suffix = '.xy'

    if var == 'Z2niX' or var == 'Z2niY' or var == 'Z':
        fname = path + '/' + fprefix + '_' + var + suffix

    else:
        fname = path + '/' + fprefix + '_' + var + '_' + time + suffix
    return fname


#-----------------------------------------------------------------------
def index_time_to_string_time(time_index):
    time_str = '%.2i' % time_index
    return time_str


#-----------------------------------------------------------------------


def construct_label(fname, k_on=True):
    '''
        Construct leg label for a given filename
        final_lab, k_lab,h_lab,B_lab  = construct_label(fname)
    '''
    s = re.search('2D(?P<k_num>\d)', fname)
    if s:
        k_lab = s.group('k_num')
        kpre = 'k='
    else:
        k_lab = ' '
        kpre = ' '
        speckle = re.search('tc(?P<t_coherence>\d+)', fname)
        if speckle:
            kpre = '\tau_c  = '
            k_lab = speckle.group('t_coherence')

    b = re.search('noB', fname)
    if b:
        B_lab = 'no$ $B'
    else:
        B_lab = 'B'

    premag = re.search('_B(?P<premag>[0-9_]+)', fname)
    if premag:
        mag = premag.group('premag')[:-1]
        mag2 = re.sub('_', '.', mag)
        B_lab = r' \chi =' + mag2

    premag2 = re.search('_(?P<premag>\d+)T', fname)
    if premag2 and fname[:2] == 'r5':
        mag = premag2.group('premag')
        B_lab = mag + '\,T'

    hydro = re.search('static', fname)
    if hydro:
        h_lab = 'static'
    else:
        hydro_Cy = re.search('noCy', fname)
        if hydro_Cy:
            h_lab = 'no$ $C_y'
        else:
            h_lab = 'hydro'

    maxw = re.search('MaxW', fname)
    if maxw:
        m_lab = ' '
        B_lab = B_lab + '$ $MW'
    else:
        m_lab = ' '

    hi = re.search('hydro4', fname)
    if hi:
        m_lab = ' '
        B_lab = B_lab + ' +T_i'

    if k_on != True:
        final_lab = h_lab + ' + ' + B_lab
    else:
        final_lab = kpre + k_lab + '$ $' + h_lab + ' + ' + B_lab
    return final_lab, k_lab, h_lab, B_lab


#-----------------------------------------------------------------------


def get_tloop(path, var='Te'):
    '''
        t_index_l,t_col_l = get_tloop(path)
    '''
    list_files = os.listdir(path)
    t_index_l = []
    t_col_l = []
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            t_index_l.append(t_l)
            t_col = get_time(path + '/' + fname)
            t_col_l.append(t_col)
        return t_index_l, t_col_l


def search_path(regexp, path):
    '''
        out_name = search_path(regexp,path)
    '''
    out_name = ''
    for pp in os.listdir(path):
        rr = re.search(regexp, pp)
        if rr:
            print 'path = ', pp
            out_name = pp
            break
    return out_name


#-----------------------------------------------------------------------


def get_startline(fname):
    '''
        get the line where the data starts
    '''
    f = open(fname)
    f_list = f.readlines()
    out_line = 9
    for i in range(len(f_list)):
        if re.match(r'\s*\n', f_list[i]):
            out_line = i
            break

    return out_line


#-----------------------------------------------------------------------


def retrieve_path_tag(path):
    s = re.search('2D(?P<k_num>\d)', path)
    if s:
        k_lab = s.group('k_num')
    else:
        k_lab = ' '

    b = re.search('noB', path)
    if b:
        B_lab = 'noB'
    else:
        B_lab = 'B'

    hydro = re.search('static', path)
    if hydro:
        h_lab = 'static'
    else:
        h_lab = 'hydro'

    perc = re.search('_(?P<perc>\d+)pp', path)
    if perc:
        perc_lab = '_' + str(perc.group('perc'))
    else:
        perc_lab = ''

    final_lab = 'k' + k_lab + '_' + h_lab + '_' + B_lab + perc_lab
    return final_lab


#-----------------------------------------------------------------------
def interp_data(x_in, y_in, x_data_smooth):
    '''
        SI.PchipInterpolater
    '''
    f = SI.PchipInterpolator(x_in, y_in, extrapolate=True)
    y_data_smooth = f(x_data_smooth)    #linear

    return y_data_smooth


#-----------------------------------------------------------------------


def extract_power(x):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return b


#-----------------------------------------------------------------------


def fmt(x, pos=0.0):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


#-----------------------------------------------------------------------


def fmt_ord(x, pos):

    return r'${}$'.format(x)


#--- limiter function --------------------------------------------------
def limit_by(array, lim):
    '''
        lim can either be scalar or array the same shape as array
    '''
    array = np.where(np.abs(array) >= lim, np.sign(array) * lim, array)
    return array


#-----------------------------------------------------------------------


def trim_array(array, nx, ny):
    '''
        takes the central (nx,ny) slab of any array
        out_array = trim_array(array,nx,ny)
    '''

    if len(np.shape(array)) == 2:
        nx_a, ny_a = np.shape(array)
    else:
        nx_a = np.shape(array)
        if len(array) == (nx):
            return array
        else:
            out_array = np.zeros((nx))
            offset_x = (len(array[:]) % nx) / 2
            out_array = array[offset_x:-offset_x]
            return out_array

    if np.shape(array) == (nx, ny):
        return array
    else:
        out_array = np.zeros((nx, ny))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -nx, -ny

        if nx_a != nx:
            offset_x = ((len(array[:, 0]) % nx) / 2)
            offset_x_rh = (len(array[:, 0]) % nx) / 2
            offset_x += int(2.0 * np.abs((len(array[:, 0]) % nx) / 2.0 -
                                         (len(array[:, 0]) % nx) / 2))

        if ny_a != ny:
            offset_y = (len(array[0, :]) % ny) / 2
            offset_y_rh = (len(array[0, :]) % ny) / 2
            offset_y += int(2.0 * np.abs((len(array[0, :]) % ny) / 2.0 -
                                         (len(array[0, :]) % ny) / 2))

        #--------------

        if ny == 1:
            offset_y = 1
            offset_y_rh = 1
        if offset_x_rh != 0 and offset_y_rh != 0:
            out_array = array[offset_x:-offset_x_rh, offset_y:-offset_y_rh]
        elif offset_x_rh != 0 and offset_y_rh == 0:
            out_array = array[offset_x:-offset_x_rh, offset_y:]
        elif offset_x_rh == 0 and offset_y_rh != 0:
            out_array = array[offset_x:, offset_y:-offset_y_rh]
        else:
            out_array = array[offset_x:, offset_y:]
        if np.shape(out_array) != (nx, ny):
            print ' ERROR ARRAY NOT TRIMMED PROPERLY '
            print ' SHAPE ARRAY IN = ', np.shape(array)
            print ' SHAPE ARRAY OUT = ', np.shape(out_array), ' shape asked for = (%i,%i)' % (nx,
                                                                                              ny)
            sys.exit()
        return out_array


#-----------------------------------------------------------------------


def trim_array_1D(array, nx, ny):
    '''
        takes the central (nx,ny) slab of any array
        out_array = trim_array(array,nx,ny)
    '''
    if ny == 0 and len(np.shape(array)) == 1:
        out_array = np.zeros((nx))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -nx, -ny
        nx_a = len(array)
        if nx_a != nx:
            offset_x = ((len(array) % nx) / 2)
            offset_x_rh = (len(array) % nx) / 2
            offset_x += int(2.0 * np.abs((len(array) % nx) / 2.0 - (len(array) % nx) / 2))

        if offset_x_rh != 0:
            out_array = array[offset_x:-offset_x_rh]
        elif offset_x_rh == 0:
            out_array = array[offset_x:]
        else:
            out_array = array[offset_x:]
        if len(out_array) != nx:
            print ' ERROR ARRAY NOT TRIMMED PROPERLY - trim 1D'
            print ' SHAPE ARRAY IN = ', np.shape(array)
            print ' SHAPE ARRAY OUT = ', np.shape(out_array), ' shape asked for = (%i,%i)' % (nx,
                                                                                              ny)
            sys.exit()
            return array
        #print ' np.shape(out_array) = ', np.shape(out_array)
        return out_array
    #print('---> shape in', np.shape(array))
    # 2D stuff --->
    #======================================================?
    if ny == 3 and np.shape(array)[1] == 1:
        #print(' shape = ', np.shape(array))
        if len(np.shape(array)) == 3:
            nv, ny_in, nx_in = np.shape(array)
            out_array = np.zeros((nv, ny, nx))

        else:
            ny_in, nx_in = np.shape(array)

            if nx_in == 1:
                arr_temp = array[:, 0]
            else:
                arr_temp = array
            out_array = np.zeros((nx, ny))
            for iy in range(ny):
                out_array[:, iy] = arr_temp
            return out_array

    if len(np.shape(array)) == 2:
        nx_a, ny_a = np.shape(array)
    else:
        nx_a = np.shape(array)
        if len(array) == (nx):
            return array
        else:
            out_array = np.zeros((nx))
            offset_x = (len(array[:]) % nx) / 2
            out_array = array[offset_x:-offset_x]
            return out_array

    if np.shape(array) == (nx, ny):
        return array
    else:
        out_array = np.zeros((nx, ny))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -nx, -ny

        if nx_a != nx:
            offset_x = ((len(array[:, 0]) % nx) / 2)
            offset_x_rh = (len(array[:, 0]) % nx) / 2
            offset_x += int(2.0 * np.abs((len(array[:, 0]) % nx) / 2.0 -
                                         (len(array[:, 0]) % nx) / 2))
        if ny_a != ny:
            offset_y = (len(array[0, :]) % ny) / 2
            offset_y_rh = (len(array[0, :]) % ny) / 2
            offset_y += int(2.0 * np.abs((len(array[0, :]) % ny) / 2.0 -
                                         (len(array[0, :]) % ny) / 2))

        if ny == 1:
            offset_y = 1
            offset_y_rh = 1
        #print ' nx, ny, nx_a, ny_a = ', nx, ny, nx_a,ny_a
        #print ' offset_x = ', offset_x, offset_x_rh, offset_y, offset_y_rh
        if offset_x_rh != 0 and offset_y_rh != 0:
            out_array = array[offset_x:-offset_x_rh, offset_y:-offset_y_rh]
        elif offset_x_rh != 0 and offset_y_rh == 0:
            out_array = array[offset_x:-offset_x_rh, offset_y:]
        elif offset_x_rh == 0 and offset_y_rh != 0:
            out_array = array[offset_x:, offset_y:-offset_y_rh]
        else:
            out_array = array[offset_x:, offset_y:]
        if np.shape(out_array) != (nx, ny):
            print ' ERROR ARRAY NOT TRIMMED PROPERLY - trim 1D'
            print ' SHAPE ARRAY IN = ', np.shape(array)
            print ' SHAPE ARRAY OUT = ', np.shape(out_array), ' shape asked for = (%i,%i)' % (nx,
                                                                                              ny)
            sys.exit()
            return array
        #print ' np.shape(out_array) = ', np.shape(out_array)
        #sys.exit()
        return out_array


def get_l_num_lmfp(fname):
    '''
        l_num = get_l_num_lmfp(fname)
    '''
    s = re.search('(\d+)lmfp', fname)
    if s:
        l_num = s.group()[:-4]
    else:
        l_num = 40
    return float(l_num)


#-----------------------------------------------------------------------
def get_it_list(path):
    '''                                                                                                                                 
    last_dump_index,last_time =  get_last_dump(path)                                                                                    
    '''
    list_files = os.listdir(path)
    t_list = []
    tc_list = []
    for fname in list_files:
        var = 'Te'
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            if not int(t_l) in t_list:
                t_list.append(int(t_l))
    t_list = np.sort(t_list)
    print t_list
    return t_list


#-----------------------------------------------------------------------
def get_last_dump(path):
    '''                                                                                                                                 
    last_dump_index,last_time =  get_last_dump(path)                                                                                    
    '''
    list_files = os.listdir(path)
    t_list = []
    tc_list = []
    for fname in list_files:
        var = 'Te'
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            t_list.append(int(t_l))
            #t_col = get_time(path + '/' + fname)
            #tc_list.append(float(t_col))
            print t_list
            #fprefix = s.group('fpre')
    last_dump_index = np.max(t_list)
    return last_dump_index


#-----------------------------------------------------------------------
def get_tau_n(ne0, Te0):
    '''
        tau_n = get_tau_n(ne0,Te0)
    '''
    taun = 0.75 * (np.pi**0.5)    # time controls
    tau_n = taun * ((2.0 * Te0)**1.5) / ne0
    return tau_n


#-----------------------------------------------------------------------
def fpg_get_info(fname):
    '''
        Gets the IMPACT header info
        time,ndims,dim_array = fpg_get_info(fname)
        dim_Aray = [nv,ny,nx] or [ny,nx]
    '''

    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    mat = np.loadtxt(fname, skiprows=10)

    #print mat_info[0]
    time = float(mat_info[1])
    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])

        nx = int(mat_info[5])
        x_grid = np.array(mat_info[6])
        dim_array = [ny, nx]
        #put this in a dictionary?
        #grid_array = [
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])
        dim_array = [nv, ny, nx]
    return time, ndims, dim_array


#-----------------------------------------------------------------------
def get_fprefix(path, var='Te'):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
        fprefix = get_fprefix(path,var='Te')
    '''
    t_list = []
    tc_list = []
    list_files = os.listdir(path)
    # print ' path = ', path
    if path[-1] == '/':
        path = path[:-1]
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            t_list.append(t_l)
            t_col = get_time(path + '/' + fname)
            tc_list.append(t_col)
            fprefix = s.group('fpre')
            break
    return fprefix


#-----------------------------------------------------------------------
def get_t_list(path, var='Te'):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
        t_list, tc_list = get_t_list(path,var='Te')
    '''
    t_list = []
    tc_list = []
    list_files = os.listdir(path)
    # print ' path = ', path
    if path[-1] == '/':
        path = path[:-1]
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            t_list.append(t_l)
            t_col = get_time(path + '/' + fname)
            tc_list.append(t_col)
            fprefix = s.group('fpre')
    return t_list, tc_list


#-----------------------------------------------------------------------
def load_path(path, var):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''

    #dict = MyDict()
    #mat = np.loadtxt(fname,skiprows=out_l)
    t_list = []
    tc_list = []
    list_files = os.listdir(path)
    print ' path = ', path
    if path[-1] == '/':
        path = path[:-1]
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var + '_(?P<time>\d\d).xy.*?', fname)
        if s:
            t_l = s.group('time')
            t_list.append(t_l)
            t_col = get_time(path + '/' + fname)
            tc_list.append(t_col)
            print t_list
            fprefix = s.group('fpre')

    #--- get grids --->
    nt = len(t_list)
    time = t_list[0]
    fname = construct_fname(path, fprefix, var, time)
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    out_l = get_startline(fname)
    mat = np.loadtxt(fname, skiprows=out_l)
    dict = {}
    if var == 'Cx' or var == 'Cy':
        mat = np.loadtxt(fname, skiprows=out_l)
        dict['mat'] = mat
        return dict

    ##print mat_info
    time = float(mat_info[1])

    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])
        mat = np.loadtxt(fname, skiprows=out_l)
        grid = [y_grid, x_grid]
        hmat_final = np.transpose(mat)
        # may also need to multiply by -1 if a vector quantity (
        # we are flipping x axis)?
        dim_array = np.array([ny, nx])
        if var[-1] == 'z' or var[-1] == 'Y' or var[-1] == 'X':
            print '---- warning this routine should not be used with vector quantities---'
            print ' we have performed a transpose, may need to multiply by -1'

        v_grid = 0.0
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])

        mat = np.loadtxt(fname, skiprows=out_l)
        dim_array = np.array([nv, ny, nx])
        dim_reverse = [nx, ny, nv]

        mat = mat.reshape(dim_reverse)
        hmat_final = np.transpose(mat)
        grid = [v_grid, y_grid, x_grid]    #nv,ny,nx
        dict['v_grid'] = v_grid
        dict['nv'] = nv
    #-------------------------
    #if ndims==2:
    if ndims == 2:
        hmat_t = np.zeros((nt, ny, nx))
    else:
        hmat_t = np.zeros((nt, nv, ny, nx))
    for time in t_list:
        fname = construct_fname(path, fprefix, var, time)
        if ndims == 2:
            mat = np.loadtxt(fname, skiprows=out_l)
            print ' t_list index = ', t_list.index(time), np.shape(mat), np.shape(hmat_t)
            hmat_t[t_list.index(time), :, :] = np.transpose(mat)

        else:
            mat = mat.reshape(dim_reverse)
            hmat_final = np.transpose(mat)
            hmat_t[t_list.index(time), :, :, :] = mat
    dict['time'] = tc_list
    dict['y_grid'] = y_grid
    dict['x_grid'] = x_grid
    dict['nx'] = nx
    dict['ny'] = ny
    dict['mat'] = hmat_t

    return dict


#-----------------------------------------------------------------------
def load_dict(path, fprefix, var, time):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''
    fname = construct_fname(path, fprefix, var, time)

    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    out_l = get_startline(fname)
    mat = np.loadtxt(fname, skiprows=out_l)
    dict = {}

    time = float(mat_info[1])

    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])

        mat = np.loadtxt(fname, skiprows=out_l)
        grid = [y_grid, x_grid]
        hmat_final = mat
        v_grid = 0.0
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])

        mat = np.loadtxt(fname, skiprows=out_l)
        dim_array = [nv, ny, nx]
        dim_reverse = [nx, ny, nv]

        mat = mat.reshape(dim_reverse)
        hmat_final = np.transpose(mat)
        grid = [v_grid, y_grid, x_grid]    #nv,ny,nx

        dict['v_grid'] = v_grid
        dict['nv'] = nv
    dict['time'] = time
    dict['y_grid'] = y_grid
    dict['x_grid'] = x_grid
    dict['nx'] = nx
    dict['ny'] = ny
    dict['mat'] = hmat_final

    return dict


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def load_dict_1D(path, fprefix, var, time):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
        
        This is a special loading routine for 1D analysis which ensures all arrays are atleast 2D for use 
        in kinetic/classical heat flow reconstructions.
    '''
    fname = construct_fname(path, fprefix, var, time)

    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    out_l = get_startline(fname)
    mat = np.loadtxt(fname, skiprows=out_l)
    dict = {}
    time = float(mat_info[1])

    ndims = int(mat_info[2])
    mat = np.loadtxt(fname, skiprows=out_l)
    if (len(np.shape(mat))) == 1:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        #y_grid = np.array(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])

        mat = np.loadtxt(fname, skiprows=out_l)
        grid = [y_grid, x_grid]
        hmat_final = np.zeros((len(mat), 3))
        for ii in range(3):
            hmat_final[:, ii] = mat
        v_grid = 0.0

    elif ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])

        mat = np.loadtxt(fname, skiprows=out_l)

        grid = [y_grid, x_grid]
        hmat_final = mat
        v_grid = 0.0
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])

        mat = np.loadtxt(fname, skiprows=out_l)
        dim_array = [nv, ny, nx]
        dim_reverse = [nx, ny, nv]

        mat = mat.reshape(dim_reverse)
        hmat_3 = np.transpose(mat)

        if ny == 1:
            hmat_final = np.zeros((nv, 3, nx))
            for iy in range(3):
                hmat_final[:, iy, :] = hmat_3[:, 0, :]
            ny = 3
        else:
            hmat_final = hmat_3

        grid = [v_grid, y_grid, x_grid]    #nv,ny,nx

        dict['v_grid'] = v_grid
        dict['nv'] = nv
    dict['time'] = time
    dict['y_grid'] = y_grid
    dict['x_grid'] = x_grid
    dict['nx'] = nx
    dict['ny'] = ny
    dict['mat'] = hmat_final

    return dict


#-----------------------------------------------------------------------


def load_dict_OLD(path, fprefix, var, time):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''
    fname = construct_fname(path, fprefix, var, time)
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    out_l = get_startline(fname)

    dict = MyDict()

    mat = np.loadtxt(fname, skiprows=out_l)
    dict['mat'] = mat

    if var == 'Cx' or var == 'Cy':
        return dict
    #print mat_info[0]
    time = float(mat_info[1])
    ndims = int(mat_info[2])

    dict['time'] = time
    dict['ndims'] = ndims
    if ndims == 2:

        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])
        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])
        dim_array = [ny, nx]

    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])
        dict['nv'] = nv
        dict['v_grid'] = v_grid
        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])
        dim_array = [nv, ny, nx]

    dict['nx'] = nx
    dict['ny'] = ny
    dict['x_grid'] = x_grid
    dict['y_grid'] = y_grid

    return dict


#-----------------------------------------------------------------------
def fpg_load(fname):
    '''
        loads the file
        returns hmat,v_grid,y_grid,x_grid = fpg_load(fname)
        now can deal with fo when halo cells are also saved
    '''
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    mat = np.loadtxt(fname, skiprows=10)

    ##print mat_info
    time = float(mat_info[1])
    ndims = int(mat_info[2])
    if ndims == 2:
        ny = int(mat_info[3])
        y_grid = list_to_float(mat_info[4])

        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])

        mat = np.loadtxt(fname, skiprows=8)
        grid = [y_grid, x_grid]
        hmat_final = mat
        v_grid = 0.0
    else:

        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])

        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])

        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])

        mat = np.loadtxt(fname, skiprows=10)
        dim_array = [nv, ny, nx]
        dim_reverse = [nx, ny, nv]

        mat = mat.reshape(dim_reverse)
        hmat_final = np.transpose(mat)
        grid = [v_grid, y_grid, x_grid]    #nv,ny,nx

    return hmat_final, v_grid, y_grid, x_grid


#-----------------------------------------------------------------------


def get_index(ny, nx):
    index = ny + 2 + 1
    return index


#-----------------------------------------------------------------------


def get_time(fname):
    f = open(fname)
    data = f.readlines()
    t = float(data[1])
    return t


#-----------------------------------------------------------------------


def get_hallparam(Bz_norm, vte_norm, Z2ni_norm):
    chi = Bz_norm * (vte_norm**3) * (0.75 * (np.pi**0.5)) * (1.0 / Z2ni_norm)
    return chi


def calc_norms(var, sample=0.0, forced_power=[], normal_class=conv_factors_custom):
    '''
        norm_const, ylab = calc_norms(var)
    '''
    ##print var, np.shape(var)
    c_fmt = '%3.2f'
    q_e = constants.q_e
    m_e = constants.m_e
    m_p = constants.m_p
    k_b = constants.k_b
    epsilon0 = constants.epsilon0

    c_index = normal_class.c_index
    SI_on = normal_class.SI_on
    v_te = normal_class.v_te
    tau_ei = normal_class.tau_ei
    nu_ei = normal_class.nu_ei
    lambda_mfp = normal_class.lambda_mfp

    n0, T0 = normal_class.n0, normal_class.T0
    #-----------------

    if SI_on:
        if var == 'Cx' or var == 'Cy':

            norm_const = v_te * 1e-3
            title = r'$' + var[0] + '_' + var[1] + r'$ [ $  kms^{-1}$ ]'

        elif var == 'n':
            power = extract_power(n0)
            mod = r'$ 10^{' + str(power) + '} $'

            norm_const = n0 * (10**-power)
            title = r'$n_e$ [' + mod + r'$\si{cm^{-3}}$' + r']'

        elif var == 'Te':
            norm_const = 2.0 * T0
            title = r'$T_e$ [ $eV$ ]'
            c_fmt = '%3.0f'
        elif var == 'Ui':
            norm_const = (2.0 / 3.0) * (2.0 * T0)
            title = r'$T_i$ [ $eV$ ]'
            c_fmt = '%3.0f'

        elif var == 'Bz':
            norm_const = (m_e / (q_e * tau_ei))

            power = extract_power(sample * norm_const)
            norm_const *= (10**-power)
            ##print ' sample = ', sample, 'power = ', power

            var_name = '$B_z$'
            if power == 0:
                mod = ''
                units = r'[$si{T}$]'
                title = r'$B_z$ [$\si{T}$ ]'

            else:
                mod = r'$ 10^{' + str(power) + '} $'

                units = r'[' + mod + '$\si{T}$]'
                title = r'$B_z$ [' + mod + '$\si{T}$ ]'
            c_fmt = '%1.1f'

        elif var == 'wt':
            if len(forced_power) != 0:
                power = forced_power[0]
            else:
                power = extract_power(sample)

            if power != 0:
                mod = r'$ 10^{' + str(power) + '} $'
            else:
                mod = r''
            ##mod = r'$ 10^{' +str(power)  + '} $'

            norm_const = 1.0 * (10**(-power))

            title = mod + r'$\omega \tau_{ei}$'
            c_fmt = '%1.1f'

        elif var[0] == 'E':
            #print 'DOING E-field - min sample = ', sample
            norm_const = (lambda_mfp / (tau_ei**2)) * (m_e / q_e)
            power = extract_power(norm_const * sample)
            #print ' power extracted = ', power
            c_fmt = '%1.1f'
            mod = r'$ 10^{' + str(power) + '}$'
            norm_const = norm_const * (10**-power)
            if power == 0:
                mod = ''
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$V/m$ ]'

        elif var[0] == 'q':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            if power == 0:
                mod = ''
            #c_fmt = '%1.0e'
            #norm_const = m_e*(v_te**3)*n_04
            title = r'$' + var[0] + '_' + var[
                1] + r'$ [ ' + mod + '$q_0$ ]'    # + r' [ $Jms^{-1}$ ]'

        elif var[0] == 'j':
            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            if power == 0:
                mod = ''
            #norm_const = 1.0#q_e*n0*v_te

            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$j_0$ ]'    #r' [ $C/m^2 s$ ]'
        elif var == 'U':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            #c_fmt = '%1.0e'
            #norm_const = m_e*(v_te**3)*n_04
            title = r'$' + var[0] + r'$ [ ' + mod + '$m_e v_n^2 n_0$ ]'    # + r' [ $Jms^{-1}$ ]'

    else:
        norm_const = 1.0
        title = var

        if var == 'Cx' or var == 'Cy':

            norm_const = 1.0
            title = r'$' + var[0] + '_' + var[1] + r'$ [ $  v_{th}$ ]'

        elif var == 'n':
            norm_const = n0 * 1e-22
            title = r'$n_e$ [ $n_0$ ]'

        elif var == 'Te':
            norm_const = 2.0
            title = r'$T_e$ [ $T_{e,0}$ ]'
            c_fmt = '%3.0f'

        elif var == 'Bz':
            norm_const = 1.0

            power = extract_power(sample * norm_const)
            norm_const *= (10**-power)
            var_name = '$B_z$'
            if power == 0:
                mod = ''
                units = r'[$si{T}$]'
                title = r'$B_z$ [$B_{z,n}$ ]'

            else:
                mod = r'$ 10^{' + str(power) + '} $'

                units = r'[' + mod + '$B_{z,n}$]'
                title = r'$B_z$ [' + mod + '$B_{z,n}$ ]'
            c_fmt = '%1.1f'

        elif var == 'wt':
            if len(forced_power) != 0:
                power = forced_power[0]
            else:
                power = extract_power(sample)
            if power != 0:
                mod = r'$ 10^{' + str(power) + '} $'
            else:
                mod = r''
            norm_const = 1.0 * (10**(-power))

            title = mod + r'$\omega \tau_{ei}$'
            c_fmt = '%1.1f'

        elif var[0] == 'E':
            #print 'DOING E-field - min sample = ', sample
            norm_const = (lambda_mfp / (tau_ei**2)) * (m_e / q_e)
            power = extract_power(norm_const * sample)
            c_fmt = '%1.1f'
            mod = r'$ 10^{' + str(power) + '}$'
            norm_const = norm_const * (10**-power)

            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$V/m$ ]'

        elif var[0] == 'q':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'

            title = r'$' + var[0] + '_' + var[
                1] + r'$ [ ' + mod + '$q_0$ ]'    # + r' [ $Jms^{-1}$ ]'

        elif var[0] == 'j':
            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'

            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$j_0$ ]'    #r' [ $C/m^2 s$ ]'
        elif var[0] == 'U':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            title = r'$' + var[0] + r'$ [ ' + mod + '$m_e v_n^2 n_0$ ]'    # + r' [ $Jms^{-1}$ ]'

    #print 'var: ', var, title, c_fmt
    return norm_const, title, c_fmt


def calc_norms_2(var, sample=0.0, SI_on=True, normal_class=conv_factors_custom(), forced_power=[]):
    '''
        dict = calc_norms(var)
        
    '''

    q_e = constants.q_e
    m_e = constants.m_e
    m_p = constants.m_p
    k_b = constants.k_b
    epsilon0 = constants.epsilon0
    v_te = normal_class.v_te
    c_index = normal_class.c_index
    SI_on = normal_class.SI_on
    tau_ei = normal_class.tau_ei
    nu_ei = normal_class.nu_ei
    lambda_mfp = normal_class.lambda_mfp
    xstep_factor = normal_class.xstep_factor
    tstep_factor = normal_class.tstep_factor
    xlab = normal_class.xlab
    ylab = normal_class.ylab
    n0, T0 = normal_class.n0, normal_class.T0
    Z0 = normal_class.Z
    c_fmt = '%3.2f'

    def get_mod(forced_power, sample_in=0.0):

        if len(forced_power) != 0:
            power = forced_power[0]
        else:
            power = extract_power(sample_in)
        power = int(power)

        if power != 0:
            mod = r'$ 10^{' + str(power) + '} $'
            mult = 10**-power
        else:
            mod = r''
            mult = 1.0
        return mod, mult

    if SI_on:
        if var == 'Cx' or var == 'Cy':

            norm_const = v_te * 1e-3
            title = r'$' + var[0] + '_' + var[1] + r'$ [ $  kms^{-1}$ ]'
            var_name = '$' + var[0] + '_' + var[1] + r'$ '
            units = r'$ [ $  kms^{-1}$ ]'
        elif var == 'n':
            norm_const = n0 * 1e-22
            var_name = r'$n_e$'
            units = r'[ $10^{22}$ $cm^-3$ ]'
            title = r'$n_e$ [ $10^{22}$ $cm^-3$ ]'

        elif var == 'Te':
            norm_const = 2.0 * T0
            title = r'$T_e$ [ $\si{\eV}$ ]'
            units = r'[ $\si{\eV}$ ]'
            var_name = r'$T_e$'
            c_fmt = '%3.0f'
        elif var == 'Te keV':
            norm_const = 2.0 * T0 * 1e-3
            title = r'$T_e$ [ $\si{\kilo\eV}$ ]'
            units = r'[ $\si{\kilo\eV}$ ]'
            var_name = r'$T_e$'
            c_fmt = '%3.0f'
        elif var == 'Z':
            norm_const = Z0
            var_name = r'Z'
            c_fmat = '%3.0f'
            units = r''
            title = 'Z'
        elif var == 'Bz':
            norm_const = (m_e / (q_e * tau_ei))    # prof_Bz_ave**-1

            #power = extract_power(sample*norm_const)
            power = 0
            norm_const *= (10**-power)
            print ' sample = ', sample, 'power = ', power

            var_name = '$B_z$'
            if power == 0:
                mod = ''
                units = r'[$\si{T}$]'
                title = r'$B_z$ [$\si{T}$ ]'

            else:
                mod = r'$ 10^{' + str(power) + '} $'

                units = r'[' + mod + '$\si{T}$]'
                title = r'$B_z$ [' + mod + '$\si{T}$ ]'
            c_fmt = '%1.1f'

        elif var == 'wt':
            if len(forced_power) != 0:
                power = forced_power[0]
            else:
                power = extract_power(sample)

            if power != 0:
                mod = r'$ 10^{' + str(power) + '} $'
            else:
                mod = r''
            ##mod = r'$ 10^{' +str(power)  + '} $'

            norm_const = 1.0 * (10**(-power))
            var_name = r'$\omega \tau_{ei}$'
            units = ''
            title = mod + r'$\omega \tau_{ei}$'
            c_fmt = '%1.1f'

        elif var[0] == 'E':
            #print 'DOING E-field - min sample = ', sample
            norm_const = (lambda_mfp / (tau_ei**2)) * (m_e / q_e)
            power = extract_power(norm_const * sample)
            #print ' power extracted = ', power
            c_fmt = '%1.1f'
            mod = r'$ 10^{' + str(power) + '}$'
            norm_const = norm_const * (10**-power)
            units = r'$ [ ' + mod + '$V/m$ ]'
            var_name = r'$' + var[0] + '_' + var[1]
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$V/m$ ]'

        elif var[0] == 'q':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            #c_fmt = '%1.0e'
            #norm_const = m_e*(v_te**3)*n_04
            var_name = r'$' + var[0] + '_' + var[1]
            units = r'$ [ ' + mod + '$q_0$ ]'
            title = r'$' + var[0] + '_' + var[
                1] + r'$ [ ' + mod + '$q_0$ ]'    # + r' [ $Jms^{-1}$ ]'

        elif var[0] == 'P':

            c_fmt = '%1.1f'
            norm_const = (2.0 / 3.0) * m_e * (v_te**2) * (n0 * 1e6) * 1e-5

            mod, mult = get_mod(forced_power, sample * norm_const)
            norm_const *= mult

            var_name = r'$' + var[0] + '$'
            units = r'[ ' + mod + '$bar$ ]'
            title = r'$' + var[0] + r'$ [ ' + mod + 'bar ]'    # + r' [ $Jms^{-1}$ ]'

        elif var[0] == 'U':

            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            #c_fmt = '%1.0e'
            #norm_const = m_e*(v_te**3)*n_04
            var_name = r'$' + var[0] + '$'
            units = r'[ ' + mod + '$m_e v_n^2 n_0$ ]'
            title = r'$' + var[0] + r'$ [ ' + mod + '$m_e v_n^2 n_0$ ]'    # + r' [ $Jms^{-1}$ ]'

        elif var[0] == 'j':
            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0 * (10**-power)
            mod = r'$ 10^{' + str(power) + '} $'
            #norm_const = 1.0#q_e*n0*v_te
            var_name = r'$' + var[0] + '_' + var[1] + r'$'
            units = r'$ [ ' + mod + '$j_0$ ]'
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$j_0$ ]'    #r' [ $C/m^2 s$ ]'
    else:
        norm_const = 1.0
        title = var

    dict = {}
    dict['norm_const'] = norm_const
    dict['title'] = title
    dict['c_fmt'] = c_fmt
    dict['var_name'] = var_name
    dict['units'] = units
    #print 'var: ', var, title
    return dict


#-----------------------------------------------------------------------
def get_grad_1d(x_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata = get_grad_1d(x_grid,data)
    '''
    dxT = np.zeros((len(T_data)))
    dx = x_grid[2:] - x_grid[:-2]
    dxT[1:-1] = (T_data[2:] - T_data[:-2]) / dx
    dxT[0] = dxT[1]
    dxT[-1] = dxT[-2]
    return dxT


def get_grad_dx2_1d(x_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata = get_grad_1d(x_grid,data)
    '''
    dxT = np.zeros((len(T_data)))
    dx = x_grid[2:] - x_grid[:-2]
    dxT[1:-1] = (T_data[2:] - 2.0 * T_data[1:-1] + T_data[:-2]) / (2.0 * (dx**2))
    dxT[0] = dxT[1]
    dxT[-1] = dxT[-2]
    return dxT


def add_grid_bcs(grid_in):
    grid = np.zeros((len(grid_in) + 2))
    grid[1:-1] = grid_in
    grid[0] = grid_in[0] - (grid_in[1] - grid_in[0])
    grid[-1] = grid_in[-1] + (grid_in[-1] - grid_in[-2])
    return grid


#-----------------------------------------------------------------------
def get_grad(x_grid, y_grid, T_data, bc='none'):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata,dydata = get_grad(x_grid,y_grid,data)
    '''
    if bc == 'periodic':
        nx, ny = np.shape(T_data)
        T_data_glb = np.zeros((nx + 2, ny + 2))

        T_data_glb[1:-1, 1:-1] = T_data.copy()

        T_data_glb[0, 1:-1] = T_data[-1, :]
        T_data_glb[-1, 1:-1] = T_data[0, :]

        T_data_glb[1:-1, 0] = T_data[:, -1]
        T_data_glb[1:-1, -1] = T_data[:, 0]
        x_grid_glb = add_grid_bcs(x_grid)
        y_grid_glb = add_grid_bcs(y_grid)

        X, Y = np.meshgrid(y_grid_glb, x_grid_glb)
        dy = Y[2:, :] - Y[:-2, :]
        dx = X[:, 2:] - X[:, :-2]
        dxT = (T_data_glb[:, 2:] - T_data_glb[:, :-2]) / dx
        dyT = (T_data_glb[2:, :] - T_data_glb[:-2, :]) / dy
        return dyT, dxT
    else:
        if len(np.shape(T_data)) == 1:
            dx = x_grid[2:] - x_grid[:-2]
            dxT = (T_data[2:] - T_data[:-2]) / dx
            return dxT
        else:
            dx = x_grid[2:] - x_grid[:-2]
            dy = y_grid[2:] - y_grid[:-2]
            dxT = np.zeros((len(dx), len(dy)))
            dyT = np.zeros((len(dx), len(dy)))
            ldx, ldy = len(dx), len(dy)
            #print 'SHAPE DX DY = ', ldx,ldy
            for ii in range(len(dx)):
                #print 'ii=  ', ii
                dyT[ii, :] = (T_data[ii, 2:] - T_data[ii, :-2]) / dy

            for ii in range(len(dy)):
                #print 'ix = ', ii, np.shape((T_data[2:,ii] - T_data[:-2,ii])), np.shape(dxT), np.shape(dx)
                dxT[:, ii] = (T_data[2:, ii] - T_data[:-2, ii]) / dx
            return dxT, dyT


#-----------------------------------------------------------------------


def get_grad_varZ(x_grid, y_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata,dydata = get_grad(x_grid,y_grid,data)
    '''

    if len(np.shape(T_data)) == 1:
        dx = x_grid[2:] - x_grid[:-2]
        dxT = (T_data[2:] - T_data[:-2]) / dx
        return dxT
    else:
        nx, ny = np.shape(T_data)
        if len(x_grid) > len(y_grid):

            dx = x_grid[2:] - x_grid[:-2]

            dxT = np.zeros((len(dx), ny))
            dyT = np.zeros((len(dx), ny))
            for ii in range(ny):
                #print 'ix = ', ii, np.shape((T_data[2:,ii] - T_data[:-2,ii])), np.shape(dxT), np.shape(dx)
                dxT[:, ii] = (T_data[2:, ii] - T_data[:-2, ii]) / dx
            return dxT, dyT
        else:
            dy = y_grid[2:] - y_grid[:-2]
            dxT = np.zeros((nx, len(dy)))
            dyT = np.zeros((nx, len(dy)))
            for ii in range(nx):
                dyT[ii, :] = (T_data[ii, 2:] - T_data[ii, :-2]) / dy
            return dxT, dyT


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------


def load_data_all(path, fprefix, time):
    '''
        dict  = load_data_all_1D(path,fprefix,time)
    '''

    dict_Bz = load_dict(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = load_dict(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = load_dict(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']

    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = load_dict(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = load_dict(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    wt = get_wt(path, time)
    dict_ne = load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']

    dict_Cx = load_dict(path, fprefix, 'Cx', time)
    Cx = dict_Cx['mat']

    ##Z2ni = ne

    dict_te = load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']

    grid = dict_fo
    nv, ny, nx = np.shape(fo)

    grid['x_grid'] = trim_array_1D(x_grid, nx, 0)
    grid['y_grid'] = trim_array_1D(y_grid, ny, 0)
    v_grid = grid['v_grid']

    dict_Z = load_dict(path, fprefix, 'Z', time)
    prof_Z = dict_Z['mat']
    Z2ni = np.transpose(trim_array(ne, nx, ny)) * np.transpose(trim_array(prof_Z, nx, ny))
    dict_U = load_dict(path, fprefix, 'U', time)
    Ue = dict_U['mat']

    dxB, dyB = get_grad(x_grid, y_grid, Bz)

    Bz = np.transpose(trim_array(Bz, nx, ny))
    jx = np.transpose(trim_array(jx_c, nx, ny))
    jy = np.transpose(trim_array(jy_c, nx, ny))
    qx = np.transpose(trim_array(qx_c, nx, ny))
    qy = np.transpose(trim_array(qy_c, nx, ny))
    wt = np.transpose(trim_array(wt, nx, ny))
    prof_Z = np.transpose(trim_array(prof_Z, nx, ny))
    Ue = np.transpose(trim_array(Ue, nx, ny))

    #-------
    #print('---> loading all data', nx,ny)
    #pdb.set_trace()
    dxT, dyT = get_grad(x_grid, y_grid, Te)
    ne = np.transpose(trim_array(ne, nx, ny))
    Cx = np.transpose(trim_array(Cx, nx, ny))

    Te = np.transpose(trim_array(Te, nx, ny))
    dxT = np.transpose(trim_array(dxT, nx, ny))
    dyT = np.transpose(trim_array(dyT, nx, ny))
    dxB = np.transpose(trim_array(dxB, nx, ny))
    dyB = np.transpose(trim_array(dyB, nx, ny))

    rA = (Z2ni)**-1

    dict = {}

    dict['ne'] = ne
    dict['Te'] = Te
    dict['U'] = Ue

    dict['Cx'] = Cx
    dict['qx'] = qx
    dict['qy'] = qy
    dict['jx'] = jx
    dict['jy'] = jy
    dict['dxT'] = dxT
    dict['dyT'] = dyT
    dict['Bz'] = Bz
    dict['dxB'] = dxB
    dict['dyB'] = dyB
    dict['x_grid'] = x_grid
    dict['y_grid'] = y_grid
    dict['v_grid'] = v_grid

    dict['grid'] = grid
    dict['Z2ni'] = Z2ni
    dict['wt'] = wt
    dict['Z'] = prof_Z
    dict['fo'] = fo
    dict['nv'] = nv
    dict['ny'] = ny
    dict['nx'] = nx
    dict['time'] = dict_te['time']

    return dict


#-----------------------------------------------------------------------
def load_mat_1D(path_in, fprefix, time, var, nx, ny):
    dict_te = load_dict_1D(path_in, fprefix, var, time)
    mat = dict_te['mat']
    return np.transpose(trim_array_1D(mat, nx, ny))


def load_data_all_1D(path, fprefix, time):
    '''
        dict  = load_data_all_1D(path,fprefix,time)
    '''
    dict_Bz = load_dict_1D(path, fprefix, 'Bz', time)
    Bz = dict_Bz['mat']

    dict_jxX = load_dict_1D(path, fprefix, 'jxX', time)
    jxX = dict_jxX['mat']
    dict_jyY = load_dict_1D(path, fprefix, 'jyY', time)
    jyY = dict_jyY['mat']

    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])

    dict_qxX = load_dict_1D(path, fprefix, 'qxX', time)
    qxX = dict_qxX['mat']
    dict_qyY = load_dict_1D(path, fprefix, 'qyY', time)
    qyY = dict_qyY['mat']
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])

    wt = get_wt_1D(path, time)

    dict_ne = load_dict_1D(path, fprefix, 'n', time)
    ne = dict_ne['mat']

    ##Z2ni = ne

    dict_te = load_dict_1D(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = load_dict_1D(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    #print '-------------------------------------------'
    #print ' \n\ndict ----- time ===== ', dict_fo['time']
    #print '-------------------------------------------'
    grid = dict_fo
    nv, ny, nx = np.shape(fo)
    ny = 3
    grid['x_grid'] = trim_array_1D(x_grid, nx, 0)
    grid['y_grid'] = trim_array_1D(y_grid, ny, 0)
    grid['ny'] = 3
    #nv,ny,nx = np.shape(fo)
    dict_Z = load_dict_1D(path, fprefix, 'Z', time)
    prof_Z = dict_Z['mat']
    Z2ni = np.transpose(trim_array_1D(ne, nx, ny)) * np.transpose(trim_array_1D(prof_Z, nx, ny))
    dict_U = load_dict_1D(path, fprefix, 'U', time)
    Ue = dict_U['mat']

    jx = np.transpose(trim_array_1D(jx_c, nx, ny))
    jy = np.transpose(trim_array_1D(jy_c, nx, ny))
    qx = np.transpose(trim_array_1D(qx_c, nx, ny))
    qy = np.transpose(trim_array_1D(qy_c, nx, ny))
    wt = np.transpose(trim_array_1D(wt, nx, ny))
    prof_Z = np.transpose(trim_array_1D(prof_Z, nx, ny))
    Ue = np.transpose(trim_array_1D(Ue, nx, ny))

    #-------

    dxT, dyT = get_grad_varZ(x_grid, y_grid, Te)
    ne = np.transpose(trim_array_1D(ne, nx, ny))
    Te = np.transpose(trim_array_1D(Te, nx, ny))
    dxT = np.transpose(trim_array_1D(dxT, nx, ny))
    dyT = np.transpose(trim_array_1D(dyT, nx, ny))

    dxB, dyB = get_grad_varZ(x_grid, y_grid, Bz)
    Bz = np.transpose(trim_array_1D(Bz, nx, ny))
    dxB = np.transpose(trim_array_1D(dxB, nx, ny))
    dyB = np.transpose(trim_array_1D(dyB, nx, ny))

    rA = (Z2ni)**-1

    dict = {}

    dict['ne'] = ne
    dict['Te'] = Te
    dict['Cx'] = load_mat_1D(path, fprefix, time, 'Cx', nx, ny)
    dict['U'] = Ue

    dict['qx'] = qx
    dict['qy'] = qy
    dict['jx'] = jx
    dict['jy'] = jy
    dict['dxT'] = dxT
    dict['dyT'] = dyT
    dict['Bz'] = Bz
    dict['dxB'] = dxB
    dict['dyB'] = dyB
    dict['x_grid'] = x_grid
    dict['y_grid'] = y_grid
    dict['v_grid'] = grid['v_grid']
    dict['grid'] = grid
    dict['Z2ni'] = Z2ni
    dict['wt'] = wt
    dict['Z'] = prof_Z
    dict['fo'] = fo
    dict['nv'] = nv
    dict['ny'] = ny
    dict['nx'] = nx
    dict['time'] = dict_te['time']

    return dict


#-----------------------------------------------------------------------
def get_grad_3d(grid, data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata,dydata = get_grad(x_grid,y_grid,data)
    '''
    if len(np.shape(data)) == 1:
        x_grid = grid['x_grid']
        dx = x_grid[2:] - x_grid[:-2]
        dxT = (data[2:] - data[:-2]) / dx
        return dxT
    elif len(np.shape(data)) == 2:
        x_grid = grid['x_grid']
        y_grid = grid['y_grid']
        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        dxT = np.zeros(nv, (len(dx), len(dy)))
        dyT = np.zeros((len(dx), len(dy)))
        ldx, ldy = len(dx), len(dy)
        #print 'SHAPE DX DY = ', ldx,ldy
        for ii in range(len(dx)):
            #print 'ii=  ', ii
            dyT[ii, :] = (data[ii, 2:] - data[ii, :-2]) / dy

        for ii in range(len(dy)):
            #print 'ix = ', ii, np.shape((data[2:,ii] - data[:-2,ii])), np.shape(dxT), np.shape(dx)
            dxT[:, ii] = (data[2:, ii] - data[:-2, ii]) / dx
    else:
        nv, ny, nx = np.shape(data)
        x_grid = grid['x_grid']
        y_grid = grid['y_grid']

        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        x, y, z = np.meshgrid(np.ones(ny), np.ones(nv), np.ones(nx))

        dy_3d, y, z = np.meshgrid(dy, np.ones(nv), np.ones(nx))
        x, y, dx_3d = np.meshgrid(np.ones(ny), np.ones(nv), dx)
        dxT = np.zeros(np.shape(data))
        dyT = np.zeros(np.shape(data))

        dxT[:, :, 1:-1] = (data[:, :, 2:] - data[:, :, :-2]) / dx_3d
        dyT[:, 1:-1, :] = (data[:, 2:, :] - data[:, :-2, :]) / dy_3d
        # Boundary conditions
        dxT[:, :, 0] = dxT[:, :, 1]
        dxT[:, :, -1] = dxT[:, :, -2]
        if ny != 1:
            dxT[:, 0, :] = dxT[:, 1, :]
            dxT[:, -1, :] = dxT[:, -2, :]

        return dxT, dyT


#-----------------------------------------------------------------------
def get_grad_3d_varZ(grid, data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata,dydata = get_grad(x_grid,y_grid,data)
    '''
    if len(np.shape(data)) == 1:
        x_grid = grid['x_grid']
        dx = x_grid[2:] - x_grid[:-2]
        dxT = (data[2:] - data[:-2]) / dx
        return dxT
    elif len(np.shape(data)) == 2:
        x_grid = grid['x_grid']
        y_grid = grid['y_grid']
        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        dxT = np.zeros(nv, (len(dx), len(dy)))
        dyT = np.zeros((len(dx), len(dy)))
        ldx, ldy = len(dx), len(dy)
        #print 'SHAPE DX DY = ', ldx,ldy
        for ii in range(len(dx)):
            #print 'ii=  ', ii
            dyT[ii, :] = (data[ii, 2:] - data[ii, :-2]) / dy

        for ii in range(len(dy)):
            #print 'ix = ', ii, np.shape((data[2:,ii] - data[:-2,ii])), np.shape(dxT), np.shape(dx)
            dxT[:, ii] = (data[2:, ii] - data[:-2, ii]) / dx
    else:
        nv, ny, nx = np.shape(data)
        x_grid = grid['x_grid']
        y_grid = grid['y_grid']
        dx = x_grid[2:] - x_grid[:-2]
        dy = y_grid[2:] - y_grid[:-2]
        x, y, z = np.meshgrid(np.ones(ny), np.ones(nv), np.ones(nx))
        dy_3d, y, z = np.meshgrid(dy, np.ones(nv), np.ones(nx))
        x, y, dx_3d = np.meshgrid(np.ones(ny), np.ones(nv), dx)
        dxT = np.zeros(np.shape(data))
        dyT = np.zeros(np.shape(data))

        if ny > 3:
            dyT[:, 1:-1, :] = (data[:, 2:, :] - data[:, :-2, :]) / dy_3d

            dyT[:, 0, :] = dyT[:, 1, :]
            dyT[:, -1, :] = dyT[:, -2, :]

        if nx > 3:
            dxT[:, :, 1:-1] = (data[:, :, 2:] - data[:, :, :-2]) / dx_3d
            # Boundary conditions
            dxT[:, :, 0] = dxT[:, :, 1]
            dxT[:, :, -1] = dxT[:, :, -2]
        return dxT, dyT


#-----------------------------------------------------------------------


def calc_q_T(qx_c, qy_c, jx_c, jy_c, ne, Cx, Cy, Ue):
    '''
        
    '''
    c2 = Cx**2 + Cy**2
    #qx_T = qx_c + 0.5*ne*((c2 + 2.0*Ue)*Cx - c2*jx_c - 2.0*Cx*(Cx*jx_c + Cy*jy_c))
    #qx_T = qx_c + 0.5*ne*((c2 + 2.0*Ue)*Cx - c2*jx_c - 2.0*Cx*(Cx*jx_c + Cy*jy_c))
    qx_T = qx_c + 0.5 * ne * (c2 * Cx) + Ue * Cx - (c2 * jx_c + 2.0 * Cx *
                                                    (Cx * jx_c + Cy * jy_c)) * 0.5

    qy_T = qy_c + 0.5 * ne * (c2 * Cy) + Ue * Cy - 0.5 * (c2 * jy_c + 2.0 * Cy *
                                                          (Cx * jx_c + Cy * jy_c))
    return qx_T, qy_T


#---------------------------------------------------


def get_q_T(path, fprefix, time):
    '''
        
    '''
    Cx_dict = load_dict(path, fprefix, 'Cx', time)
    Cy_dict = load_dict(path, fprefix, 'Cy', time)
    n_dict = load_dict(path, fprefix, 'n', time)
    Ue_dict = load_dict(path, fprefix, 'U', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    qx_dict = load_dict(path, fprefix, 'qxX', time)
    qy_dict = load_dict(path, fprefix, 'qyY', time)
    jx_dict = load_dict(path, fprefix, 'jxX', time)
    jy_dict = load_dict(path, fprefix, 'jyY', time)

    Te = Te_dict['mat']
    Cx = Cx_dict['mat']
    Cy = Cy_dict['mat']
    ne = n_dict['mat']
    Ue = Ue_dict['mat']
    jxX = jx_dict['mat']
    jyY = jy_dict['mat']
    qxX = qx_dict['mat']
    qyY = qy_dict['mat']

    #Cx,Cy,ne,Ue,qxX,qyY,jxX,jyY = Cx_dict['mat'],Cy_dict['mat'],n_dict['mat'], Ue_dict['mat'], qx_dict['mat'], qy_dict['mat'], jx_dict['mat'], jy_dict['mat']

    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])
    jx_c = 0.5 * (jxX[1:, :] + jxX[:-1, :])
    jy_c = 0.5 * (jyY[:, 1:] + jyY[:, :-1])

    nx, ny = np.shape(jx_c)
    #print 'ARRAY SHAPES === ', np.shape(Cx), np.shape(Cy),np.shape(ne),np.shape(Ue), np.shape(qx_c), np.shape(qy_c), np.shape(jx_c),np.shape(jy_c)

    Cx = trim_array(Cx, nx, ny)
    Cy = trim_array(Cy, nx, ny)
    ne = trim_array(ne, nx, ny)
    Ue = trim_array(Ue, nx, ny)
    Te = trim_array(Te, nx, ny)
    Diff = (1.5 * Te * ne) - Ue
    #Cx = trim_array(Cx,nx,ny)
    #print ' DIFF ARRAY = ', Diff,'\n max dif = ', np.max(Diff)

    #print 'ARRAY SHAPES === ', np.shape(Cx), np.shape(Cy),np.shape(ne),np.shape(Ue), np.shape(qx_c), np.shape(qy_c), np.shape(jx_c),np.shape(jy_c)
    qxT, qyT = calc_q_T(qx_c, qy_c, jx_c, jy_c, ne, Cx, Cy, Ue)
    return qxT, qyT


#-----------------------------------------------------------------------


def plot_ax(ax1, x_grid, data, norm_const, c='b', tlab='00', ls='-'):
    if len(np.shape(data)) > 1:
        if len(data[0, :]) < 3:
            x_gr = x_grid[:len(data[:, 0])]
            ax1.plot(x_gr,
                     data[:, 1] * norm_const,
                     c=c,
                     linestyle=ls,
                     label=r'\textit{$' + str(tlab) + '$}')
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        else:
            x_gr = x_grid[:len(data[:, 0])]
            ax1.plot(x_gr,
                     data[:, yi] * norm_const,
                     c=c,
                     linestyle=ls,
                     label=r'\textit{$' + str(tlab) + '$}')
            ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
    else:
        x_gr = x_grid[:len(data)]
        ax1.plot(x_gr,
                 data[:, yi] * norm_const,
                 c=c,
                 linestyle=ls,
                 label=r'\textit{$' + str(tlab) + '$}')

    return


#-----------------------------------------------------------------------
def plot_2D_general(ax,
                    data,
                    label,
                    middle=None,
                    colormap='jet',
                    limits=None,
                    xlabel=None,
                    ylabel=None,
                    clabel=None,
                    vmin=None,
                    vmax=None):
    if middle != None:
        norm = MidPointNorm(middle)
    else:
        norm = None

    if limits:
        im = ax.imshow(data,
                       aspect='auto',
                       cmap=colormap,
                       norm=norm,
                       extent=limits,
                       vmin=vmin,
                       vmax=vmax)
        plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)

    else:
        im = ax.imshow(data, aspect='auto', cmap=colormap, norm=norm, vmin=vmin, vmax=vmax)
        plt.colorbar(im, ax=ax, aspect='auto', norm=norm, label=clabel)
    ax.set_title(label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    #print ' plotted ---- ', label


#-----------------------------------------------------------------------


def plot_xhlines(ax,
                 x_grid_im,
                 lim=[],
                 lineout_list=[20, 40, 50, 73],
                 cmap=cm.bone,
                 linestyle='--',
                 dashes=(4, 2)):
    color_lineout = cmap(np.linspace(0.0, 1, len(lineout_list) + 1))
    color_lineout = color_lineout[:-1]

    if len(lineout_list) == 0:
        return
    for xl in range(len(lineout_list)):

        ax.axhline(x_grid_im[lineout_list[xl]],
                   c=color_lineout[xl],
                   linestyle=linestyle,
                   dashes=dashes)


#-----------------------------------------------------------------------
def plot_xvlines(ax,
                 x_grid_im,
                 lim,
                 lineout_list=[20, 40, 50, 73],
                 color_lineout=[],
                 cmap=cm.bone,
                 custom_cmap=[],
                 linestyle='--',
                 dashes=(4, 2)):
    if len(color_lineout) == 0:
        if len(custom_cmap) == 0:
            color_lineout = cmap(np.linspace(0.0, 1, len(lineout_list) + 1))
            color_lineout = color_lineout[:-1]
        else:
            color_lineout = custom_cmap

    if len(lineout_list) == 0:
        return
    for xl in range(len(lineout_list)):

        ax.axvline(x_grid_im[lineout_list[xl]],
                   c=color_lineout[xl],
                   linewidth=1.2,
                   linestyle=linestyle,
                   dashes=dashes)


#-----------------------------------------------------------------------


def load_labeldict():
    '''
       dict = cf.load_labeldict()
    '''
    dict = {}
    var_list = ['q_SH_x', 'q_SH_y', 'q_RL_x', 'q_RL_y', 'q_TE_x', 'q_TE_y', 'q_E_x', 'q_E_y']

    for var in var_list:
        dict[var] = {}
        qS = re.search('\w_(?P<type>\w+)_(?P<component>\w+)', var)
        type = qS.group('type')
        component = qS.group('component')
        dict[var]['title'] = r'$q_{' + type + ',' + component + '}$'
        dict[var]['unit'] = r'[$n_0 m_e v_n^3$]'
    return dict


#-----------------------------------------------------------------------


def produce_2D_fig(x_grid, data, lab_dict, lineout_list=[20, 40, 50, 73]):

    figsize = lab_dict['figsize']
    lims = lab_dict['lims']
    title = lab_dict['title']
    cbar_title = lab_dict['cbar_title']
    middle = lab_dict['middle']
    colormap = lab_dict['colormap']
    xlab, ylab = lab_dict['xlab'], lab_dict['ylab']
    save_name = lab_dict['save_name']
    if lab_dict.has_key('vmin') and lab_dict.has_key('vmax'):
        vmin = lab_dict['vmin']
        vmax = lab_dict['vmax']
    else:
        vmin = np.min(data)
        vmax = np.max(data)
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(1, 1, 1)
    plot_2D_general(ax1,
                    data,
                    title,
                    middle,
                    colormap,
                    lims,
                    ylab,
                    xlab,
                    cbar_title,
                    vmin=vmin,
                    vmax=vmax)
    lim_x = ax1.get_xlim()
    plot_xhlines(ax1, x_grid, lim_x, lineout_list)

    fig.savefig(save_name, dpi=400)
    plt.close(fig)


#-----------------------------------------------------------------------


def plot_twodim(ax1, path, fprefix, var, time='00', cmap='RdBu_r'):
    norm_const, c_title, c_fmt = calc_norms(var)
    fname = construct_fname(path, fprefix, var, time)
    time_f = get_time(fname)
    out_l = get_startline(fname)
    leg_list.append(get_time(fname))
    data = np.loadtxt(fname, skiprows=out_l)
    data = np.where(abs(data) <= 1e-99, 0.0, data)

    neg_list = {}
    min_data = np.min(data)
    max_data = np.max(data)
    norm_const, c_title, c_fmt = calc_norms(var, min_data)

    if var in neg_list:
        norm = MidPointNorm(midpoint=0.0)
        colormap = cmap
    elif min_data == max_data:
        norm = None
    else:
        middle = 0.5 * (min_data + max_data) * norm_const
        norm = MidPointNorm(midpoint=middle)

        colormap = cmap

    if (var != 'Cx') and (var != 'Cy'):
        dict1 = fpg_get_info(fname)

    if len(dict1['x_grid']) != np.shape(data)[0]:
        x_c_temp = dict1['x_grid'][1:-1] * xstep_factor
        #print 'Mod SHAPE data: ', np.shape(data), np.shape(x_c_temp)
    else:
        x_c_temp = dict1['x_grid'] * xstep_factor
        #print ' SHAPE data: ', np.shape(data), np.shape(x_c_temp)

    x_c_grid = x_c_temp[::-1]
    X, Y = np.meshgrid(dict1['y_grid'] * xstep_factor, x_c_grid)
    lims = [
        dict1['y_grid'][0] * xstep_factor, dict1['y_grid'][-1] * xstep_factor,
        dict1['x_grid'][0] * xstep_factor, dict1['x_grid'][-1] * xstep_factor
    ]

    im = ax1.imshow(data * norm_const,
                    cmap=colormap,
                    aspect=asp,
                    norm=norm,
                    vmin=min_data * norm_const,
                    vmax=max_data * norm_const,
                    extent=lims)

    title = r'$ ' + var + '$'
    ax1.set_xlabel(ylab)
    ax1.set_ylabel(xlab)

    if c_fmt[-1] == 'e':
        fmt_cbar = ticker.FuncFormatter(fmt)
    else:
        fmt_cbar = ticker.FuncFormatter(fmt_ord)

    plt.colorbar(im, format=fmt_cbar, shrink=0.9, label=c_title)


#-----------------------------------------------------------------------
def plot_twodim_dict(ax1, dict, var, colormap='RdBu_r'):
    norm_const, c_title, c_fmt = calc_norms(var)
    leg_list = []

    time_f = dict['time']
    leg_list.append(time_f)

    data = dict['mat']
    x_grid = dict['x_grid']
    y_grid = dict['y_grid']
    ny, nx = len(y_grid), len(x_grid)
    neg_list = {}
    min_data = np.min(data)
    max_data = np.max(data)
    norm_const, c_title, c_fmt = calc_norms(var, min_data)

    if var in neg_list:
        norm = MidPointNorm(midpoint=0.0)
        #colormap=cmap#plt.cm.seismic
    elif min_data == max_data:
        norm = None
    else:
        middle = 0.5 * (min_data + max_data) * norm_const
        norm = MidPointNorm(midpoint=middle)

        #colormap=cmap

    if (var != 'Cx') and (var != 'Cy'):
        dict = fpg_get_info(fname)

    #print ' np.len(x_grid) = ', np.shape(x_grid)
    #print ' np.shape(data) = ', np.shape(data)
    if len(x_grid) != np.shape(data)[0]:
        if len(x_grid) > np.shape(data)[0] and np.shape(data)[-1] >= ny:
            x_c_temp = x_grid[:len(data[:, 0])] * xstep_factor
            y_grid = y_grid[:len(data[0, :])]
        else:
            data = trim_array(data, nx, ny)
        #print 'Mod SHAPE data: ', np.shape(data), np.shape(x_c_temp)
    else:
        x_c_temp = x_grid * xstep_factor

    y_grid = y_grid * xstep_factor
    lims = [y_grid[0], y_grid[-1], x_grid[0], x_grid[-1]]

    x_c_grid = x_c_temp[::-1]
    X, Y = np.meshgrid(y_grid, x_grid)
    im = ax1.imshow(data * norm_const,
                    cmap=colormap,
                    aspect=asp,
                    norm=norm,
                    vmin=min_data * norm_const,
                    vmax=max_data * norm_const,
                    extent=lims)

    title = r'$ ' + var + '$'
    ax1.set_xlabel(ylab)
    ax1.set_ylabel(xlab)

    if c_fmt[-1] == 'e':
        fmt_cbar = ticker.FuncFormatter(fmt)
    else:
        fmt_cbar = ticker.FuncFormatter(fmt_ord)

    plt.colorbar(im, format=fmt_cbar, shrink=0.9, label=c_title)


#-----------------------------------------------------------------------


def get_path_style(final_lab, h_lab, B_lab):
    '''
        lstyle,marker = get_path_style(h_lab,B_lab)
    '''
    lstyle = '-'
    marker = None
    if B_lab == 'no B':
        marker = 'o'

    if B_lab[-1] == 'H':
        if B_lab[0] == '0':
            marker = None
            lstyle = '-'

        elif B_lab[0] == '1':
            print len(B_lab)
            marker = 'x'
            if len(B_lab) > 2:
                marker = 'o'
    else:
        marker = None
    if h_lab == 'static':
        lstyle = ':'
    elif h_lab == 'no$ $C_y':
        lstyle = ':'
        marker = 'x'
    elif B_lab[-2:] == 'MW':
        lstyle = '--'    # Maxwellian run
    return lstyle, marker


#-----------------------------------------------------------------------


def get_path_style_path(fprefix):
    '''
        lstyle,marker = get_path_style(h_lab,B_lab)
    '''
    lstyle = '-'
    marker = None
    if B_lab == 'no B':
        marker = 'o'
    #if B_lab[-1] == 'T':
    #    marker = 'x'

    else:
        marker = None
    if h_lab == 'static':
        lstyle = ':'
    elif h_lab == 'no$ $C_y':
        lstyle = ':'
        marker = 'x'
    elif B_lab[-2:] == 'MW':
        lstyle = '--'    # Maxwellian run
    return lstyle, marker


#-----------------------------------------------------------------------


def get_q_over_qFS(path, fprefix, time):
    '''
        qx_over_qFS,qy_over_qFS = get_q_over_qFS(path,fprefix,time)
        
        '''
    ne_dict = load_dict(path, fprefix, 'n', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    qx_dict = load_dict(path, fprefix, 'qxX', time)
    qy_dict = load_dict(path, fprefix, 'qyY', time)
    ne = ne_dict['mat']
    Te = Te_dict['mat']
    qxX = qx_dict['mat']
    qyY = qy_dict['mat']
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])

    nx, ny = np.shape(qx_c)
    ne = trim_array(ne, nx, ny)
    Te = trim_array(Te, nx, ny)
    qx_over_qFS = 2.0 * qx_c * ((ne * (2.0 * Te)**1.5)**-1)
    qy_over_qFS = 2.0 * qy_c * ((ne * (2.0 * Te)**1.5)**-1)
    return qx_over_qFS, qy_over_qFS


#-----------------------------------------------------------------------


def get_Ptotal(path, fprefix, time, Z=6.51):
    '''
   
    P_total = ne Te + ni Ti
            = 2/3 Ue + (ne/Z)(2/3)*Ui
            Note Ui = (3/2)Ti 
    ptotal = get_Ptotal(path,fprefix,time,Z=6.51)
    '''

    ne_dict = load_dict(path, fprefix, 'n', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    ue_dict = load_dict(path, fprefix, 'U', time)
    try:
        print ' path_prefix = ', path, fprefix, 'calc, Pi'
        Ui_dict = load_dict(path, fprefix, 'Ui', time)
        Ui = Ui_dict['mat']
        ue = ue_dict['mat']
        ne = ne_dict['mat']
        nx, ny = np.shape(ne)

        Ui = trim_array(Ui, nx, ny)
        ue = trim_array(ue, nx, ny)

        Pi = (2.0 / 3.0) * (Ui * ne) * (1.0 / Z)
    except ValueError:
        print ' path_prefix = ', path, fprefix, '==========skipping Pi'

        ue = ue_dict['mat']
        Pi = np.zeros((np.shape(ue)))
    ptotal = ((2.0 / 3.0) * ue) + Pi
    return ptotal


def get_qFS(path, fprefix, time):
    '''
        qFS =  get_qFS(path,fprefix,time)
        
        '''

    ne_dict = load_dict(path, fprefix, 'n', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    Te_dict = load_dict(path, fprefix, 'Te', time)
    qx_dict = load_dict(path, fprefix, 'qxX', time)
    qy_dict = load_dict(path, fprefix, 'qyY', time)

    #print ' the time of ', path, ' is: ', ne_dict['time']
    ne = ne_dict['mat']
    Te = Te_dict['mat']
    qxX = qx_dict['mat']
    qyY = qy_dict['mat']
    qx_c = 0.5 * (qxX[1:, :] + qxX[:-1, :])
    qy_c = 0.5 * (qyY[:, 1:] + qyY[:, :-1])

    nx, ny = np.shape(qx_c)
    ne = trim_array(ne, nx, ny)
    Te = trim_array(Te, nx, ny)
    qFS = 0.5 * ne * (2.0 * Te)**1.5
    return qFS


def get_speckle_fname_detail(name):
    '''
        tcoh,rad,B = get_speckle_fname_detail(name)
    '''
    s = re.search('tc(?P<tcoh>\d+.\d+)ps_d(?P<rad>\d+.\d+)um', name)
    if s:
        tcoh = s.group('tcoh')
        rad = s.group('rad')
        print 'tcoh = ', tcoh, 'rad = ', rad
    noB = re.search('noB', name)
    B = True
    if noB:
        B = False
        print ' no B'
    return tcoh, rad, B


#-----------------------------------------------------------------------
def get_avgx(mat, ax=1):
    '''
        ax  = [n,..,2,1,0]
        ax = 1 x axis
        ax = 0 y axis
    '''

    avg = np.average(mat, axis=ax)
    return avg


#-----------------------------------------------------------------------
def get_U_dev(U_data, lim=0.01):
    '''
        Gets a measure of the fractional deviation from average
       saturates at lim
    '''
    avgx_U = get_avgx(U_data)
    x_dummy, avgxy_U = np.meshgrid(np.ones((len(U_data[0, :]))), avgx_U)
    U_dev = U_data / avgxy_U
    U_dev = np.abs(U_data - avgxy_U)
    if lim > 0.0:
        U_dev = np.where(U_dev <= lim, U_dev, lim)
    return U_dev


#-----------------------------------------------------------------------
def get_U_dev_abs(U_data):
    '''
        Gets a measure of the fractional deviation from average
       saturates at lim - the perturbation amplitude
       U_dev = get_U_dev_abs(U_data)
    '''
    avgx_U = get_avgx(U_data)
    x_dummy, avgxy_U = np.meshgrid(np.ones((len(U_data[0, :]))), avgx_U)
    #U_dev = U_data/avgxy_U
    U_dev = (U_data - avgxy_U)
    #if lim>0.0:
    #    U_dev = np.where(U_dev<=lim,U_dev,lim)
    return U_dev


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def get_U_dev_frac(U_data):
    '''
        Gets a measure of the fractional deviation from average
       saturates at lim - the perturbation amplitude
       U_dev = get_U_dev_abs(U_data)
    '''
    avgx_U = get_avgx(U_data)
    x_dummy, avgxy_U = np.meshgrid(np.ones((len(U_data[0, :]))), avgx_U)
    #U_dev = U_data/avgxy_U
    U_dev = (U_data - avgxy_U) / avgxy_U
    #if lim>0.0:
    #    U_dev = np.where(U_dev<=lim,U_dev,lim)
    return U_dev


#-----------------------------------------------------------------------
def get_sigma_rms(U_data, reduced=True):
    '''
        Gets a measure of the fractional deviation from average
       saturates at lim
    '''
    if reduced:
        U_data = U_data[:, 1:-1]

    avgx_U = get_avgx(U_data)
    x_dummy, avgxy_U = np.meshgrid(np.ones((len(U_data[0, :]))), avgx_U)
    summy = np.sqrt(np.sum((U_data - avgxy_U)**2, axis=1))
    U_dev = summy / avgxy_U[:, 0]

    return summy


#-----------------------------------------------------------------------


def annotate_time(ax, lett='(a)', dx_mult=1.0, dy_mult=1.0, loc='top', fontsize=0):
    '''
        puts time in the top right corner
    '''
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    nunits = 20.0
    dx = (xmax - xmin) / nunits
    dy = (ymax - ymin) / nunits
    if loc == 'top':
        print 'x = ', xmin - dx, 'y = ', ymax
        x_coord = xmax - 10.0 * dx * dx_mult
        y_coord = ymax - 5.0 * dy * dy_mult
    else:
        print 'x = ', xmin - dx, 'y = ', ymax
        x_coord = xmax - 10.0 * dx * dx_mult
        y_coord = ymin + (dy * dy_mult)
    if fontsize != 0:
        ax.text(x_coord, y_coord, lett, fontsize=fontsize)

    else:
        ax.text(x_coord, y_coord, lett)
    return


def get_atomic_Z(fname):
    '''
        Searches fname - which should be either fort.10 or an IMPACT rundeck to retrieve the atomic number
        19/02/19 - generated
    '''
    f = open(fname, 'r')
    data = f.readlines()
    Z_fl = 0.0
    for line in data:
        a = re.search('atomic_Z\s+=\s+(?P<atomicZ>[0-9.d]+)', line)
        if a:
            Z_fl = float(re.sub('d', 'e', a.group('atomicZ')))
    f.close()
    return Z_fl


#------------------------------------------------------------------------------
def get_dims(fname):
    '''
        Searches fname - which should be either fort.10 or an IMPACT rundeck to retrieve the run dimensions
        19/02/19 - generated
        
        nv,ny,nx = get_dims(fname)
    '''
    f = open(fname, 'r')
    data = f.readlines()
    nv, ny, nx = 0, 0, 0
    for line in data:
        a_nx = re.search('nx\s+=\s+(?P<nx>[0-9]+)', line)
        a_ny = re.search('ny\s+=\s+(?P<ny>[0-9]+)', line)
        a_nv = re.search('nv\s+=\s+(?P<nv>[0-9]+)', line)

        if a_nx:
            nx = int(a_nx.group('nx'))
        if a_ny:
            ny = int(a_ny.group('ny'))
        if a_nv:
            nv = int(a_nv.group('nv'))
    f.close()
    return nv, ny, nx


#------=================================================================
#-----------------------------------------------------------------------
def set_ylim_max(ax_in, x_grid, data, xlim, y_mult=[1.0, 1.0]):
    '''
        This funciton sets the max and min y values of a particular plot
        to sensible values. for the given data
        (ymin-0.05(ymax-ymin) and ymax + 0.05(ymax-ymin)
        
    '''
    xmin, xmax = xlim
    ymin = np.min(data[(x_grid < xmax) * (x_grid >= xmin)])
    ymax = np.max(data[(x_grid < xmax) * (x_grid >= xmin)])
    dy = 0.05 * np.abs(ymax - ymin)
    if ymin != ymax:

        ax_in.set_ylim(ymin - dy * y_mult[0], ymax + dy * y_mult[1])
    else:
        print('error ymin = ymax data all the same!')
    print('ymin = %4.4e ymax = %4.4e ' % (ymin, ymax))

    return


#------------------------------------------------------------------------------
def data_intersect(x_grid,
                   data,
                   data_pos=np.array([100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0])):
    '''
        i_out = data_intersect(data,data_pos)
        Find the points where data is closest to values in data_pos
    '''
    i_out = np.zeros((len(data_pos)), dtype=int)
    for itt, data_loc in enumerate(data_pos):
        i_closest = np.argsort(np.abs(data[:-100] - data_loc))[0]
        i_out[itt] = i_closest
    return i_out


#-----------------------------------------------------------------------
def get_wt(path, time):
    '''
        CURRENT IMPACT RUNS HAVE INCORRECT HALL PARAMETER OUTPUTS
        Missing a 1/Z2ni factor, where 1/Z2ni = 1/<z^2>ni = 1/(fixed_prof_Z*ne)
        
        
        Note fixed_prof_Z = <Z^2>/<Z> in these new IMPACT runs with TF EOS.
        Use of this function should still work for other IMPACT versions (With older fixed_prof_Z but should be checked).
        
        
    '''
    fprefix = fpre(path)
    dict_wt = load_dict(path, fprefix, 'wt', time)

    dict_fo = load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']
    nv, ny, nx = np.shape(fo)

    dict_ne = load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    dict_Z = load_dict(path, fprefix, 'Z', time)
    prof_Z = dict_Z['mat']

    Z2ni = trim_array(ne, nx, ny) * trim_array(prof_Z, nx, ny)
    wte_uncorrected = trim_array(dict_wt['mat'], nx, ny)

    return wte_uncorrected / Z2ni


#-----------------------------------------------------------------------
def get_wt_1D(path, time):
    '''
        CURRENT IMPACT RUNS HAVE INCORRECT HALL PARAMETER OUTPUTS
        Missing a 1/Z2ni factor, where 1/Z2ni = 1/<z^2>ni = 1/(fixed_prof_Z*ne)
        
        
        Note fixed_prof_Z = <Z^2>/<Z> in these new IMPACT runs with TF EOS.
        Use of this function should still work for other IMPACT versions (With older fixed_prof_Z but should be checked).
        
        
    '''
    fprefix = fpre(path)
    dict_wt = load_dict_1D(path, fprefix, 'wt', time)
    dict_ne = load_dict_1D(path, fprefix, 'n', time)
    ne = dict_ne['mat']
    dict_Z = load_dict_1D(path, fprefix, 'Z', time)
    prof_Z = dict_Z['mat']

    Z2ni = ne * prof_Z

    return dict_wt['mat'] / Z2ni
