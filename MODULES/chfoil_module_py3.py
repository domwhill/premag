import numpy as np, re, os, sys
from pylab import *
import impact_norms_py3 as INN


class conv_factors_custom(object):
    # --- arguments given to the class are passed straight to init...
    #
    #def __init__(self,norm_path,Z=3.5,Ar=6.51):
    def __init__(self, norm_path, Z=5.954127788543701172, Ar=6.51):
        self.norm_name = self.search_path_class('norm.txt', norm_path)
        if norm_path[-1] != '/':
            norm_path = norm_path + '/'
        if os.path.exists(norm_path + 'indices.txt'):
            icl, icu = np.loadtxt(norm_path + 'indices.txt')
            self.cl_index = int(icl)    #80#
            self.c_index = int(icu)    #207#
        else:
            self.cl_index = 0
            self.c_index = 100
        self.SI_on = True
        self.Z = Z
        self.Ar = Ar
        [self.T0, self.n0, self.logLambda] = np.loadtxt(norm_path + self.norm_name)
        dict_norm = INN.impact_inputs(self.n0, self.T0, Z, 1.0, Ar)
        self.v_te = dict_norm['vte']    # m/s
        self.tau_ei = dict_norm['tau_ei']    # s
        self.nu_ei = 1.0 / self.tau_ei    # s^-1
        self.lambda_mfp = dict_norm['lambda_mfp']    # m
        self.Bz_ref = dict_norm['Bz_norm']
        self.yi = 10
        self.lineout_list = [20, 40, 50, 73]
        self.lineout_list = np.linspace(self.cl_index, self.c_index, 3, dtype=int)

        interval = self.c_index - self.cl_index
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
        return

    #-------------------------------------------------------------------
    def search_path_class(self, regexp, path):
        '''
            out_name = search_path(regexp,path)
        '''
        out_name = ''
        for pp in os.listdir(path):
            rr = re.search(regexp, pp)
            if rr:
                out_name = pp
                break
        return out_name


#-------------------------------------------------------------------====


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


def list_to_float(input):
    list = input.split()
    arr = np.zeros((len(list)), dtype=float)
    arr[:] = list[:]
    return arr


def construct_fname(path, fprefix, var, time):

    if var == 'fo' or var == 'fxX' or var == 'fxY' or var == 'fyX' or var == 'fyY':
        suffix = '.xyv'
    else:
        suffix = '.xy'

    #fprefix = 'thydro_hi'
    #time = '00'
    fname = path + '/' + fprefix + '_' + var + '_' + time + suffix
    return fname


def load_dict(path, fprefix, var, time):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''
    fname = construct_fname(path, fprefix, var, time)

    #dict = MyDict()
    #mat = np.loadtxt(fname,skiprows=out_l)

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
            offset_x = int((len(array[:]) % nx) / 2)
            out_array = array[offset_x:-offset_x]
            return out_array

    if np.shape(array) == (nx, ny):
        return array
    else:
        out_array = np.zeros((nx, ny))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -nx, -ny

        if nx_a != nx:
            offset_x = int((len(array[:, 0]) % nx) / 2)
            offset_x_rh = int((len(array[:, 0]) % nx) / 2)
            offset_x += int(2.0 * np.abs((len(array[:, 0]) % nx) / 2.0 -
                                         (len(array[:, 0]) % nx) / 2))

        if ny_a != ny:
            offset_y = int((len(array[0, :]) % ny) / 2)
            offset_y_rh = int((len(array[0, :]) % ny) / 2)
            offset_y += int(2.0 * np.abs((len(array[0, :]) % ny) / 2.0 -
                                         (len(array[0, :]) % ny) / 2))
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
            print(' ERROR ARRAY NOT TRIMMED PROPERLY ')
            print(' SHAPE ARRAY IN = ', np.shape(array))
            print(' SHAPE ARRAY OUT = ', np.shape(out_array),
                  ' shape asked for = (%i,%i)' % (nx, ny))
            sys.exit()
            return array
        return out_array


def get_grad(x_grid, y_grid, T_data):
    '''
        ONLY FOR CC cells - centred differencing
        dxdata,dydata = get_grad(x_grid,y_grid,data)
    '''

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
        for ii in range(len(dx)):
            dyT[ii, :] = (T_data[ii, 2:] - T_data[ii, :-2]) / dy

        for ii in range(len(dy)):
            dxT[:, ii] = (T_data[2:, ii] - T_data[:-2, ii]) / dx
        return dxT, dyT


def trim_array_1D(array, nx, ny):
    '''
        takes the central (nx,ny) slab of any array
        out_array = trim_array(array,nx,ny)
    '''
    if ny == 0 and len(np.shape(array)) == 1:
        out_array = np.zeros((nx))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -int(nx), -int(ny)
        nx_a = len(array)

        if nx_a != nx:
            offset_x = int(((len(array) % nx) / 2))
            offset_x_rh = int((len(array) % nx) / 2)
            offset_x += int(2.0 * np.abs((len(array) % nx) / 2.0 - (len(array) % nx) / 2))

        if offset_x_rh != 0:
            out_array = array[offset_x:-offset_x_rh]
        elif offset_x_rh == 0:
            out_array = array[offset_x:]
        else:
            out_array = array[offset_x:]
        if len(out_array) != nx:
            print(' ERROR ARRAY NOT TRIMMED PROPERLY - trim 1D')
            print(' SHAPE ARRAY IN = ', np.shape(array))
            print(' SHAPE ARRAY OUT = ', np.shape(out_array),
                  ' shape asked for = (%i,%i)' % (nx, ny))
            sys.exit()
            return array
        return out_array
    # 2D stuff --->
    #======================================================?
    if ny == 3 and np.shape(array)[1] == 1:
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
            offset_x = int((len(array[:]) % nx) / 2)
            out_array = array[offset_x:-offset_x]
            return out_array

    if np.shape(array) == (nx, ny):
        return array
    else:
        out_array = np.zeros((nx, ny))
        offset_x, offset_y = 0, 0
        offset_x_rh, offset_y_rh = -nx, -ny

        if nx_a != nx:
            offset_x = int((len(array[:, 0]) % nx) / 2)
            offset_x_rh = int((len(array[:, 0]) % nx) / 2)
            offset_x += int(2.0 * np.abs((len(array[:, 0]) % nx) / 2.0 -
                                         (len(array[:, 0]) % nx) / 2))
        if ny_a != ny:
            offset_y = int((len(array[0, :]) % ny) / 2)
            offset_y_rh = int((len(array[0, :]) % ny) / 2)
            offset_y += int(2.0 * np.abs((len(array[0, :]) % ny) / 2.0 -
                                         (len(array[0, :]) % ny) / 2))

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
            print(' ERROR ARRAY NOT TRIMMED PROPERLY - trim 1D')
            print(' SHAPE ARRAY IN = ', np.shape(array))
            print(' SHAPE ARRAY OUT = ', np.shape(out_array),
                  ' shape asked for = (%i,%i)' % (nx, ny))
            sys.exit()
            return array
        return out_array


def load_data_all(path, fprefix, time):
    '''
        dict  = load_data_all_1D(path,fprefix,time)
    '''
    print(' === dict load all ==== ')
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
    dict_wt = load_dict(path, fprefix, 'wt', time)
    wt = dict_wt['mat']
    dict_ne = load_dict(path, fprefix, 'n', time)
    ne = dict_ne['mat']

    ##Z2ni = ne

    dict_te = load_dict(path, fprefix, 'Te', time)
    Te = dict_te['mat']
    x_grid = dict_te['x_grid']
    y_grid = dict_te['y_grid']

    dict_fo = load_dict(path, fprefix, 'fo', time)
    fo = dict_fo['mat']

    grid = dict_fo
    nv, ny, nx = np.shape(fo)
    print('--- loading grid ---', nv, ny, nx)
    grid['x_grid'] = trim_array_1D(x_grid, nx, 0)
    grid['y_grid'] = trim_array_1D(y_grid, ny, 0)
    print(' y_Grid = ', y_grid)
    #grid['ny'] = 3
    #nv,ny,nx = np.shape(fo)
    dict_Z = load_dict(path, fprefix, 'Z', time)
    prof_Z = dict_Z['mat']
    Z2ni = np.transpose(trim_array(ne, nx, ny)) * np.transpose(trim_array(prof_Z, nx, ny))
    dict_U = load_dict(path, fprefix, 'U', time)
    Ue = dict_U['mat']

    jx = np.transpose(trim_array(jx_c, nx, ny))
    jy = np.transpose(trim_array(jy_c, nx, ny))
    qx = np.transpose(trim_array(qx_c, nx, ny))
    qy = np.transpose(trim_array(qy_c, nx, ny))
    wt = np.transpose(trim_array(wt, nx, ny))
    prof_Z = np.transpose(trim_array(prof_Z, nx, ny))
    Ue = np.transpose(trim_array(Ue, nx, ny))

    #-------

    dxT, dyT = get_grad(x_grid, y_grid, Te)
    dxn, dyn = get_grad(x_grid, y_grid, ne)
    dxBz, dyBz = get_grad(x_grid, y_grid, Bz)

    ne = np.transpose(trim_array(ne, nx, ny))
    Te = np.transpose(trim_array(Te, nx, ny))
    Bz = np.transpose(trim_array(Bz, nx, ny))

    dxT = np.transpose(trim_array(dxT, nx, ny))
    dyT = np.transpose(trim_array(dyT, nx, ny))

    dxn = np.transpose(trim_array(dxn, nx, ny))
    dyn = np.transpose(trim_array(dyn, nx, ny))

    dxBz = np.transpose(trim_array(dxBz, nx, ny))
    dyBz = np.transpose(trim_array(dyBz, nx, ny))

    rA = (Z2ni)**-1

    dict = {}

    dict['ne'] = ne
    dict['Te'] = Te
    dict['U'] = Ue

    dict['qx'] = qx
    dict['qy'] = qy
    dict['jx'] = jx
    dict['jy'] = jy
    dict['dxT'] = dxT
    dict['dyT'] = dyT
    dict['dxn'] = dxn
    dict['dyn'] = dyn

    dict['dxBz'] = dxBz
    dict['dyBz'] = dyBz
    dict['Bz'] = Bz
    dict['x_grid'] = x_grid
    dict['y_grid'] = y_grid
    dict['grid'] = grid
    dict['Z2ni'] = Z2ni
    dict['wt'] = wt * (Z2ni**-1)
    dict['Z'] = prof_Z
    dict['fo'] = fo
    dict['nv'] = nv
    dict['ny'] = ny
    dict['nx'] = nx
    dict['time'] = dict_te['time']

    return dict
