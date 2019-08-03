
'''
    
    Copy this file into the dir of run.
    Requires 'chfoil_default5_norm.txt' in directory of run - text file containing T0,n0, lnLambda
    This file will make  a directory and put 2D image map pictures of the variables at each time step in it.
    
'''
import os,site,getpass,re,sys
import subprocess as sp
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import re
import matplotlib.ticker as ticker
userid = getpass.getuser()
site.addsitedir('/Users/' + userid + '/Dropbox/IMPACT_dir/chfoil_d5/MODULES')
import chfoil_module as cf

#import IMPACT_mod as IM

import matplotlib as mpl
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.8
mpl.rcParams['xtick.major.width'] = 1.8
mpl.rcParams['xtick.major.pad'] = 3
#mpl.rcParams['axes.xmargin'] = 0.5
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 3.5
mpl.rcParams['ytick.minor.width'] = 2.5
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18

mpl.rcParams['legend.fontsize']=26
mpl.rcParams['legend.frameon']=False
mpl.rcParams['axes.titlesize']=24
mpl.rcParams['axes.labelsize']='x-large'
mpl.rcParams['axes.linewidth']=2
#mpl.rcParams['axes.labelpad']= 8

mpl.rcParams['figure.subplot.wspace']=0.55
mpl.rcParams['figure.subplot.hspace']=0.4
mpl.rcParams['figure.subplot.bottom']=0.05
#mpl.rcParams['figure.subplot.top']=0.05
mpl.rcParams['figure.dpi']=120
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rcParams['figure.figsize'] = 30,25
asp = 'auto'#0.50

q_e = 1.602e-19
m_e = 9.11e-31
m_p = 1.67e-27
k_b = 1.38e-23
epsilon0 = 8.854e-12
#---------------
#n0 = (1.87e22) # [cm^-3]
#T0 = 850 #eV
# v_te = 1.729e7 # ms^-1
# tau_ei = 1.0578e-14 # s
# nu_ei = 9.45e-13 # s^-1
# lambda_mfp = 1.829e-7 # m


cwdpath = os.getcwd()
path2 = cwdpath.split('/')[:-1]
init_path = '/'.join(path2)
print 'init_path = ', init_path
cd5 = cf.conv_factors_premag()
Z0 = cd5.Z
T0 = cd5.T0
n0 = cd5.n0
logLambda = cd5.logLambda

v_te = cd5.v_te # ms^-1
tau_ei = cd5.tau_ei # s
nu_ei = 1.0/tau_ei # s^-1
lambda_mfp = cd5.lambda_mfp # m




SI_on = True
mult = 6
yi = 4
if SI_on:
    xstep_factor = lambda_mfp*1e6
    tstep_factor = tau_ei*1e12
    xlab = r'$x$ [ $\mu m$ ]'
    ylab = r'$y$ [ $\mu m$ ]'
    leg_title = r'time [ $ps$ ]'

else:
    xstep_factor = 1.0
    tstep_factor = 1.0
    xlab = r'\textit{x [ $\lambda_{mfp}$ ]}'
    xlab = r'\textit{x [ $\lambda_{mfp}$ ]}'
    leg_title = r'\text{$time$ [t$_{col}$ ]}'
	

#-----------------------------------------------------------------------
class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
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
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

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
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

#-----------------------------------------------------------------------
class MyDict(dict):
    pass

#-----------------------------------------------------------------------
def construct_fname(path,fprefix,var,time):
    
    if var=='fo' or var=='fxX' or var=='fxY' or var=='fyX' or var=='fyY':
        suffix = '.xyv'
    else:
        suffix = '.xy'
        
    #fprefix = 'thydro_hi'
    #time = '00'
    fname = path + '/' + fprefix +'_' + var + '_' + time + suffix
    return fname

#-----------------------------------------------------------------------
def get_startline(fname):
    '''
        get the line where the data starts
    '''
    f = open(fname)
    f_list = f.readlines()
    for i in range(len(f_list)):
        if re.match(r'\n',f_list[i]):
            #print 'found match on line: ', i
            out_line = i
        else:
            out_line = 10
            continue
    return out_line
    
def list_to_float(input):
    list = input.split()
    arr = np.zeros((len(list)),dtype=float)
    arr[:] = list[:]
    return arr
#-----------------------------------------------------------------------    
def fpg_get_info(fname):
    '''
        Gets the IMPACT header info
        dict = fpg_get_info(fname)
    '''
    info = open(fname, 'r')
    data = info.readlines()
    mat_info = data[:10]
    mat = np.loadtxt(fname,skiprows=10)

    print mat_info[0]
    time = float(mat_info[1])
    ndims = int(mat_info[2])
    
    dict = MyDict()

    dict['time'] = time
    dict['ndims'] = ndims
    if ndims ==2:

        ny = int(mat_info[3])        
        y_grid = list_to_float(mat_info[4])
        nx = int(mat_info[5])
        x_grid = list_to_float(mat_info[6])
        dim_array = [ny,nx]
        
    else:
        
        nv = int(mat_info[3])
        v_grid = list_to_float(mat_info[4])
        dict['nv'] = nv
        dict['v_grid'] = v_grid     
        ny = int(mat_info[5])
        y_grid = list_to_float(mat_info[6])
    
        nx = int(mat_info[7])
        x_grid = list_to_float(mat_info[8])
        dim_array = [nv,ny,nx]
    
    dict['nx'] = nx
    dict['ny'] = ny
    dict['x_grid'] = x_grid
    dict['y_grid'] = y_grid

    return dict

#-----------------------------------------------------------------------
def get_index(ny,nx):
       index = ny + 2 + 1
       return index

#-----------------------------------------------------------------------       
def get_time(fname):
       f = open(fname)
       data = f.readlines()
       t = float(data[1])
       return t    
def extract_power(x):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return b
    
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)
def fmt_ord(x, pos):

    return r'${}$'.format(x)
           
def calc_norms(var,sample =0.0):
    '''
        norm_const, ylab = calc_norms(var)
    '''
    #print var, np.shape(var)
    c_fmt = '%3.2f'
    
    if SI_on:
        if var== 'Cx' or var=='Cy':
            
            norm_const = v_te*1e-3
            title = r'$' + var[0] + '_' + var[1] + r'$ [ $  kms^{-1}$ ]'
            
        elif var=='n':
            norm_const = n0*1e-22
            title = r'$n_e$ [ $10^{22}$ $cm^-3$ ]'

        elif var=='Te':
            norm_const = 2.0*T0
            title = r'$T_e$ [ $eV$ ]'
            c_fmt = '%3.0f'
        
        elif var=='Bz':
            #power = extract_power(sample, pos)
            norm_const = (m_e/(q_e*tau_ei))
            title = r'$B_z$ [ $T$ ]'
            
        elif var=='wt':
            power = extract_power(sample)
            mod = r'$ 10^{' +str(power)  + '} $'
            
            norm_const = 1.0*(10**(-power))
            
            title = mod + r'$\omega \tau_{ei}$'
            c_fmt = '%1.1f'
            
        elif var[0]=='E':
            print 'DOING E-field - min sample = ', sample
            norm_const = (lambda_mfp/(tau_ei**2))*(m_e/q_e)
            power = extract_power(norm_const*sample)
            print ' power extracted = ', power
            c_fmt = '%1.1f'
            mod = r'$ 10^{' +str(power)  + '}$'
            norm_const = norm_const*(10**-power)
            
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$V/m$ ]'
        
        elif var[0] =='q':
            
            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0*(10**-power)
            mod = r'$ 10^{' +str(power)  + '} $'
            #c_fmt = '%1.0e'
            #norm_const = m_e*(v_te**3)*n_04
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' + mod + '$q_0$ ]'# + r' [ $Jms^{-1}$ ]'

        elif var[0] =='j':
            c_fmt = '%1.1f'
            power = extract_power(sample)
            norm_const = 1.0*(10**-power)
            mod = r'$ 10^{' +str(power)  + '} $'
            #norm_const = 1.0#q_e*n0*v_te
            
            title = r'$' + var[0] + '_' + var[1] + r'$ [ ' +mod + '$j_0$ ]'#r' [ $C/m^2 s$ ]'    
    else:
        norm_const = 1.0
        title = var
    
    print 'var: ', var, title        
    return norm_const, title, c_fmt
#-----------------------------------------------------------------------
if __name__=="__main__":
    basedir = "$BASEDIR"
    rundir = "$RUN"

    path = os.getcwd() #'chfoil_default5_heat11_2D3_noB--cx1/'#
    os.system('gunzip -vf ' +  path + '/*.gz')
    #var_list = ['ExX','EyY','qxX','qyY','Bz','wt']
    var_list = ['n','Te','Cx','Cy','Bz','ExX','ExY','EyX','EyY','qxX','qxY','qyX','qyY','jxX','jxY','jyX','jyY']
    #var_list = ['n','Te','Cx','Cy','Bz','ExX','EyY','qxX','qyY','jxX','jyY']

    #var_list = ['n','Te']
    
    fsize = 17
    params = {
             'axes.labelsize': fsize,
              'text.fontsize': fsize+4,
              'figure.titlesize': fsize,
              'legend.fontsize': fsize,
             'text.usetex': False,
             'xtick.minor.pad':10,
             'figure.subplot.left': 0.1,
            'figure.subplot.right': 0.96,   # the right side of the subplots of fig
            'figure.subplot.top': 0.93,
            'figure.subplot.wspace': 0.4,  
            'figure.subplot.hspace': 0.25}


    #t_list = ['00', '01', '02']
    t_list = []
    tc_list = []
    #print var_list
    
    list_files = os.listdir(path)
    for fname in list_files:
        s = re.search('(?P<fpre>.*?)_' + var_list[0] + '_(?P<time>\d\d).xy.*?',fname)
        if s:
            t_l = s.group('time')
            t_list.append(t_l)
            t_col = get_time(path +'/' +fname)
            tc_list.append(t_col)
            print t_list
            fprefix = s.group('fpre')
    fprefix = path.split('/')[-1]        
    fname1 = 'thydro_hi--andy_wt_00.xy'
    fname2 = 'thydro_hi--andy_wt_00.xy'
    
    
    # lineouts
    index = 3
    print 't_list: \n',t_list
    #out_l = get_startline(fname1)
    #data1 = np.loadtxt(fname1,skiprows=out_l)
    t_points = [0,4,5,6,7,8]
    #subplots_adjust(hspace=0.001)
    number_of_subplots=len(var_list)
    leg_list = []
    x_index = 3
    y_index = number_of_subplots/x_index
    print 'len var(list) = ', len(var_list),' X_INDEX Y_INDEX = ', x_index,y_index
    if (x_index*y_index)< len(var_list):
        x_index+=1
    
    #x_index,y_index = 2,1
    #t_list = ['00','03','05','06','07','08']
    # select time you want to plot
    save_path = 'twodim_temp_pert'
    os.system('mkdir ' + save_path)
    
    #t_list = ['15']
    #-------- time loop 
    for tt in range(len(t_list)):
        time = t_list[tt]
        fig = plt.figure()
        #--------------sub plot loop---------------------------
        for i,v in enumerate(xrange(number_of_subplots)):
            v = v+1
            ax1 = fig.add_subplot(y_index,x_index,v)
            var = var_list[i]
            
            norm_const, c_title,c_fmt = calc_norms(var)
            print 'time: ', time, 'tt: ', tt
            fname = construct_fname(path,fprefix,var,time)
            time_f= get_time(fname)
            out_l = get_startline(fname)
            leg_list.append(get_time(fname))
            data = np.loadtxt(fname,skiprows=out_l)
            data = np.where(abs(data)<=1e-99,0.0,data)
            print 'VAR: ', var, 'shape: ', np.shape(data)
            print ' xstep_factor = ', xstep_factor
            #neg_list = {'ExX','EyY','qxX','qyY'}
            neg_list = {}
            min_data = np.min(data)
            max_data = np.max(data)
            norm_const, c_title,c_fmt = calc_norms(var, min_data)
            print ' min_data = ', min_data, 'max_data = ', max_data
            if var in neg_list:
                norm = MidPointNorm(midpoint=0.0)
                colormap='RdBu_r'#plt.cm.seismic
                #if var == 'Bz':
                #norm = colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                #                         vmin=-1.0, vmax=40.0),
            elif min_data == max_data:
                norm = None
            else:
                middle = 0.5*(min_data+max_data)*norm_const
                #if var=='Te':
                #    middle = 0.332*norm_const
                norm = MidPointNorm(midpoint=middle)

                #norm = 0.3
                colormap='RdBu_r'
                #colormap=cm.Spectral_r#'RdBu_r'#'jet'#cm.seismic#'RdBu_r'
            # headers on Cx and Cy files don't work yet
            if (var!='Cx') and (var!='Cy'):
                dict1 = fpg_get_info(fname)
        
            #dict1 = fpg_get_info(fname)
            #lims = [dict1['x_grid'][0],dict1['x_grid'][-1],dict1['y_grid'][0],dict1['y_grid'][-1]]
            if len(dict1['x_grid']) != np.shape(data)[0]:
                x_c_temp =dict1['x_grid'][1:-1]*xstep_factor
                print 'Mod SHAPE data: ', np.shape(data), np.shape(x_c_temp)
            else:
                x_c_temp =dict1['x_grid']*xstep_factor
                print ' SHAPE data: ', np.shape(data), np.shape(x_c_temp)
            x_c_grid = x_c_temp[::-1]
            X,Y = np.meshgrid(dict1['y_grid']*xstep_factor,x_c_grid)
            lims = [dict1['y_grid'][0]*xstep_factor,dict1['y_grid'][-1]*xstep_factor,dict1['x_grid'][0]*xstep_factor,dict1['x_grid'][-1]*xstep_factor]
            data_dev = cf.get_U_dev_abs(data)
            im = ax1.imshow(np.log(np.abs(data_dev)),
                            cmap=colormap,
                            aspect=asp,
                            extent=lims)
        
            #CS= ax1.contour(X,Y,data*norm_const,colors=('k'),extent=lims)
            title = r'$ ' + var +'$'
            #ax1.set_title(title)
            ax1.set_xlabel(ylab)
            ax1.set_ylabel(xlab)
        
            if c_fmt[-1] == 'e':
                fmt_cbar = ticker.FuncFormatter(fmt)
            else:
                fmt_cbar = ticker.FuncFormatter(fmt_ord)
        
            plt.colorbar(im,
                            shrink=0.9,
                            label=var)
            #ax1.colorbar()
            #cbar = plt.colorbar(ax=ax1)
        plt.suptitle(t_list[0] + ' time: ' + str(time_f) + ' t_col')
        rcParams.update(params)
        save_name = save_path + '/' + fprefix + time + '.png'
        fig.savefig(save_name, dpi=120)
        #plt.clf()
        print '----- saved time ', time,' as ',save_name
    #plt.colorbar()
    #plt.legend(leg_list,ncol=ncolly,bbox_to_anchor=[3.0, .4], loc='center right')
    #plt.tight_layout()
    #plt.show()
    plt.clf()
    plt.close(fig)

