'''
# 3 subplots side by side --> (ms wave format)
fig = fprl.newfig_generic_2yscale(1.0,scale_width=1.2,scale_ratio=0.6)
plt.subplots_adjust(wspace=0.15,bottom=0.2)

#---> else with log axes
        fig = fprl.newfig_generic_2yscale(1.0,scale_width=1.2,scale_ratio=0.5)
        #plt.subplots_adjust(wspace=0.55,left=0.2,bottom=0.2,right=0.95)
        plt.subplots_adjust(wspace=0.45,left=0.1,bottom=0.15,right=0.95)
<---
'''

import numpy as np, getpass,sys,re
import matplotlib as mpl
from pylab import *
userid = getpass.getuser()
mpl.use('pgf')
class plotting_params():
    '''
    cmap = fprl.plotting_params.lineouts_cmap
    '''
    lineouts_cmap = cm.plasma

def thesis_picpath():
    '''
        pic_path = thesis_picpath()
    '''
    pic_path = '/Users/' + userid +'/Dropbox/PhD/Reports/THESIS/THESIS_TEX/PICS_THESIS'
    return pic_path

def figsize(scale,scale_width=1.0,scale_ratio=1.0):
    '''
     fig_size = figsize(scale)
    '''
    fig_width_pt = 483.497#483.497#246.0 # column width#510.0 full page#483.497#469.755                  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale*scale_width    # width in inches
    fig_height = fig_width*golden_mean*scale_ratio              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

def figsize_custom(width_scale,height_scale):
    '''
     fig_size = figsize(scale)
    '''
    fig_width_pt = 483.497*width_scale#469.755                  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt    # width in inches
    fig_height = fig_width*height_scale # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 12,               # LaTeX default is 10pt font.
    "text.fontsize": 12,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    'figure.subplot.bottom': 0.12,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "lines.linewidth": 1.5,
    "axes.linewidth": 1.5,
    'xtick.major.width':1.5,
    'ytick.major.width':1.5,
    'xtick.minor.width':1.0,
    'ytick.minor.width':1.0,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}",
        ],
    'text.latex.preamble': [r'\usepackage{siunitx}']

    }
pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 12,               # LaTeX default is 10pt font.
    "text.fontsize": 12,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    'figure.subplot.bottom': 0.12,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "lines.linewidth": 1.0,
    "axes.linewidth": 1.0,
    'xtick.major.width':1.0,
    'ytick.major.width':1.0,
    'xtick.minor.width':0.5,
    'ytick.minor.width':0.5,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}",
        ],
    'text.latex.preamble': [r'\usepackage{siunitx}']

    }
maxlinewidth=1.0
maxwidth = 1
minwidth = 0.7
minorsize=2
majorsize=3

pgf_with_2subplots = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 12,               # LaTeX default is 10pt font.
    "axes.linewidth": maxlinewidth,
    "xtick.major.size": majorsize,
    "ytick.major.size": majorsize,
    "xtick.minor.size": minorsize,
    "ytick.minor.size": minorsize,
    "xtick.major.width": maxwidth,
    "ytick.major.width": maxwidth,
    "xtick.minor.width": minwidth,
    "ytick.minor.width": minwidth,
    "lines.linewidth": maxlinewidth,
    "text.fontsize": 12,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    'figure.subplot.bottom': 0.12,
    'figure.subplot.wspace': 0.35,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "lines.linewidth": 1.5,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}",
        ],
    'text.latex.preamble': [r'\usepackage{siunitx}']

    }

pgf_with_latex_presentation = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 18,               # LaTeX default is 10pt font.
    "text.fontsize": 18,
    "legend.fontsize": 14,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    'figure.subplot.bottom': 0.12,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "lines.linewidth": 2.0,
    "axes.linewidth": 2.0,
    'xtick.major.width':2.0,
    'ytick.major.width':2.0,
    'xtick.minor.width':1.5,
    'ytick.minor.width':1.5,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}",
        ],
    'text.latex.preamble': [r'\usepackage{siunitx}']
    }

mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width):
    '''
      fig, ax = newfig(width)  
    '''
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    ax.tick_params(which='both',direction='in')
    

    return fig, ax
def newfig_generic(width):
    '''
      fig = newfig_generic(width)  
    '''
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    #ax = fig.add_subplot(111)
    

    return fig
def newfig_generic_presentation(width_scale,height_scale=1.0):
    '''
      fig = newfig_generic(width)  
    '''
    mpl.rcParams.update(pgf_with_latex_presentation)
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize_custom(width_scale,height_scale))
    #ax.tick_params(which='both',direction='in')

    #ax = fig.add_subplot(111)
    

    return fig
def newfig_generic_2yscale(width,scale_width=1.0,scale_ratio=0.8,clearon=True):
    '''
      fig, ax = newfig(width)  
      has extra rspace for ysclae
    '''
    mpl.rcParams.update(pgf_with_2subplots)
    import matplotlib.pyplot as plt
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width,scale_width,scale_ratio),forward=True)
    print ' fig size = ', figsize(width,scale_width)

    return fig

def newfig_generic_3yscale(width,scale_width=1.0,scale_ratio=0.5,clearon=True):
    '''
      fig, ax = newfig(width)  
      has extra rspace for ysclae
    '''
    mpl.rcParams.update(pgf_with_2subplots)
    import matplotlib.pyplot as plt
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width,scale_width,scale_ratio),forward=True)
    print ' fig size = ', figsize(width,scale_width)

    return fig



def newfig_custom(width_scale,height_scale):
    '''
      fig= newfig_custom(width_scale,height_scale) 
    '''
    plt.clf()
    fig = plt.figure(figsize=figsize_custom(width_scale,height_scale))
    #ax = fig.add_subplot(111)
    

    return fig

def annotate_axis(ax,lett='(a)',dx_mult=1.0,dy_mult=1.0,fontsize=0):
    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    nunits =20.0
    dx = (xmax - xmin)/nunits
    dy = (ymax-ymin)/nunits
    print 'x = ',xmin-dx,'y = ', ymax
    x_coord = xmin-5.0*dx*dx_mult
    y_coord = ymax - dy*dy_mult
    if fontsize!=0:
        ax.text(x_coord,y_coord,lett,fontsize=fontsize)

    else:
        ax.text(x_coord,y_coord,lett)
    return

def set_ticks(ax):
    ax.tick_params(which='both',direction='in')


def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))

def savefig_thesis(filename):
    pic_path = thesis_picpath()
    fname = pic_path + '/' + filename
    print 'saving as ===', fname
    print ' sys arg 0 = ', sys.argv
    
    f = open(pic_path + '/pic_path.log','r')
    lines = f.readlines()
    line_present = False
    for line in lines:
        r=re.search(fname, line)
        if r:
            il = lines.index(line)
            line_present = True
            f.close()
            break
    
    if not line_present:
        f.close()
        f =open(pic_path + '/pic_path.log','a') 
        f.write(sys.argv[0] + ' ' + fname)
        f.close()
    
    plt.savefig('{}.pgf'.format(fname))
    plt.savefig('{}.pdf'.format(fname))
    
def fix_ax(ax):
    ax.tick_params(which='both',direction='in')
    ax.grid(color='0.5',linestyle='-') 
'''

# Simple plot
fig, ax  = newfig(0.6)

def ema(y, a):
    s = []
    s.append(y[0])
    for t in range(1, len(y)):
        s.append(a * y[t] + (1-a) * s[t-1])
    return np.array(s)
y = [0]*200
y.extend([20]*(1000-len(y)))
s = ema(y, 0.01)

ax.plot(s)
ax.set_xlabel('X Label')
ax.set_ylabel('EMA')

savefig('ema')
'''