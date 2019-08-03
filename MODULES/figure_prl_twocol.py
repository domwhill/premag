import numpy as np, getpass
import matplotlib as mpl
import pylab
userid = getpass.getuser()
mpl.use('pgf')

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
    fig_width_pt = 246.0#483.497#246.0 # column width#510.0 full page#483.497#469.755                  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale*scale_width    # width in inches
    fig_height = fig_width*golden_mean*scale_ratio              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

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

maxlinewidth=1.0
maxwidth = 1
minwidth = 0.7
minorsize=1
majorsize=2
'''
    plt.rcParams.update(
        {'lines.linewidth':1,
         'legend.fontsize':8,
         'axes.titlesize':8,
         'axes.linewidth':1,
         'axes.labelsize':8,
         'axes.labelpad':1,
         'xtick.major.size':2,
         'ytick.major.size':2,
         'xtick.major.width':1,
         'ytick.major.width':1,
         'ytick.minor.width':0.7,
         'ytick.minor.size':1,
         'xtick.labelsize':7,
         'ytick.labelsize':7})
'''


pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    #"axes.labelsize": 12,               # LaTeX default is 10pt font.
    #"text.fontsize": 12,
    #"legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    #"xtick.labelsize": 10,
    #"ytick.labelsize": 10, 
    "axes.labelsize": 7,            # LaTeX default is 10pt font.
    "axes.linewidth": 1.0,
    "text.fontsize": 10,
    "legend.fontsize": 6,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.major.size": majorsize,
    "ytick.major.size": majorsize,
    "xtick.minor.size": minorsize,
    "ytick.minor.size": minorsize,
    "xtick.major.width": maxwidth,
    "ytick.major.width": maxwidth,
    "xtick.minor.width": minwidth,
    "ytick.minor.width": minwidth,
    "lines.linewidth": maxlinewidth,
    'figure.subplot.bottom': 0.15,
    'figure.subplot.top': 0.95,
    'figure.subplot.left': 0.15,
    'figure.subplot.right': 0.95,
    #'figure.subplot.bottom': 0.12,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "figure.dpi":600,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}"
        ],
    'text.latex.preamble': [r'\usepackage{siunitx}']
    }

pgf_with_latex_2scale = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    #"axes.labelsize": 12,               # LaTeX default is 10pt font.
    #"text.fontsize": 12,
    #"legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    #"xtick.labelsize": 10,
    #"ytick.labelsize": 10, 
    "axes.labelsize": 7,            # LaTeX default is 10pt font.
    "axes.linewidth": 1.0,
    "text.fontsize": 10,
    "legend.fontsize": 5,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    #"xtick.major.size": 20,
    #"ytick.major.size": 20,
    #"xtick.minor.size": 10,
    #"ytick.minor.size": 10,
    "xtick.major.width": maxwidth,
    "ytick.major.width": maxwidth,
    "xtick.minor.width": minwidth,
    "ytick.minor.width": minwidth,
    "lines.linewidth": maxlinewidth,
    'figure.subplot.bottom': 0.15,
    'figure.subplot.top': 0.95,
    'figure.subplot.left': 0.15,
    'figure.subplot.right': 0.85,
    #'figure.subplot.bottom': 0.12,
    "figure.figsize": figsize(1.0),     # default fig size of 0.9 textwidth
    "figure.dpi":600,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}"
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

def newfig_generic(width,scale_width=1.0,scale_ratio=1.0):
    '''
      fig, ax = newfig(width)  
    '''
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width,scale_width,scale_ratio),forward=True)
    print ' fig size = ', figsize(width,scale_width)
    #pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    #pylab.axes([0.125,0.3,0.95-0.125,0.95-0.3])

    return fig
    

def newfig_generic_twinx(width,scale_width=1.0,scale_ratio=1.0):
    '''
        newfig_generic but with increased right hand spacing 
        to account for y axis title labels on both 
        left and right hand side y axes. Used for 1 subplot only.
      fig = newfig_generic_twinx(width)  
    '''
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    plt.subplots_adjust(left=0.16,right=0.85,bottom=0.16)
    fig.set_size_inches(figsize(width,scale_width,scale_ratio),forward=True)
    print ' fig size = ', figsize(width,scale_width)
    #pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    #pylab.axes([0.125,0.3,0.95-0.125,0.95-0.3])

    return fig
    

def newfig_generic_2yscale(width,scale_width=1.0,scale_ratio=1.0):
    '''
      fig, ax = newfig(width)  
      has extra rspace for ysclae
    '''
    mpl.rcParams.update(pgf_with_latex_2scale)
    import matplotlib.pyplot as plt
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width,scale_width,scale_ratio),forward=True)
    print ' fig size = ', figsize(width,scale_width)
    #pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    #pylab.axes([0.125,0.3,0.95-0.125,0.95-0.3])

    return fig

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))

def savefig_thesis(filename):
    pic_path = thesis_picpath()
    fname = pic_path + '/' + filename
    plt.savefig('{}.pgf'.format(fname))
    plt.savefig('{}.pdf'.format(fname))
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
