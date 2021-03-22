import numpy as np, getpass
import matplotlib as mpl
import pylab
import matplotlib.gridspec as GS
from matplotlib import ticker
import re
import pdb
userid = getpass.getuser()
mpl.use('pgf')


def thesis_picpath():
    '''
        pic_path = thesis_picpath()
    '''
    pic_path = '/Users/' + userid + '/Dropbox/PhD/Reports/THESIS/THESIS_TEX/PICS_THESIS'
    return pic_path


def figsize(scale, scale_width=1.0, scale_ratio=1.0):
    '''
     fig_size = figsize(scale)
    '''
    fig_width_pt = 246.0    #483.497#246.0 # column width#510.0 full page#483.497#469.755                  # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0 / 72.27    # Convert pt to inch
    golden_mean = (np.sqrt(5.0) - 1.0) / 2.0    # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt * inches_per_pt * scale * scale_width    # width in inches
    fig_height = fig_width * golden_mean * scale_ratio    # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size


def annotate_axis(ax, lett='(a)', dx_mult=1.0, dy_mult=1.0, fontsize=0):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    nunits = 20.0
    dx = (xmax - xmin) / nunits
    dy = (ymax - ymin) / nunits
    print 'x = ', xmin - dx, 'y = ', ymax
    x_coord = xmin - 5.0 * dx * dx_mult
    y_coord = ymax - dy * dy_mult
    if fontsize != 0:
        ax.text(x_coord, y_coord, lett, fontsize=fontsize)

    else:
        ax.text(x_coord, y_coord, lett)
    return


maxlinewidth = 1.0
maxwidth = 1
minwidth = 0.7
minorsize = 1
majorsize = 2
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

pgf_with_latex = {    # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",    # change this if using xetex or lautex
    "text.usetex": True,    # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],    # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 7,    # LaTeX default is 10pt font.
    "axes.linewidth": 1.0,
    #"text.fontsize": 10,
    "legend.fontsize": 6,    # Make the legend/label fonts a little smaller
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
    "figure.figsize": figsize(1.0),    # default fig size of 0.9 textwidth
    "figure.dpi": 600,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",    # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}"
    ],
    'text.latex.preamble': [r'\usepackage{siunitx}']
}

pgf_with_latex_2scale = {    # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",    # change this if using xetex or lautex
    "text.usetex": True,    # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],    # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 7,    # LaTeX default is 10pt font.
    "axes.linewidth": 1.0,
    #"text.fontsize": 10,
    "legend.fontsize": 5,    # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.major.width": maxwidth,
    "ytick.major.width": maxwidth,
    "xtick.minor.width": minwidth,
    "ytick.minor.width": minwidth,
    "lines.linewidth": maxlinewidth,
    'figure.subplot.bottom': 0.15,
    'figure.subplot.top': 0.95,
    'figure.subplot.left': 0.15,
    'figure.subplot.right': 0.85,
    "figure.figsize": figsize(1.0),    # default fig size of 0.9 textwidth
    "figure.dpi": 600,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",    # plots will be generated using this preamble
        r"\usepackage[detect-all,locale=DE]{siunitx}"
    ],
    'text.latex.preamble': [r'\usepackage{siunitx}']
}
pgf_with_latex_nofontsize = {    # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",    # change this if using xetex or lautex
    "text.usetex": True,    # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],    # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    # default fig size of 0.9 textwidth
    "figure.dpi": 600,
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",    # plots will be generated using this preamble
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
    ax.tick_params(which='both', direction='in')

    return fig, ax


def newfig_generic(width, scale_width=1.0, scale_ratio=1.0):
    '''
      fig, ax = newfig(width)  
    '''
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width, scale_width, scale_ratio), forward=True)
    print ' fig size = ', figsize(width, scale_width)
    #pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    #pylab.axes([0.125,0.3,0.95-0.125,0.95-0.3])

    return fig


def newfig_generic_twinx(width, scale_width=1.0, scale_ratio=1.0):
    '''
        newfig_generic but with increased right hand spacing 
        to account for y axis title labels on both 
        left and right hand side y axes. Used for 1 subplot only.
      fig = newfig_generic_twinx(width)  
    '''
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    plt.subplots_adjust(left=0.16, right=0.85, bottom=0.16)
    fig.set_size_inches(figsize(width, scale_width, scale_ratio), forward=True)
    print ' fig size = ', figsize(width, scale_width)
    #pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    #pylab.axes([0.125,0.3,0.95-0.125,0.95-0.3])

    return fig


def newfig_generic_2yscale(width, scale_width=1.0, scale_ratio=1.0):
    '''
      fig, ax = newfig(width)  
      has extra rspace for ysclae
    '''
    mpl.rcParams.update(pgf_with_latex_2scale)
    import matplotlib.pyplot as plt
    #plt.clf()
    #fig = plt.figure(figsize=figsize(width,scale_width))
    fig = plt.figure()
    fig.set_size_inches(figsize(width, scale_width, scale_ratio), forward=True)
    print ' fig size = ', figsize(width, scale_width)

    return fig


def set_ax(ax, lim):
    ax.set_xlim(lim[0], lim[1])
    ax.set_ylim(lim[2], lim[3])


def clear_xax(ax):
    ax.set_xticklabels([])
    ax.set_xticks([])


def clear_yax(ax):
    ax.set_yticklabels([])
    ax.set_yticks([])


def sci_fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def fmt(x, pos):
    return r'%1.1e' % x


def custom_fmt(x, pow_in=0):
    '''
        This function takes first 2 sigfig of number and returns them as a
        1.2 where 1 = first number and 2 = second number
    :param x:
    :param pos:
    :return:
    '''
    init_f, pow = ('%1.1e' % (x)).split('e')
    if pow_in != 0:
        init_f = '%1.1f' % (float(init_f) * 10**(pow_in - int(pow)))
    return init_f


def custom_fmt_lab(x, lab):
    '''
        This function takes first 2 sigfig of number and returns them as a
        1.2 where 1 = first number and 2 = second number
    :param x:
    :param pos:
    :return:
    '''
    pow = int(('%1.1e' % (x)).split('e')[-1])
    #pdb.set_trace()
    a = re.search(r'\[\$', lab)
    b = re.search(r'\[', lab)

    if a:
        delim = '[$'
    elif b:
        delim = '['
    else:
        return 0, lab
    lab_sep = lab.split(delim)
    if pow != 0:
        mod_lab = '10 ^ {%i}' % (pow)
    else:
        mod_lab = ''
    lab_out = lab_sep[0] + delim + mod_lab + lab_sep[-1]
    return pow, lab_out


def newfig_generic_3lineouts(width=1.0):
    fig = newfig_generic_2yscale(1.2, scale_width=1.0, scale_ratio=0.5)
    return fig


def newfig_4implot(grid, data_in, out_file='im.png', cmap='RdBu_r', vmin=0, vmax=0, **kwargs):
    '''

    '''
    data1 = data_in['data1']
    data2 = data_in['data2']
    data3 = data_in['data3']
    data4 = data_in['data4']

    x_grid = grid['x_grid']
    y_grid = grid['y_grid']

    gs1 = GS.GridSpec(2, 3, width_ratios=[1.0, 1.0, 0.04])
    fig = newfig_generic(1.0)
    plt.subplots_adjust(hspace=0.12, wspace=0.12, right=0.8)

    ax1 = plt.subplot(gs1[0, 0])
    ax2 = plt.subplot(gs1[0, 1])
    ax3 = plt.subplot(gs1[1, 0])
    ax4 = plt.subplot(gs1[1, 1])
    cax = plt.subplot(gs1[:, 2])

    lims = [x_grid[0], x_grid[-1], y_grid[0], y_grid[-1]]
    if vmin == vmax:
        vmin = np.min(data4)
        vmax = np.max(data4)
    im = ax1.imshow(data1, aspect='auto', extent=lims, cmap=cmap, vmin=vmin, vmax=vmax)
    im = ax2.imshow(data2, aspect='auto', extent=lims, cmap=cmap, vmin=vmin, vmax=vmax)
    im = ax3.imshow(data3, aspect='auto', extent=lims, cmap=cmap, vmin=vmin, vmax=vmax)
    im = ax4.imshow(data4, aspect='auto', extent=lims, cmap=cmap, vmin=vmin, vmax=vmax)

    if 'lim' in data_in:
        lim = data_in['lim']
        set_ax(ax1, lim)
        set_ax(ax2, lim)
        set_ax(ax3, lim)
        set_ax(ax4, lim)

    if 'lab' in data_in:
        lab = data_in['lab']
    else:
        lab = ''

    c2 = fig.colorbar(im, cax=cax, ax=cax, format=ticker.FuncFormatter(sci_fmt), label=lab)
    if 'ylab' in data_in:
        ylab = data_in['ylab']
        ax1.set_ylabel(ylab)
        ax3.set_ylabel(ylab)
    clear_yax(ax2)
    clear_yax(ax4)

    clear_xax(ax1)
    clear_xax(ax2)
    if 'xlab' in data_in:
        xlab = data_in['xlab']
        ax3.set_xlabel(xlab)
        ax4.set_xlabel(xlab)
    return fig, ax1, ax2, ax3, ax4, cax


def sci_fmt_colorbar(fig, im, ax, cax='', lab=[]):
    '''
        NOTE: This function CANNOT be used with c2.update_ticks() (afterwards) as
        this will reset the tick labels to the previous format
        :param fig:
        :param im:
        :param ax:
        :param cax:
        :param lab:
        :return:
        '''
    # --- get clabels
    if cax == '':

        c2 = fig.colorbar(
            im,
            ax=ax,    # ticker.FuncFormatter(custom_fmt),
            label=lab)
    else:
        c2 = fig.colorbar(im, cax=cax, ax=ax, format=ticker.FuncFormatter(sci_fmt), label=lab)

    return c2


def custom_colorbar(fig, im, ax, cax='', lab=[]):
    '''
    NOTE: This function CANNOT be used with c2.update_ticks() (afterwards) as
    this will reset the tick labels to the previous format
    :param fig:
    :param im:
    :param ax:
    :param cax:
    :param lab:
    :return:
    '''
    # --- get clabels
    if cax == '':

        c2 = fig.colorbar(
            im,
            ax=ax,    #ticker.FuncFormatter(custom_fmt),
            label=lab)
    else:
        c2 = fig.colorbar(
            im,
            cax=cax,
            ax=ax,    #ticker.FuncFormatter(custom_fmt),
            label=lab)

    # --> now swap y ticks
    conv_tick = lambda a: float(re.sub('\$', '', a.get_text()))
    yticks = c2.ax.get_yticklabels()
    new_yticks = []
    pow_max = 0
    c_lab_max = lab

    tick_locator = ticker.MaxNLocator(nbins=3)
    c2.locator = tick_locator
    c2.update_ticks()

    for tick in yticks:
        tick_float = conv_tick(tick)
        pow, c_lab = custom_fmt_lab(tick_float, lab)
        if np.abs(int(pow)) > pow_max:
            pow_max = int(pow)
            c_lab_max = c_lab
    for tick in yticks:
        tick_float = conv_tick(tick)
        new_tick = '$' + custom_fmt(tick_float, pow_max) + '$'
        new_yticks.append(new_tick)
    c2.ax.set_yticklabels(new_yticks)
    c2.set_label(c_lab_max)

    return c2


def get_lab_power(c2):
    label_list = []
    pwd = []
    for i in c2.ax.get_yticklabels():
        a = (i.get_text())
        label_list.append(a)
        asearch = re.search('10\^\{(?P<pw>\d+)\}', i)
        if asearch:
            pwd.append(int(a.group('pw')))

    return pwd


def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))


def savefig_thesis(filename):
    pic_path = thesis_picpath()
    fname = pic_path + '/' + filename
    plt.savefig('{}.pgf'.format(fname))
    plt.savefig('{}.pdf'.format(fname))


def fmt_axticker(ax, axis_type='y', fmt='%i'):
    '''
        set the format for the labels 
    '''
    if axis_type == 'y':
        ax_in = ax.yaxis
    else:
        ax_in = ax.xaxis
    ax_in.set_major_formatter(ticker.FormatStrFormatter(fmt))
    return


def set_ylim_max(ax_in, grid, data, y_mult=[1.0, 1.0], xlim=[-5.0, 20.0]):
    xmin, xmax = xlim[0], xlim[1]
    ymin = np.min(data[(grid < xmax) * (grid >= xmin)])
    ymax = np.max(data[(grid < xmax) * (grid >= xmin)])
    dy = 0.05 * np.abs(ymax - ymin)
    if ymin != ymax:

        ax_in.set_ylim(ymin - dy * y_mult[0], ymax + dy * y_mult[1])
    else:
        print('error ymin = ymax data all the same!')
    print('ymin = %4.4e ymax = %4.4e ' % (ymin, ymax))

    return


def get_ylim_max(grid, data, y_mult=[1.0, 1.0], xlim=[-5.0, 20.0]):
    xmin, xmax = xlim[0], xlim[1]
    ymin = np.min(data[(grid < xmax) * (grid >= xmin)])
    ymax = np.max(data[(grid < xmax) * (grid >= xmin)])
    dy = 0.05 * np.abs(ymax - ymin)
    if ymin != ymax:

        ymin_out = ymin - dy * y_mult[0]
        ymax_out = ymax + dy * y_mult[1]
        return ymin_out, ymax_out

    else:
        return


def get_fig_3plots(**kwargs):
    fig = newfig_generic_2yscale(
        1.4, scale_width=1.2,
        scale_ratio=0.5)    #(1.1,scale_width=1.5,scale_ratio=0.5)#plt.figure()
    fig.subplots_adjust(left=0.1, right=0.9, wspace=0.6, top=0.9, bottom=0.28)
    return fig


def get_fig_2plots(**kwargs):
    fig = newfig_generic(1.4, 0.9)
    plt.subplots_adjust(left=0.15, right=0.9, bottom=0.18, top=0.9, wspace=0.1)
    return fig
