import system
from matplotlib import pyplot
from matplotlib import backend_bases
from matplotlib import transforms
from matplotlib import rcParams
import time
import eight_by_two_plot as other_plots

import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

def gen_fig_and_subplots():
    L=6
    figure, axes = pyplot.subplots()
    # figure.set_figheight(L)
    figure.set_figwidth(L)
    # axes.set_aspect('equal')
    # axes.set_aspect(1.0/axes.get_data_ratio(), adjustable='box')
    return figure, axes

def plot_eigen_curves(poly=True):
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    sub_label = 'a b c d a b c d'.split(' ')
    pyplot.rc('xtick', labelsize=14)
    pyplot.rc('ytick', labelsize=14)
    # figure, subplots = pyplot.subplots(4,2)
    # figure.set_figheight(24)
    # figure.set_figwidth(12)
    # tick marks
    #https://matplotlib.org/3.1.1/gallery/ticks_and_spines/tick-locators.html#sphx-glr-gallery-ticks-and-spines-tick-locators-py
    square_layout = False
    figure, subplots = gen_fig_and_subplots()


    if poly:
        print('poly mode')
        prefix = 'poly'
        x_axis_label = '$\\alpha$'
        y_axis_off = -.08 #how shifted is the y axis label
        y_axis_0 = 0.03
        x_axis_off =  -0.035 #how shifted is the x axis label
        eig_val_yaxis = '$\\omega_i^2 \\Omega^{-2}$'
        eig_val_title = '$\\omega^2 \\Omega^{-2}$ vs. $\\alpha$'
        get_values = other_plots.load_poly_values
        eig_limit = 10
        limit_sets = lambda x: (0, x[-1])
        annotate = poly_annotate
    else:
        print('mono mode')
        prefix = 'mono'
        x_axis_label = '$\\alpha^{-1}$'
        y_axis_off = -.13 #how shifted is the y axis label
        y_axis_0 = -0.00
        x_axis_off = -.075 #how shifted is the x axis label
        eig_val_yaxis = '$\\omega^2 \\alpha^{-1} \\Omega^{-2}$'
        eig_val_title = '$\\omega^2 \\alpha^{-1}$ vs. $\\alpha^{-1}$'
        eig_limit = 4
        limit_sets = lambda x: (x[-1], x[0])
        # limit_sets = lambda x: (1.0, 2.0)
        get_values = other_plots.load_mono_values
        annotate = mono_annotate

    x_values, eigen_vals, eigen_vecs = get_values()
    #x_values, 1-D list of alpha or 1/alpha
    #eigen values: same length as x_values, each entry 7 eigen values at that x-value
    # eigen_vecs: 7D list, one for each mode, x-value number of entries, of 7 vectors each

    x_min, x_max = limit_sets(x_values)
    yticks = (-1., -.5, 0, .5, 1.)
    fontsize=20
    subplots.set_xlabel(x_axis_label, fontsize=fontsize)
    subplots.set_ylabel(eig_val_yaxis, fontsize=fontsize)
    subplots.set_ylim(0, eig_limit)
    subplots.set_xlim(x_min, x_max)
    subplots.xaxis.set_label_coords(.5, x_axis_off)
    subplots.yaxis.set_label_coords(y_axis_off+y_axis_0, .5)
    subplots.tick_params(axis='both', labelsize=14)
    subplots.set_title(
        label = eig_val_title,
        loc='right', fontsize=fontsize)

    for i in range(7):
        y_curve = [entry[i] for entry in eigen_vals]
        subplots.plot(x_values, y_curve)
        # x_values = numpy.array(x_values)
        # y_curve = numpy.array(y_curve)
        # subplots.plot(1.0/x_values, y_curve/x_values)

    annotate(subplots)

    subplots.set_aspect(1.0/subplots.get_data_ratio(), adjustable='box')
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+('-%s-eigen_values.pdf' % prefix))



def poly_annotate(subplots):
    fsize = 16
    subplots.annotate(text='$1$', xy=(1.0,.6**2),fontsize=fsize)
    subplots.annotate(text='$2$', xy=(.75,1.52**2),fontsize=fsize)
    subplots.annotate(text='$3$', xy=(1.25,1.75**2),fontsize=fsize)
    subplots.annotate(text='$4$', xy=(.95,1.87**2),fontsize=fsize)
    subplots.annotate(text='$5$', xy=(1.0,2.2**2),fontsize=fsize)
    subplots.annotate(text='$6$', xy=(1.0,2.55**2),fontsize=fsize)
    subplots.annotate(text='$7$', xy=(.90,2.92**2),fontsize=fsize)

def mono_annotate(subplots):
    fsize = 16
    subplots.annotate(text='$1$', xy=(.05,.26),fontsize=fsize)
    subplots.annotate(text='$2$', xy=(.15,.57),fontsize=fsize)
    subplots.annotate(text='$3$', xy=(.15,1.156),fontsize=fsize)
    subplots.annotate(text='$4$', xy=(.275,2.1),fontsize=fsize)
    subplots.annotate(text='$5$', xy=(.5,1.67),fontsize=fsize)
    subplots.annotate(text='$6$', xy=(.025,2.1),fontsize=fsize)
    subplots.annotate(text='$7 \\Uparrow$', xy=(.10,3.5),fontsize=fsize)

if __name__ == '__main__':
    # plot_eigen_curves()
    plot_eigen_curves(poly=False)
    