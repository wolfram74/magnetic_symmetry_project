import system
from matplotlib import pyplot
import time

import numpy
pi = numpy.pi

magnets = system.System()

vect_temp = '%.03f,'*7

def vectors_strip(array):
    #array of w^2 vector pairs
    #want 7 arrays of w^2 coord pairs
    output = [[] for i in range(8)]
    for entry in array:
        output[0].append(entry[0])

        for i in range(7):
            output[i+1].append(entry[i])
    return output

def color_vect(j, max_val=6):
    norm = float(j)/max_val
    # r = 1-norm**2
    r = 1-norm
    g = -(norm)*(norm-1)*4
    # b = norm**2
    b = norm
    return (r,g,b,1)


def eigen_vec_drift_plot():
    eigen_drifting = [[] for i in range(7)]
    alphas = []
    figure, subplots = pyplot.subplots(7)
    figure.set_figheight(28)
    figure.set_figwidth(6)
    magnets.load_state(.01)
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    a_max = 2.01
    # order preserved until 1.8 I think
    while magnets.alpha <a_max:
        print(magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        # print(w2)
        for i in range(7):
            eigen_drifting[i].append(tidy_eigs[i+1][1])
            #eigen drifting has vec for mode i (1,7)
        magnets.load_state(magnets.alpha+.011)
    #element in eigen drifting is array
        #conisting of w^2, vector pairs

    for i in range(7):
        subplots[i].set(title='$\\omega_%d$' % (i+1))
        subplots[i].set_ylim(-1,1)
        subplots[i].set_xlim(0,a_max)
        subplots[i].axhline(y=0,color='k',linestyle='--')
        subplots[i].grid(which='both', axis='x')

        lines = vectors_strip(eigen_drifting[i])
        for j in range(7):
            line = lines[j]
            subplots[i].plot(
                alphas, lines[j+1],
                color=colors[j],
                linestyle=styles[j],
                marker=markers[j],
                markevery=5
                )
    time_label = ("%d" % time.time())[-4:]
    pyplot.savefig(time_label+'-vectors.png')


if __name__ == '__main__':
    eigen_vec_drift_plot()
