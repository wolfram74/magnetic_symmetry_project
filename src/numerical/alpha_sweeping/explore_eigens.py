import system
from matplotlib import pyplot
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import time

import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

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
    a_max = 2.48
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
        subplots[i].set_title(label='$\\omega_%d$' % (i+1), loc='right')
        subplots[i].set_ylim(-1,1)
        subplots[i].set_xlim(0,a_max)
        subplots[i].set_xlabel('$\\alpha$', fontsize=16)
        subplots[i].set_ylabel('$\\delta \\phi$', fontsize=16)
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
                markevery=3
                )
    subplots[0].annotate(s='$\\phi_0$', xy=(.05,.9))
    subplots[0].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,.25))
    subplots[0].annotate(s='$\\phi_1, \\phi_4$', xy=(1.0,0.1))
    subplots[0].annotate(s='$\\phi_2, \\phi_3$', xy=(.5,-.3))
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-vectors.png')

def eigen_val_drift_plot():
    eigen_drifting = []
    alphas = []
    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)

    figure.set_figwidth(6)
    subplots.set_ylim(0,10)
    subplots.set_xlabel('$\\alpha$', fontsize=16)
    subplots.set_ylabel('$\\omega^2$', fontsize=16)
    magnets.load_state(.01)
    a_max=2.48
    while magnets.alpha <a_max:
        print(magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        # print(w2)
        e_vals = [tidy_eigs[i+1][0] for i in range(7)]
        eigen_drifting.append(e_vals)
        magnets.load_state(magnets.alpha+.011)
    for i in range(7):
        curve = [vals[i] for vals in eigen_drifting]
        # subplots.plot(alphas, numpy.log(curve))
        subplots.plot(alphas, curve)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs.png')

def vec_drift_mono():
    eigen_drifting = [[] for i in range(7)]
    alphas = []
    figure, subplots = pyplot.subplots(7)
    figure.set_figheight(28)
    figure.set_figwidth(6)
    magnets.load_state(3.5)
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    a_min = 1.15
    # a_min = 1.15 # goal
    # order preserved until 1.8 I think
    while magnets.alpha > a_min:
        print(magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra_mono()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        # print(w2)
        for i in range(7):
            try:
                eigen_drifting[i].append(tidy_eigs[i+1][1])
            except:
                print(alphas[-1])
                print(i)
                print(tidy_eigs)
                return
            #eigen drifting has vec for mode i (1,7)
        mono_specific= magnets.alpha<2.5
        magnets.load_state(magnets.alpha-.009, down=mono_specific)
    #element in eigen drifting is array
        #conisting of w^2, vector pairs

    for i in range(7):
        subplots[i].set_title(label='$\\omega_%d$' % (i+1), loc='right')
        subplots[i].set_ylim(-1,1)
        subplots[i].set_xlim(a_min,3.5)
        subplots[i].set_xlabel('$\\alpha$', fontsize=16)
        subplots[i].set_ylabel('$\\delta \\phi$', fontsize=16)
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
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-vectors_mono.png')

def val_drift_mono():
    eigen_drifting = []
    alphas = []
    a_min=1.15
    magnets.load_state(3.5)

    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)
    figure.set_figwidth(6)
    subplots.set_xlabel('$\\alpha$', fontsize=16)
    subplots.set_ylabel('$\\omega^2$', fontsize=16)
    subplots.set_xlim(a_min,3.5)
    subplots.set_ylim(0,10)
    while magnets.alpha >a_min:
        print(magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra_mono()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        # print(w2)
        e_vals = [tidy_eigs[i+1][0] for i in range(7)]
        eigen_drifting.append(e_vals)
        mono_specific= magnets.alpha<2.5
        magnets.load_state(magnets.alpha-.009, down=mono_specific)
    for i in range(7):
        curve = [vals[i] for vals in eigen_drifting]
        subplots.plot(alphas, curve)
        # subplots.plot(alphas, numpy.log(curve))
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs_mono.png')

def eig_vec_schematic():
    figure, subplots = pyplot.subplots(7)
    figure.set_figheight(28)
    figure.set_figwidth(6)
    magnets.load_state(.75)
    named_eigs = magnets.labeled_spectra()
    equilib = magnets.state[:7]
    for index in range(7):
        mode_ID = index+1
        mode = named_eigs[mode_ID]
        single_mode(
            subplots[index], equilib, mode,
            magnets.alpha, mode_ID
            )
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-eig_schematic.png')
    return

def single_mode(frame, gam_0, mode, alpha, mode_ID):
    arrow_patches = []
    colors = []
    title = (mode_ID, alpha)
    frame.set_ylim(-1.75,1.75)
    frame.set_xlim(-1.75,1.75)
    for ind in range(7):
        phi0 = gam_0[ind]
        phi = mode[1][ind]
        tht0i = 2*pi*(ind-1)/6.
        if ind ==0:
            ci = [0.,0.]
        else:
            ci = [cos(tht0i), sin(tht0i)]
        arrow_patches.append(
            mpatches.Rectangle(ci, .5, .005, 180*phi0/pi, fill=False)
            )
        arrow_patches.append(
            mpatches.Arrow(
                ci[0], ci[1], cos(phi0+phi)/2.2, sin(phi0+phi)/2.2,
                width=.1, zorder = 1.,
                # capstyle='round', fill=False,
                # alpha=.5, hatch='o'
                )
            )
        arrow_patches.append(
            mpatches.Circle(
                ci, radius=.5, alpha=.1, zorder = .5
                )
            )
        colors.append((.1,.5,.1, .5))
        colors.append((.1,.1,.5))
        colors.append((.5,.1,.1,.1))

    shapes = PatchCollection(
        arrow_patches,
        facecolors=colors
        )
    frame.set_aspect('equal')
    frame.add_collection(shapes)
    frame.set_title('$\\omega_%d$, $\\alpha=%.3f$' % title)
    return

if __name__ == '__main__':
    eigen_vec_drift_plot()
    # eigen_val_drift_plot()
    # vec_drift_mono()
    # val_drift_mono()
    # eig_vec_schematic()
