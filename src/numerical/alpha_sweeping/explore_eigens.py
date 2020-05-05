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

def cache_poly_load(source='poly'):
    template = './saved_eigenmodes/'+source+'_%d_mode.txt'
    alphas = []
    e_vals = []
    eigen_drifting = [[] for i in range(7)]
    eigen_vals = [[] for i in range(7)]
    for i in range(7):
        address = template % (i+1)
        cache = open(address, 'r')
        for line in cache:
            split = line.rstrip().split(',')
            alpha = float(split[0])
            val = float(split[1])
            # print(split[2:])
            vector = [float(el) for el in split[2:-1]]
            if i==0:
                alphas.append(alpha)
            eigen_vals[i].append(val)
            eigen_drifting[i].append(vector)
        cache.close()
    for ind in range(len(eigen_vals[0])):
        e_vals.append([mode[ind] for mode in eigen_vals])
    return alphas, eigen_drifting, e_vals

def eigen_vec_drift_plot(cached = True):
    eigen_drifting = [[] for i in range(7)]
    alphas = []
    figure, subplots = pyplot.subplots(7)
    figure.set_figheight(28)
    figure.set_figwidth(6)
    magnets.load_state(.01)
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    marker_size = []
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    a_max = 2.48
    if cached:
        alphas, eigen_drifting, freqs= cache_poly_load()
    while magnets.alpha <a_max:
        if cached:
            break
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
    sub_label = 'a b c d e f g'.split(' ')
    for i in range(7):
        subplots[i].set_title(label='(%s) $\\omega_%d$' % (sub_label[i],(i+1)), loc='right')
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
                markevery=3, markersize=3
                )
    vec_poly_annotate(subplots)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-vectors.png')

def vec_poly_annotate(subplots):
    subplots[0].annotate(s='$\\phi_0$', xy=(.05,.9))
    subplots[0].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,.25))
    subplots[0].annotate(s='$\\phi_1, \\phi_4$', xy=(1.0,0.1))
    subplots[0].annotate(s='$\\phi_2, \\phi_3$', xy=(.5,-.3))
    # subplots[0].annotate(s='a', xy=(2.5,.9))

    subplots[1].annotate(s='$\\phi_0$', xy=(.1,.05))
    subplots[1].annotate(s='$\\phi_5$', xy=(.5,.6))
    subplots[1].annotate(s='$\\phi_4$', xy=(1.0,.40))
    subplots[1].annotate(s='$\\phi_2$', xy=(.5,.2))
    subplots[1].annotate(s='$\\phi_3$', xy=(.5,-.2))
    subplots[1].annotate(s='$\\phi_1$', xy=(1.0,-.40))
    subplots[1].annotate(s='$\\phi_6$', xy=(.5,-.65))

    subplots[2].annotate(s='$\\phi_0$', xy=(.9,.05))
    subplots[2].annotate(s='$\\phi_2$', xy=(.25,.60))
    subplots[2].annotate(s='$\\phi_6$', xy=(.25,.35))
    subplots[2].annotate(s='$\\phi_4$', xy=(.25,.15))
    subplots[2].annotate(s='$\\phi_1$', xy=(.25,-.20))
    subplots[2].annotate(s='$\\phi_5$', xy=(.25,-.40))
    subplots[2].annotate(s='$\\phi_3$', xy=(.25,-.7))

    subplots[3].annotate(s='$\\phi_0$', xy=(.05,.075))
    subplots[3].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,.4))
    subplots[3].annotate(s='$\\phi_2, \\phi_3$', xy=(.25,.15))
    subplots[3].annotate(s='$\\phi_1, \\phi_4$', xy=(.25,-0.5))

    subplots[4].annotate(s='$\\phi_0$', xy=(.1,.1))
    subplots[4].annotate(s='$\\phi_1$', xy=(.25,.65))
    subplots[4].annotate(s='$\\phi_2$', xy=(.25,.40))
    subplots[4].annotate(s='$\\phi_5$', xy=(.25,.10))
    subplots[4].annotate(s='$\\phi_6$', xy=(.25,-.15))
    subplots[4].annotate(s='$\\phi_3$', xy=(.25,-.45))
    subplots[4].annotate(s='$\\phi_4$', xy=(.25,-.75))

    subplots[5].annotate(s='$\\phi_0$', xy=(.8,.5))
    subplots[5].annotate(s='$\\phi_2, \\phi_3$', xy=(.25,.5))
    subplots[5].annotate(s='$\\phi_1, \\phi_4$', xy=(.25,-0.22))
    subplots[5].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,-.65))

    subplots[6].annotate(s='$\\phi_0$', xy=(.05,.075))
    subplots[6].annotate(s='$\\phi_2, \\phi_3$', xy=(.75,.65))
    subplots[6].annotate(s='$\\phi_1, \\phi_4$', xy=(1.2,0.33))
    subplots[6].annotate(s='$\\phi_5,\\phi_6$', xy=(1.2,-.3))

def eigen_val_drift_plot(cached =False):
    eigen_drifting = []
    alphas = []
    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)

    figure.set_figwidth(6)
    subplots.set_ylim(0,4)
    subplots.set_xlabel('$\\alpha$', fontsize=16)
    subplots.set_ylabel('$\\omega$', fontsize=16)
    magnets.load_state(.01)
    a_max=2.48
    if cached:
        alphas, e_vecs, eigen_drifting = cache_poly_load()
    while magnets.alpha <a_max and not cached:
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
        subplots.plot(alphas, numpy.sqrt(curve))
    val_poly_annotate(subplots)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs.png')

def val_poly_annotate(subplot):
    subplot.annotate(s='$\\omega_1$', xy=(1.0,.4))
    subplot.annotate(s='$\\omega_2$', xy=(1.0,1.52))
    subplot.annotate(s='$\\omega_3$', xy=(1.0,1.72))
    subplot.annotate(s='$\\omega_4$', xy=(1.0,1.87))
    subplot.annotate(s='$\\omega_5$', xy=(1.0,2.2))
    subplot.annotate(s='$\\omega_6$', xy=(1.0,2.6))
    subplot.annotate(s='$\\omega_7$', xy=(1.0,2.92))

def vec_drift_mono(cached = False):
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
    if cached:
        alphas, eigen_drifting, e_vals = cache_poly_load(source='mono')
    while magnets.alpha > a_min and not cached:
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
                markevery=5, markersize = 3
                )
    vec_mono_annotate(subplots)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-vectors_mono.png')

def vec_mono_annotate(subplots):
    subplots[0].annotate(s='$\\phi_0$', xy=(1.5,.5))
    subplots[0].annotate(s='$\\phi_1, \\phi_4$', xy=(1.5,0.1))
    subplots[0].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(1.4,-.6))

    subplots[2].annotate(s='$\\phi_0, \\phi_1, \\phi_4$', xy=(1.5,.1))
    subplots[2].annotate(s='$\\phi_2, \\phi_6$', xy=(1.5,0.6))
    subplots[2].annotate(s='$\\phi_3, \\phi_5$', xy=(1.5,-.63))

    subplots[4].annotate(s='$\\phi_0$', xy=(1.75,.15))
    subplots[4].annotate(s='$\\phi_1, \\phi_4$', xy=(1.5,0.8))
    subplots[4].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(1.75,-.4))

    subplots[5].annotate(s='$\\phi_0, \\phi_1, \\phi_4$', xy=(1.5,.1))
    subplots[5].annotate(s='$\\phi_2, \\phi_3$', xy=(1.5,0.6))
    subplots[5].annotate(s='$\\phi_5, \\phi_6$', xy=(1.5,-.63))

    subplots[6].annotate(s='$\\phi_0$', xy=(1.5,.79))
    subplots[6].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(1.5,.3))
    subplots[6].annotate(s='$\\phi_1, \\phi_4$', xy=(1.5,-0.1))
    # subplots[0].annotate(s='$\\phi_2, \\phi_3$', xy=(1.6,-.6))

    subplots[1].annotate(s='$\\phi_0$', xy=(1.2,.05))
    subplots[1].annotate(s='$\\phi_2,\\phi_5$', xy=(1.5,.6))
    subplots[1].annotate(s='$\\phi_4$', xy=(1.5,.25))
    subplots[1].annotate(s='$\\phi_1$', xy=(1.5,-.30))
    subplots[1].annotate(s='$\\phi_3,\\phi_6$', xy=(1.5,-.63))

    subplots[3].annotate(s='$\\phi_0$', xy=(1.2,.05))
    subplots[3].annotate(s='$\\phi_4$', xy=(1.75,.75))
    subplots[3].annotate(s='$\\phi_3,\\phi_6$', xy=(1.8,.3))
    subplots[3].annotate(s='$\\phi_2,\\phi_5$', xy=(1.8,-.35))
    subplots[3].annotate(s='$\\phi_1$', xy=(1.75,-.60))

def val_drift_mono(cached =False):
    eigen_drifting = []
    alphas = []
    a_min=1.15
    magnets.load_state(3.5)
    log_plot = False
    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)
    figure.set_figwidth(6)
    subplots.set_xlabel('$\\alpha$', fontsize=16)
    subplots.set_ylabel('$\\omega$', fontsize=16)
    subplots.set_xlim(a_min,3.5)
    if not log_plot:
        subplots.set_ylim(0,4)
    if cached:
        alphas, e_vecs, eigen_drifting = cache_poly_load(source='mono')
    while magnets.alpha >a_min and not cached:
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
        if log_plot:
            subplots.plot(alphas, numpy.log(curve))
        else:
            # subplots.plot(alphas, curve)
            subplots.plot(alphas, numpy.sqrt(curve))
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs_mono.png')

def mono_mode1():
    magnets = system.System()
    alph_max = 2.1
    alph_min = 2.0
    del_alpha = -.001
    magnets.load_state(alph_max, down=True)
    alphas = []
    omeg2s = []

    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)
    figure.set_figwidth(6)
    subplots.set_xlabel('$\\alpha$', fontsize=16)
    subplots.set_ylabel('\t\t$\\omega^2$', fontsize=16, labelpad=-10)
    subplots.set_ylim(0.000,.002)

    # print(magnets.alpha, alph_min)
    # print(magnets.alpha >= alph_min)

    while magnets.alpha >= alph_min:
        alphas.append(magnets.alpha)
        eigs = magnets.labeled_spectra_mono()
        omeg2s.append(eigs[1][0])
        # print('%.9f\t%.9f' %(alphas[-1], omeg2s[-1]))
        magnets.shift_alpha_and_stablize(del_alpha)

    # subplots.plot(alphas, numpy.sqrt(omeg2s))
    subplots.plot(alphas, omeg2s)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-mono_mode1sqrt.png')



def eig_vec_schematic():
    figure, subplots = pyplot.subplots(7)
    figure.set_figheight(28)
    figure.set_figwidth(6)
    magnets.load_state(2.05, down=True)
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
    # eigen_vec_drift_plot()
    # eigen_val_drift_plot(cached=True)
    # vec_drift_mono(cached=True)
    # val_drift_mono(cached=True)
    # eig_vec_schematic()
    mono_mode1()
