import system
from matplotlib import pyplot
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import time

import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()

def deltA_from_deltU(alpha, deltU):
    return -deltU*alpha**2/(1.+deltU*alpha)

def val_drift_mono(cached =False):
    eigen_drifting = []
    alphas = []
    u_min= .01
    deltU = -.01
    magnets.load_state(1.15, down=True)
    u = 1/magnets.alpha

    log_plot = True
    running = True
    last = False

    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)
    figure.set_figwidth(6)
    subplots.set_xlabel('$1/\\alpha$', fontsize=16)
    subplots.set_ylabel('$\\omega^2/\\alpha$', fontsize=16)
    subplots.set_xlim(u_min, 1./magnets.alpha )

    if not log_plot:
        subplots.set_ylim(0,4)
    if cached:
        alphas, e_vecs, eigen_drifting = cache_poly_load(source='mono')
    # while 1./magnets.alpha > u_min and not cached:
    while running and not cached:
        print(magnets.alpha, 1./magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra_mono()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        # print(w2)
        e_vals = [tidy_eigs[i+1][0] for i in range(7)]
        eigen_drifting.append(e_vals)
        if u < u_min and last:
            break
        magnets.shift_alpha_and_stablize(deltA_from_deltU(magnets.alpha, deltU))
        u = 1/magnets.alpha
        if u < u_min:
            last=True
    u_vals = 1./numpy.array(alphas)
    for i in range(7):
        curve = [vals[i] for vals in eigen_drifting]
        if log_plot:
            scaled = curve*(u_vals)
            subplots.plot(u_vals, numpy.log(scaled))
        else:
            # subplots.plot(alphas, curve)
            scaled = curve*(u_vals)
            subplots.plot(u_vals, scaled)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs_mono_vs_u.png')

if __name__=='__main__':
    val_drift_mono()
