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
    u_min= 10**-3
    deltU = -.01
    magnets.load_state(1.15, down=True)
    u = 1/magnets.alpha

    log_plot = False
    running = True
    last = False
    OMEGAS = [0.1815, 1.3229, 1.3229, 1.3229, 1.9126, 2.0000, 10.5203]
    figure, subplots = pyplot.subplots(1)
    figure.set_figheight(6)
    figure.set_figwidth(6)
    # subplots.set_xlabel('$\log(1/\\alpha)$', fontsize=16)
    # subplots.set_xlim(numpy.log10(u_min), numpy.log10(1/magnets.alpha))
    subplots.set_xlabel('$1/\\alpha$', fontsize=16)
    subplots.set_xlim(0,1/magnets.alpha)
    # subplots.set_ylabel('$\\frac{\\omega^2/\\alpha}{\Omega_i^2}$', fontsize=16)
    subplots.set_ylabel('$R_i$', fontsize=16)

    if not log_plot:
        subplots.set_ylim(0,2.5)
    if cached:
        alphas, e_vecs, eigen_drifting = cache_poly_load(source='mono')
    # while 1./magnets.alpha > u_min and not cached:
    while running and not cached:
        print(magnets.alpha, 1./magnets.alpha)
        alphas.append(magnets.alpha)
        tidy_eigs = magnets.labeled_spectra_mono()
        #tidy_eigs is a dict with w^2, vec pairs
        # w2 = [val[0] for val in tidy_eigs]
        print(len(tidy_eigs.keys()))
        e_vals = [tidy_eigs[i+1][0] for i in range(7)]
        eigen_drifting.append(e_vals)
        if u < u_min and last:
            break
        if magnets.alpha < 5.:
            delta_alph = deltA_from_deltU(magnets.alpha, deltU)
        if magnets.alpha < 5. and magnets.alpha+delta_alph < 2.5:
            magnets.load_state(magnets.alpha+delta_alph, down=True)
        else:
            magnets.load_state(magnets.alpha*1.151)
        # magnets.shift_alpha_and_stablize(deltA_from_deltU(magnets.alpha, deltU))
        u = 1/magnets.alpha
        if u < u_min:
            last=True
    u_vals = 1./numpy.array(alphas)
    for i in range(7):
        curve = [vals[i] for vals in eigen_drifting]
        if log_plot:
            scaled = curve*(u_vals)
            subplots.plot(u_vals, numpy.log10(scaled))
        else:
            # subplots.plot(alphas, curve)
            scaled = curve*(u_vals)/OMEGAS[i]
            # scaled = curve*(u_vals)
            # scaled = numpy.sqrt(curve*(u_vals))
            subplots.plot(u_vals, scaled)
            # subplots.plot(numpy.log10(u_vals), scaled)
    # val_mono_annotate(subplots)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-freqs_mono_vs_u.png', bbox_inches='tight')


def val_mono_annotate(subplot):
    subplot.annotate(s='$\\omega_1$', xy=(.05,.25))
    subplot.annotate(s='$\\omega_2$', xy=(.1,.85))
    subplot.annotate(s='$\\omega_3$', xy=(.1,1.176))
    subplot.annotate(s='$\\omega_4$', xy=(.05,1.55))
    subplot.annotate(s='$\\omega_5$', xy=(.05,1.725))
    subplot.annotate(s='$\\omega_6$', xy=(.05,1.95))
    subplot.annotate(s='$\\omega_7 \Uparrow$', xy=(.10,3.5))


if __name__=='__main__':
    val_drift_mono()
