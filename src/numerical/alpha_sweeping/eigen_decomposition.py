import system
import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()

import time
from matplotlib import pyplot
from matplotlib import backend_bases
from matplotlib import transforms
from matplotlib import rcParams
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


'''
f(x+h)-f(x-h)/2h
h = .01
usable alphas, .02 through 2.46
'''

def gen_deriv_equis():
    dphi_dalph = []
    magnets.load_state(.01)
    h2 = .02
    lower = magnets.state[:7]
    magnets.load_state(magnets.alpha+.01)
    # print(lower,magnets.alpha)
    while magnets.alpha < 2.47:
        alpha = magnets.alpha
        # print(alpha)
        current = magnets.state[:7]
        magnets.load_state(magnets.alpha+.011)
        higher = magnets.state[:7]
        dphi_dalph.append( (alpha, (higher-lower)/2))
        lower = current
    return dphi_dalph

def display_derivs(values):
    template = '%.4f||'*7
    print(len(values))
    for line in values:
        # print(line[0])
        amped_vals = [el*(10**4) for el in line[1]]
        # print(amped_vals)
        print(line[0], template%tuple(amped_vals))
            # print(line[0], line[1])

def change_basis(values):
    '''
    value[i][1] = sum(c_j*eig_vec_j)
    output is the equilibrium expressed in the eigen vector basis, c_j's
    '''
    dphi_in_eigens = []
    for line in values:
        magnets.load_state(line[0])
        eigs = magnets.labeled_spectra()
        factors = []
        for label in eigs.keys():
            vec_i = eigs[label][1]
            factors.append(numpy.dot(vec_i, line[1]))
            # print(numpy.dot(vec_i, vec_i))
        dphi_in_eigens.append((line[0], factors))
    return dphi_in_eigens

def plot_curves(values):
    curves = [[] for el in range(7)]
    alphas = numpy.linspace(.02, 2.46, num=len(values))
    # alphas = numpy.linspace(.02, 2.46, num=245)
    print(alphas)
    for line in values:
        for ind in range(7):
            curves[ind].append(numpy.log(line[1][ind]**2))
    for curve in curves:
        pyplot.plot(alphas, curve)

    pyplot.ylim(-30,0.)
    pyplot.xlim(0,2.5)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+('decompositio.png'))


if __name__ == '__main__':
    derivatives = gen_deriv_equis()
    # print('derivs done')
    derivs_star = change_basis(derivatives)
    print('new basis')
    # display_derivs(derivs_star)
    plot_curves(derivs_star)