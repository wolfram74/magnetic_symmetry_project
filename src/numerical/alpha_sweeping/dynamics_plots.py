import system
from matplotlib import pyplot
import time

import numpy
pi = numpy.pi

magnets = system.System()


def shift_outers(state_in, epsi=10.**-2):
    delta = numpy.zeros(len(state_in))
    delta[1:7] = epsi
    return delta

def oscilation_path(alpha, del_state, T_end=5., down=False):
    magnets.load_state(alpha, down=down)
    magnets.elapsed = 0.
    state_0 = numpy.array(magnets.state)
    magnets.state+=del_state
    path = [(0, del_state[:7])]
    E_0 = magnets.total_PE()+magnets.total_KE()
    # print(magnets.gamma)
    # print(magnets.elapsed, T_end)
    while magnets.elapsed < T_end:
        magnets.advance_in_time()
        E_t = magnets.total_PE()+magnets.total_KE()
        delta = magnets.state[:7]-state_0[:7]
        path.append((magnets.elapsed, tuple(delta)))
    print(alpha, abs((E_t-E_0)/E_0))
    return path

def plot_path(path, alpha, title=None, ax=pyplot):
    times = [el[0] for el in path]
    for i in range(7):
        del_phis = [el[1][i] for el in path]
        ax.plot(times, del_phis)
    if not ax==pyplot:
        return
    if title:
        ax.savefig('%s.png' % title)
    else:
        ax.savefig('%d-alph_%s.png' % (time.time(), alpha))

def small_amp_alph_sweep():
    alph = .88
    del_state = shift_outers(magnets.state, 10.**-3)
    for i in range(1):
        alph += .01
        path = oscilation_path(alph, del_state)
        plot_path(path, alph)
        pyplot.clf()

def large_amp_alpha_sweep():
    alphs= [0., .5, 1., 2.4, 2.5, 3.]
    del_state = shift_outers(magnets.state, .2)
    for alph in alphs:
        path = oscilation_path(alph, del_state)
        plot_path(path, alph)
        pyplot.clf()

def eig_plots():
    alph = 1.0
    magnets.load_state(alph)
    # magnets.calc_L_mat(.01)
    # eig_vecs = numpy.linalg.eig(magnets.L_mat)[1]
    tidy_eigs = magnets.spectrum_finder()
    # for ind in range(len(eig_vecs)):
    template = "a=%.3f_w=%.3f.png"
    for eig in tidy_eigs:
        delt = numpy.zeros(2*7)
        vec = .01*eig[1]
        delt[:7] = vec
        freq = eig[0]**.5
        print(eig)
        print(delt[:7])
        # print(delt[7:])
        # continue
        path = oscilation_path(alph, delt, 4*pi/freq)
        title = template % (alph,freq)
        plot_path(path, alph, title)
        pyplot.clf()
        time.sleep(1)

def amplitdue_comparisons(alph=1., down=False):
    # alph = 1.0
    magnets.load_state(alph, down=down)
    tidy_eigs = magnets.spectrum_finder()
    tidy_eigs = sorted(tidy_eigs, key=lambda el : el[0])
    for eig in tidy_eigs:
        print(eig[0])
    amps = [10**-3, 10**-2, 10**-1]
    figure, subplots = pyplot.subplots(7,3)
    figure.set_figheight(28)
    figure.set_figwidth(18)
    faxe_path = [[0, tuple(range(7))],[1, tuple(range(7))]]
    for j in range(3):
        subplots[0][j].set(title='$|\\delta|=%03f$' % amps[j])
    for i in range(7):
        # subplots[i][0].set(yaxis='$\\omega=%f$'%tidy_eigs[i][0])
        subplots[i][0].set_ylabel('$\\omega^2=%f$'%tidy_eigs[i][0])
    for j in range(3):
        amp = amps[j]
        print(amp)
        y_rang = numpy.linspace(-.25*amp, .25*amp, 20)
        x_rang = numpy.ones(len(y_rang))
        for i in range(7):
            w2 = tidy_eigs[i][0]
            vec = tidy_eigs[i][1]
            w = w2**.5
            T = 2*pi/w
            delt = numpy.zeros(2*7)
            delt[:7] = [amp*el for el in vec]
            path = oscilation_path(alph, delt, 2*T, down=down)
            plot_path(path, alph, ax=subplots[i][j])
            subplots[i][j].plot(x_rang*T, y_rang, linestyle='--')
    time_label = ("%d" % time.time())[-4:]
    pyplot.savefig(time_label+'_%.03f_-comparative.png' % alph)



if __name__=='__main__':
    # large_amp_alpha_sweep()
    # small_amp_alph_sweep()
    # eig_plots()
    # alphs = [.5, 1., 1.5, 2.4, 3.]
    alphs = [1.17, 2.47]
    for alph in alphs:
        amplitdue_comparisons(alph,down=True)
