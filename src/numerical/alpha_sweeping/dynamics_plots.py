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

def oscilation_path(alpha):
    magnets.load_state(alpha)
    del_state = shift_outers(magnets.state, 10.**-3)
    state_0 = numpy.array(magnets.state)
    magnets.state+=del_state
    path = [(0, del_state[:7])]
    E_0 = magnets.total_PE()+magnets.total_KE()
    print(magnets.gamma)
    while magnets.elapsed < 30:
        magnets.advance_in_time()
        E_t = magnets.total_PE()+magnets.total_KE()
        print(abs((E_t-E_0)/E_0))
        delta = magnets.state[:7]-state_0[:7]
        path.append((magnets.elapsed, tuple(delta)))
    return path

def plot_path(path, alpha):
    times = [el[0] for el in path]
    for i in range(7):
        del_phis = [el[1][i] for el in path]
        pyplot.plot(times, del_phis)
    pyplot.savefig('%d-alph_%s.png' % (time.time(), alpha))

if __name__=='__main__':
    alph = .5
    path = oscilation_path(alph)
    plot_path(path, alph)
