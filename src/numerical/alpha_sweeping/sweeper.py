import system
from matplotlib import pyplot
import time

from numpy import pi


def gen_text():
    magnets = system.System()
    equilibria = gen_states(magnets)
    store = open('stored_states.txt', 'a')
    template = '%.9f, '*8+'\n'
    print(template)
    for state in equilibria:
        store.write(template % state)
    store.close()

def gen_down_text():
    magnets = system.System()
    equilibria = goin_down_states(magnets)
    store = open('stored_states_down.txt', 'w')
    template = '%.9f, '*8+'\n'
    print(template)
    for state in equilibria:
        store.write(template % state)
    store.close()

def gen_states(magnets):
    running = True
    monitor = 10.
    magnets.load_state(1243)
    forward = True
    equilibria = []
    first_step = True
    while running:
        delta_alpha = magnets.alpha*0.15
        print(delta_alpha)
        magnets.shift_alpha_and_stablize(delta_alpha)
        entry = [magnets.alpha]+list(magnets.state[:7])
        equilibria.append(tuple(entry))
        if magnets.alpha >= 1.15**105:
            running = False
    return equilibria

def goin_down_states(magnets):
    running = True
    monitor = 10.
    magnets.load_state(2.5)
    forward = True
    equilibria = []
    first_step = True
    while running:
        magnets.shift_alpha_and_stablize(-0.01)
        entry = [magnets.alpha]+list(magnets.state[:7])
        equilibria = [tuple(entry)]+equilibria
        if magnets.alpha <= .8:
            running = False
    return equilibria

def eigen_modes_poly(poly=True):
    magnets = system.System()
    if poly:
        magnets.load_state(.01)
        address_template = './saved_eigenmodes/poly_%d_mode.txt'
        amax = 2.48
        looping = lambda alph: alph<amax
        shift = lambda mags: mags.load_state(mags.alpha+.011)
        mode_find = lambda mags: mags.labeled_spectra()

    if not poly:
        magnets.load_state(3.5)
        address_template = './saved_eigenmodes/mono_%d_mode.txt'
        amin=1.15
        looping = lambda alph: alph>amin
        shift = lambda mags: mags.load_state(mags.alpha-.009, down=mono_specific)
        mode_find = lambda mags: mags.labeled_spectra_mono()

    destinations = []
    eigen_states = mode_find(magnets)
    for key in eigen_states.keys():
        destinations.append(open(address_template%key, 'w'))


    #save data is alpha, w^2, phi vector
    line_template = '%.9f,'+'%.9f,'+'%.9f,'*7+'\n'
    while looping(magnets.alpha):
        eigen_states = mode_find(magnets)
        for key in eigen_states.keys():
            data = [magnets.alpha, eigen_states[key][0]]+list(eigen_states[key][1])
            destinations[key-1].write(line_template % tuple(data))
        mono_specific = magnets.alpha < 2.5 and not poly
        shift(magnets)


if __name__ =='__main__':
    gen_text()
    # gen_down_text()
    # eigen_modes_poly()
    # eigen_modes_poly(poly=False)


# for ind in range(7):
#     PE0 = magnets.total_PE()
#     magnets.state[ind]+=.01
#     print(magnets.total_PE()-PE0)
#     magnets.state[ind]-=.01
