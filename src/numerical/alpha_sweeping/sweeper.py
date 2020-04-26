import system
from matplotlib import pyplot
import time

from numpy import pi


def gen_text():
    magnets = system.System()
    equilibria = gen_states(magnets)
    store = open('stored_states.txt', 'w')
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
    magnets.alpha = 0.
    forward = True
    equilibria = []
    first_step = True
    while running:
        magnets.shift_alpha_and_stablize(0.01)
        entry = [magnets.alpha]+list(magnets.state[:7])
        equilibria.append(tuple(entry))
        if magnets.alpha >= 4.:
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


if __name__ =='__main__':
    # gen_text()
    gen_down_text()


# for ind in range(7):
#     PE0 = magnets.total_PE()
#     magnets.state[ind]+=.01
#     print(magnets.total_PE()-PE0)
#     magnets.state[ind]-=.01
