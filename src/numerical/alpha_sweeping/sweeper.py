import system
from matplotlib import pyplot
import time

from numpy import pi


def gen_text():
    magnets = system.System()
    equilibria = gen_states(magnets)
    store = open('stored_states.txt', 'w')
    template = '%f, '*8+'\n'
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
        magnets.advance_in_time()
        if first_step:
            print(magnets.elapsed)
            first_step = False

        if magnets.elapsed > 2.:
            magnets.gamma = 2*(magnets.alpha+1)*(magnets.elapsed-2.)/magnets.elapsed

        if magnets.total_delta_sqr() > 10**-8 or magnets.elapsed<5.:
            if magnets.elapsed > monitor:
                monitor+=10
                print('monitor at', magnets.elapsed, magnets.alpha)
                print('delta sqr %f' % magnets.total_delta_sqr())
                print('\n')
            continue
        magnets.elapsed = 0.
        magnets.gamma = 0.
        print('state', magnets.state[:7])
        print('PE', magnets.total_PE())
        entry = [magnets.alpha]+list(magnets.state[:7])
        equilibria.append(tuple(entry))
        monitor = 10
        magnets.alpha+=.01
        first_step = True

        if magnets.alpha >= 4.:
            running = False
    return equilibria



if __name__ =='__main__':
    gen_text()


# for ind in range(7):
#     PE0 = magnets.total_PE()
#     magnets.state[ind]+=.01
#     print(magnets.total_PE()-PE0)
#     magnets.state[ind]-=.01
