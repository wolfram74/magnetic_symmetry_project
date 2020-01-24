import system
from matplotlib import pyplot
import time

from numpy import pi

def gen_text(start_state):
    magnets = system.System()
    equilibria = []
    running = True
    # running = False
    monitor = 10.
    magnets.state *= 0.
    magnets.alpha = 0.
    magnets.state[:7] = start_state
    print(magnets.state)
    while running:
        magnets.advance_in_time()
        if magnets.elapsed >5.:
            magnets.gamma = 2*(magnets.alpha+1)*(magnets.elapsed-5.)/magnets.elapsed
        if magnets.total_delta_sqr()>1*10**-8 or magnets.elapsed<10.:
            if magnets.elapsed > monitor:
                monitor+=10
                print('monitor at', magnets.elapsed)
                print(magnets.gamma, magnets.step_size, magnets.total_delta_sqr())
                # print(magnets.state[:7])
                # print(magnets.state[7:])
                print('\n')
            continue
        print('done at %f' % magnets.elapsed)
        print(magnets.total_delta_sqr(), magnets.total_PE())
        print(magnets.state[:7])
        print('\n')
        equilibria.append(
            [magnets.alpha,
            tuple(magnets.state[:7]),
            magnets.net_dipole_mag(), magnets.total_PE()]
            )
        # if magnets.alpha < 1.:
        magnets.alpha+=.05
        # else:
        #         magnets.alpha *= 1.075
        print(magnets.alpha)
        # print(magnets.total_force())
        magnets.elapsed = 0.
        magnets.gamma = 0.
        monitor = 10
        if magnets.alpha > 3.01:
            running = False
    angle_plotter(equilibria)
    state_plotter(equilibria)

def angle_plotter(equilibs):
    gam_0 = equilibs[0][1]
    figure, subplots = pyplot.subplots(2,1)
    alpha = [el[0] for el in equilibs]
    for ind in range(7):
        phi_i = [el[1][ind]/pi for el in equilibs]
        del_phi_i = [(gam_0[ind] - el[1][ind])/pi for el in equilibs]
        subplots[0].plot(alpha, phi_i, label='%d' %ind)
        subplots[1].plot(alpha, del_phi_i, label='%d' %ind)
    pyplot.legend()
    pyplot.xlabel('$\\alpha$', fontsize=16)
    pyplot.savefig('%d-ang.png' % time.time())
    return

def state_plotter(equilibs):
    figure, subplots = pyplot.subplots(2,1)
    alpha = [el[0] for el in equilibs]
    mu = [el[2] for el in equilibs]
    PE = [el[3] for el in equilibs]
    subplots[0].plot(alpha, mu, label='mu')
    subplots[1].plot(alpha, PE, label='PE')
    pyplot.legend()
    pyplot.xlabel('$\\alpha$', fontsize=16)
    pyplot.savefig('%d-s.png' % time.time())
    return



if __name__ =='__main__':
    # circular = [0] + [(i-1)*pi/3.+pi/2 for i in range(1,7)]
    # print(circular)
    # gen_text(circular)
    anti_radial = [0] + [(i-1)*pi/3.+pi*(1-i%2) for i in range(1,7)]
    print(anti_radial)
    gen_text(anti_radial)
