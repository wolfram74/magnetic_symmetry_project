--import system
from matplotlib import pyplot
import time

# from numpy import pi
import numpy
pi = numpy.pi

def gen_text(start_state):
    magnets = system.System()
    equilibria = []
    running = True
    # running = False
    monitor = 10.
    magnets.state *= 0.
    magnets.alpha = 0.
    magnets.state[:7] = start_state
    forward = True
    print(magnets.state)
    while running:
        magnets.advance_in_time()
        if magnets.elapsed >5.:
            magnets.gamma = 2*(magnets.alpha+1)*(magnets.elapsed-5.)/magnets.elapsed
        if magnets.total_delta_sqr()>1*10**-8 or magnets.elapsed<10.:
            if magnets.elapsed > monitor:
                monitor+=10
                print('monitor at', magnets.elapsed, magnets.alpha)
                # print(magnets.gamma, magnets.step_size, magnets.total_delta_sqr())
                # print(magnets.state[:7])
                # print(magnets.state[7:])
                print('\n')
            continue
        # print('done at %f' % magnets.elapsed)
        # print(magnets.total_delta_sqr(), magnets.total_PE())
        # print(magnets.state[:7])
        # print('\n')
        if magnets.alpha >= 0:
            equilibria.append(
                [magnets.alpha,
                tuple(magnets.state[:7]),
                magnets.lim_moment(), magnets.lim_U()
                # magnets.net_dipole_mag(), magnets.total_PE()
                ]
                )
        # if magnets.alpha < 2.459:
        #     magnets.alpha += .01
        # else:
        #     magnets.alpha += .0001
        # if forward:
        magnets.alpha+=.05
        # else:
        #     magnets.alpha-=.05
        # print(magnets.alpha)
        # print(magnets.total_force())
        magnets.elapsed = 0.
        magnets.gamma = 0.
        monitor = 10
        if magnets.alpha >= 2.8:
            # forward = False
            running = False
        # if not forward and magnets.alpha <.1:
        #     running = False
    angle_plotter(equilibria)
    # print(magnets.lim_moment(), magnets.lim_U())
    # state_plotter(equilibria)

def angle_plot_from_txt():
    state_data = open('stored_states.txt', 'r')
    equilibs = []
    for line in state_data:
        txt_list = line.split(',')
        # print(txt_list)
        txt_list.pop()
        num_list = [float(el) for el in txt_list]
        equilibs.append([num_list[0], num_list[1:]])
    angle_plotter(equilibs)

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
    ax2 = pyplot.subplot(212)
    ax2.set_ylim([.9,1.4])
    ax2.set_xlim([alpha[0],alpha[-1]])
    mu = [el[2] for el in equilibs]
    PE = [el[3] for el in equilibs]
    print(len(mu))
    subplots[0].plot(alpha, mu, label='mu')
    ax2.plot(alpha, PE, label='PE/lim_U')
    # ax2.plot(alpha, PE, label='PE')
    pyplot.legend()
    pyplot.xlabel('$\\alpha$', fontsize=16)
    pyplot.savefig('%d-s.png' % time.time())
    return

def down_slope_states():
    magnets = system.System()
    start = 1.1
    end = 1.09
    step = .00001
    magnets.load_state(start, down=True)
    print(magnets.state)
    alpha = [magnets.alpha]
    mux_vals = [magnets.net_dipole_moment()[0]]
    U_vals = [magnets.total_PE()]
    while magnets.alpha > end:
        print(magnets.alpha)
        if magnets.alpha-step <= 0.:
            break
        if magnets.alpha <= end:
            break
        if step >= .011:
            magnets.load_state(magnets.alpha-.0101, down=True)
        else:
            magnets.shift_alpha_and_stablize(-step)
        alpha.append(magnets.alpha)
        mux_vals.append(magnets.net_dipole_moment()[0])
        U_vals.append(magnets.total_PE())

    figure, subplots = pyplot.subplots(2,1)
    figure.set_figheight(8)
    figure.set_figwidth(12)
    subplots[0].plot(alpha,mux_vals, label='$\mu_x$')
    subplots[0].legend()
    subplots[1].plot(alpha,U_vals, label='$PE$')
    subplots[1].legend()
    pyplot.xlabel('$\\alpha$', fontsize=16)
    pyplot.savefig('%d-states.png' % time.time())


def dipole_components_plot():
    figure, subplots = pyplot.subplots(3,1)
    magnets = system.System()
    dipoles = [(
        magnets.alpha,
        magnets.net_dipole_moment(),
        magnets.net_dipole_mag()
        )]
    target = 0.
    while magnets.alpha < 3.0:
        target+= .01
        magnets.load_state(target)
        # magnets.shift_alpha_and_stablize(0.01)
        print(magnets.alpha)
        dipoles.append((
            magnets.alpha,
            magnets.net_dipole_moment(),
            magnets.net_dipole_mag()
        ))
    while magnets.alpha > .5:
        magnets.shift_alpha_and_stablize(-0.01)
        print(magnets.alpha)
        dipoles.append((
            magnets.alpha,
            magnets.net_dipole_moment(),
            magnets.net_dipole_mag()
        ))

    alphas = [el[0] for el in dipoles]
    mu_x = [el[1][0] for el in dipoles]
    mu_y = [el[1][1] for el in dipoles]
    mu = [el[2] for el in dipoles]
    subplots[0].plot(alphas, mu_x, label='x')
    subplots[1].plot(alphas, mu_y, label='y')
    subplots[2].plot(alphas, mu, label='mag')
    pyplot.savefig('%d-mu.png' % time.time())

    return

def eigen_value_plots():
    magnets = system.System()
    magnets.load_state(1.2, down=True)
    vals = [el[0] for el in magnets.spectrum_finder()]
    eigen_vals = [(magnets.alpha, sorted(vals))]
    print(eigen_vals[0])
    while magnets.alpha > 1.05:
        magnets.shift_alpha_and_stablize(-.0001)
        vals = [el[0] for el in magnets.spectrum_finder()]
        eigen_vals.append( (magnets.alpha, sorted(vals)) )
    alphas = [el[0] for el in eigen_vals]
    # alphas = numpy.linspace(0, 2, 50)
    # print(eigen_vals[-1])
    # print(alphas[:5])
    # print(alphas[-5:])
    figure, subplots = pyplot.subplots(1,1)
    figure.set_figheight(6)
    figure.set_figwidth(12)
    subplots.set_ylim((-.5,6.))
    for i in range(7):
        eigs = [el[1][i] for el in eigen_vals]
        subplots.plot(alphas, eigs)
    subplots.plot(alphas, numpy.zeros(len(alphas)),linestyle='--')
    pyplot.savefig('%d-eigs.png' % time.time())

def vector_drift():
    def magnitude(vec, i):
        return sum([vec[ind]*alpha_0[i][ind] for ind in range(7)])
    def massage(drift_vals):
        output = [drift_vals[0]]
        for i in range(len(drift_vals)-2):
            ind = i+1
            left_bad = drift_vals[ind]*drift_vals[0] < 0.

            if left_bad:
                output.append(drift_vals[ind]*-1.)
                continue
            output.append(drift_vals[ind])
        output.append(drift_vals[-1])
        return output
    magnets = system.System()
    # magnets.load_state(.4)
    vals = magnets.spectrum_finder()
    vals = sorted(vals, key=lambda el : el[0])
    alpha_0 = [el[1] for el in vals]
    eigen_drifts = [[] for el in range(7)]
    print(eigen_drifts)
    figure, subplots = pyplot.subplots(7,1)
    figure.set_figheight(28)
    figure.set_figwidth(12)
    alphas = []
    while magnets.alpha < 3.0:
        magnets.load_state(magnets.alpha+.0101)
        vals = magnets.spectrum_finder()
        vals = sorted(vals, key=lambda el : el[0])
        for ind in range(7):
            eigen_drifts[ind].append(vals[ind][1])
        alphas.append(magnets.alpha)
        drifts = [magnitude(vals[i][1], i) for i in range(7)]
        # drifts = vals[6][1]
        # eigen_drifts.append( (magnets.alpha, drifts) )
    # # print(eigen_drifts[-1])
    # # print(alphas[:5])
    # # print(alphas[-5:])
    for j in range(7):
        for i in range(7):
            eigs = [el[i] for el in eigen_drifts[j]]
            subplots[j].plot(alphas, eigs)
    pyplot.savefig('%d-drift.png' % time.time())

def alignment_drift():
    def magnitude(vec, i):
        return sum([vec[ind]*alpha_0[i][ind] for ind in range(7)])
    def massage(drift_vals):
        output = [drift_vals[0]]
        for i in range(len(drift_vals)-2):
            ind = i+1
            left_bad = drift_vals[ind]*drift_vals[0] < 0.

            if left_bad:
                output.append(drift_vals[ind]*-1.)
                continue
            output.append(drift_vals[ind])
        output.append(drift_vals[-1])
        return output
    magnets = system.System()
    # magnets.load_state(.4)
    vals = magnets.spectrum_finder()
    vals = sorted(vals, key=lambda el : el[0])
    alpha_0 = [el[1] for el in vals]
    eigen_drifts = []
    print(eigen_drifts)
    figure, subplots = pyplot.subplots(1,1)
    figure.set_figheight(4)
    figure.set_figwidth(12)
    alphas = []
    while magnets.alpha < 3.0:
        magnets.load_state(magnets.alpha+.0101)
        vals = magnets.spectrum_finder()
        vals = sorted(vals, key=lambda el : el[0])
        eigen_drifts.append(
            [magnitude(vals[ind][1], ind) for ind in range(7)]
            )
        alphas.append(magnets.alpha)
    for i in range(7):
        eigs = [el[i] for el in eigen_drifts]
        # print(len(alphas), len(eigs))
        subplots.plot(alphas, eigs)
    pyplot.savefig('%d-drift-align.png' % time.time())



if __name__ =='__main__':
    # circular = [0] + [(i-1)*pi/3.+pi/2 for i in range(1,7)]
    # print(circular)
    # gen_text(circular)
    # anti_radial = [0] + [(i-1)*pi/3.+pi*(1-i%2) for i in range(1,7)]
    # print(anti_radial)
    # gen_text(anti_radial)
    # angle_plot_from_txt()
    # dipole_components_plot()
    eigen_value_plots()
    # vector_drift()
    # alignment_drift()
    # down_slope_states()
