import system
import numpy

def debugging1():
    magnets = system.System()
    # print(magnets.r_vals)
    # print(magnets.theta_vals)
    print(magnets.state)
    print(magnets.total_force())
    magnets.alpha = .5
    # print(magnets.kernels)
    print(magnets.total_force())
    magnets.rk45_step()
    print(magnets.delta_4)
    print(magnets.delta_5)
    # print(magnets.step_rescale())
    print(magnets.total_KE())
    for i in range(10):
        print(magnets.advance_in_time())
        print(magnets.state)
        print(magnets.total_force()[0],magnets.total_force()[7])
        print(magnets.total_force()[0+1],magnets.total_force()[7+1])
    magnets2 = system.System()
    # magnets2.alpha = .5
    # for i in range(40):
    #     print('tot_E', magnets2.total_KE()+magnets2.total_PE())
    #     magnets2.advance_in_time()
    # print(magnets2.elapsed)
    print(magnets.theta_vals[0])
    print(magnets.theta_vals[:,0])
    print(magnets.theta_vals[0] - magnets.theta_vals[:,0])


def central_bug():
    magnets = system.System()
    magnets.alpha = 0.
    magnets.state*=0.
    PE0 = magnets.total_PE()
    print(PE0)
    for i in range(7):
        # magnets.state[i]+= .01
        magnets.alpha += 0.2
        print(magnets.total_PE()-PE0)
        # magnets.state[i]-= .01

def debugging_2():
    magnets = system.System()
    magnets.load_state(.4)
    PE0 = magnets.total_PE()
    print(magnets.alpha)
    print(magnets.state[:7])
    print(PE0)
    magnets.gamma = 1.3
    for i in range(10):
        magnets.advance_in_time()
    print(magnets.state[:7])
    print(magnets.total_PE()-PE0)

    return

def debugging_3():
    magnets = system.System()
    magnets.load_state(2.46)
    mu0 = magnets.net_dipole_mag()
    print(magnets.net_dipole_moment())
    for i in range(3):
        magnets.shift_alpha_and_stablize(0.02)
        mu1 = magnets.net_dipole_mag()
        print(mu0, mu1)
        print(magnets.net_dipole_moment())

def debugging_4():
    # derivatives
    magnets = system.System()
    magnets.load_state(.5)
    step_sizes = numpy.array([.7**num for num in range(10)])
    for step in step_sizes:
        print('h = %f' % step)
        total = 0.
        for i in range(7):
            for j in range(7):
                DE = magnets.extended_diff(i, j, step)
                DC = magnets.central_diff(i, j, step)
                delt = abs(DE-DC)
                avg = (DE+DC)/2
                total += abs(delt/avg)
                # total += DE/DC
                if DE*DC<0.:
                    print('sign mismatch')
                # print(DE,DC)
        print(total/49.)
        # print(total/48.)
        print()

def debugging_5():
    magnets = system.System()
    magnets.load_state(.5)
    step_sizes = numpy.array([.1**num for num in range(4)])
    last = False
    biggest_shift = 0.
    for step in step_sizes:
        print(step)
        magnets.calc_L_mat(step)
        # print(last)
        if not last:
            last = tuple([tuple(el) for el in magnets.L_mat])
            continue
        for i in range(7):
            for j in range(7):
                delt = last[i][j]-magnets.L_mat[i][j]
                if delt > biggest_shift:
                    biggest_shift = delt
        print(biggest_shift)
        biggest_shift = 0.
        print(last[0][0], magnets.L_mat[0][0])
        last = tuple([tuple(el) for el in magnets.L_mat])

def debugging_6():
    magnets = system.System()
    # magnets.load_state(.5)
    U_0 = magnets.total_PE()
    magnets.calc_L_mat(.01)
    # U_A = magnets.total_PE()
    # print(U_A-U_0)
    template = " %.3f :"*7
    print('state')
    print(template % tuple(magnets.state[:7]))
    print('')
    # for line in magnets.L_mat:
    #     print(template % tuple(line))
    eig_stuff = numpy.linalg.eig(magnets.L_mat)
    tidy_eigs = magnets.spectrum_finder()
    for ind in range(len(eig_stuff[0])):
        print(eig_stuff[0][ind], tidy_eigs[ind][0])
        print(template % tuple(eig_stuff[1][:,ind]))
        print(template % tuple(tidy_eigs[ind][1]))
        print()
    # U_A = magnets.total_PE()
    # print(U_A-U_0)

    tidy_eigs = magnets.spectrum_finder()
    # for mode in tidy_eigs:
    #     print(mode[0])
    #     print(template % tuple(mode[1]))

    # U_A = magnets.total_PE()
    # print(U_A-U_0)

def debugging_labeler():
    magnets = system.System()
    # magnets.load_state(1.96)
    # modes = magnets.labeled_spectra()
    # for mode_id in modes.keys():
    #     print(mode_id)
    #     print(modes[mode_id])
    checks = [1.94+.01*i for i in range(10)]
    for alph in checks:
        print(alph)
        magnets.load_state(alph)
        modes = magnets.labeled_spectra()
        # print(modes)
        vec = modes[4][1]
        neg4 = magnets.sign_compare(vec)
        neg6 = magnets.sign_compare(modes[6][1])
        print('4', neg4)
        print('6', neg6[:3], neg6[3:])
        print('')


# central_bug()
# debugging_2()
# debugging_3()
# debugging_4()
# debugging_5()
# debugging_6()
debugging_labeler()
