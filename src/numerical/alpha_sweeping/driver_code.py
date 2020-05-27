import system
import numpy

def rho_n(n,N=6):
    angle_quant = numpy.pi/N
    rho_denom = numpy.sin(angle_quant)
    return numpy.sin(angle_quant*n)/rho_denom

def stump_ZN(N):
    ZN = 0.
    angle_quant = numpy.pi/N
    rho_denom = numpy.sin(angle_quant)
    # print(angle_quant, rho_denom)
    for n in range(1, N):
        rho = rho_n(n, N=N)
        numer = 1.+numpy.cos(angle_quant*n)**2
        ZN += numer/rho**3
    return ZN

def stump_mode_i(i, N=6):
    ZN = stump_ZN(N)
    angle_quant = numpy.pi/N
    lam_i = ZN
    # print(i)
    for nu in range(1,N):
        rho = rho_n(nu, N=N)
        numer = 1.+numpy.sin(angle_quant*nu)**2
        cos_term = numpy.cos(2*angle_quant*nu*i)
        lam_i += numer*cos_term/rho**3
    return lam_i

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
    checks = [1.55-.01*i for i in range(4)]
    for alph in checks:
        print(alph)
        magnets.load_state(alph, down=True)
        modes = magnets.labeled_spectra()
        # print(modes)
        for val in modes.keys():
            print(val, modes[val])
        # vec = modes[4][1]
        # neg4 = magnets.sign_compare(vec)
        # neg6 = magnets.sign_compare(modes[6][1])
        # print('4', neg4)
        # print(vec)
        # print('6', neg6[:3], neg6[3:])
        # print(modes[6][1])
        # print('')

def writing_more_sig_figs():
    magnets = system.System()
    magnets.load_state(.5)
    eigen_states = magnets.labeled_spectra()
    print(eigen_states.keys())
    address_template = './saved_eigenmodes/poly_%d_mode.txt'
    destinations = []
    line_template = '%.9f,'+'%.9f,'+'%.9f,'*7
    for key in eigen_states.keys():
        destinations.append(open(address_template%key, 'w'))

    for key in eigen_states.keys():
        data = [magnets.alpha, eigen_states[key][0]]+list(eigen_states[key][1])
        destinations[key-1].write(line_template % tuple(data))

def low_alpha_mode1():
    magnets= system.System()
    magnets.load_state(.18)
    modes = magnets.labeled_spectra()
    print(modes[1])
    magnets.shift_alpha_and_stablize(.01)
    modes = magnets.labeled_spectra()
    print(modes[1])

def high_alpha_modes():
    magnets= system.System()
    magnets.load_state(1000.)
    modes = magnets.labeled_spectra_mono()
    print(magnets.alpha)
    for key in modes.keys():
        print(key)
        print(modes[key][0]/magnets.alpha)
        print(modes[key][1])
    # print(modes)
    # magnets.load_state(2.52,down=True)
    # modes = magnets.labeled_spectra_mono()
    # print(magnets.alpha)
    # for key in modes.keys():
    #     print(key)
    #     print(modes[key])
    # magnets.shift_alpha_and_stablize(.7)
    # modes = magnets.labeled_spectra_mono()
    # print(modes)
    # magnets.shift_alpha_and_stablize(.2)
    # modes = magnets.labeled_spectra_mono()
    # print(modes)
    # magnets.shift_alpha_and_stablize(.2)
    # modes = magnets.labeled_spectra_mono()
    # print(modes)

def mode2vs4():
    magnets= system.System()
    magnets.load_state(3.5)
    modes = magnets.labeled_spectra_mono()
    vec2 = modes[2][1]
    vec4 = modes[4][1]
    abs2 = numpy.abs(vec2)
    abs4 = numpy.abs(vec4)
    print(sum(abs2))
    print(sum(abs4))

def limit_modes():
    magnets= system.System()
    magnets.load_state(3.5)
    magnets.shift_alpha_and_stablize(6.5)
    for i in range(9):
        print(i+1)
        magnets.shift_alpha_and_stablize(10.)
    modes = magnets.labeled_spectra_mono()
    print(magnets.alpha)
    for key in modes.keys():
        print('mode %d' % key)

        omega = modes[key][0]
        print("omega^2: %.9f , omega^2/alpha: %.9f" % (omega, omega/magnets.alpha))
        vec_template = '%.6f, '*7
        print('vector:')
        print(vec_template% tuple(modes[key][1]))


def compare_with_sump():
    # for i in range(3, 11):
    #     ZN = stump_ZN(i)
    #     print(i, ZN)
    N = 6
    omeg_2 = []
    for i in range(N):
        omeg_2.append(stump_mode_i(i+1, N=N))
    omeg_2 = sorted(omeg_2)
    print(omeg_2)
    w1, w2, w3, w4 = omeg_2[0],omeg_2[1],omeg_2[4],omeg_2[5]
    ratios = [w2/w1, w3/w2, w4/w3]
    print(ratios)
    print("%04f,"*3 % tuple(ratios))
    magnets = system.System()
    naive_spectra = magnets.spectrum_finder()
    print(naive_spectra)
    e_vals = [term[0] for term in naive_spectra]
    print(e_vals)
    e_vals = sorted(e_vals)
    w1, w2, w3, w4 = e_vals[0+1],e_vals[1+1],e_vals[4+1],e_vals[5+1]
    ratios2 = [w2/w1, w3/w2, w4/w3]
    print(ratios2)
    print("%04f,"*3 % tuple(ratios2))
    deltas = [ratios2[i]-ratios[i] for i in range(3)]
    print(deltas)


# central_bug()
# debugging_2()
# debugging_3()
# debugging_4()
# debugging_5()
# debugging_6()
# debugging_labeler()
# writing_more_sig_figs()
# low_alpha_mode1()
high_alpha_modes()
# mode2vs4()
# limit_modes()
# compare_with_sump()
