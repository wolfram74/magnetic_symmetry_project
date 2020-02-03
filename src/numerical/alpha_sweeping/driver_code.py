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
    step_sizes = numpy.array([.8**num for num in range(10)])
    for step in step_sizes:
        print('h = %f' % step)

# central_bug()
# debugging_2()
# debugging_3()
debugging_4()
