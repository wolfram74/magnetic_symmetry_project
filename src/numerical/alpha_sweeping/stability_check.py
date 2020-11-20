import system
import numpy
import random
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()

def random_perturb():
    delta = numpy.zeros(14)
    for i in range(7):
        delta[i] = random.random()
    mag = numpy.linalg.norm(delta)
    # print(mag)
    delta = delta/mag
    return delta

def try_nudge():
    # magnets.load_state(2.05)
    # gam_p = tuple(magnets.state)
    magnets.load_state(2.05, down=True)
    gam0 = tuple(magnets.state)
    #how big is the difference between poly and mono phase at this alpha?
    # state_del = numpy.array(gam_p)-numpy.array(gam0)
    # mag_state_delta = numpy.linalg.norm(state_del)
    # print(state_del)
    # print(mag_state_delta) ~2.656
    nudge = random_perturb()
    magnets.state+=nudge*.2
    magnets.shift_alpha_and_stablize(0.0)
    deviation = numpy.linalg.norm(gam0-magnets.state)
    print(deviation)
    return deviation

if __name__=='__main__':
    # delta = random_perturb()
    # print(delta)
    deviations = []
    for i in range(100):
        deviations.append(try_nudge())
    sorted_devs = sorted(deviations)
    print(sorted_devs[-10:])

