import system
import matplotlib.animation as animation
import matplotlib.pyplot as pyplot
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import time

from numpy import pi, cos, sin

def frame_gen(state):
    axes.clear()
    data = state[0]
    config = state[1]
    arrow_patches = []
    colors = []
    for ind in range(7):
        phi = config[ind]
        tht0i = 2*pi*(ind-1)/6.
        if ind ==0:
            ci = [0.,0.]
        else:
            ci = [cos(tht0i), sin(tht0i)]
        arrow_patches.append(
            mpatches.Arrow(
                ci[0], ci[1], cos(phi)/2, sin(phi)/2,
                width=.2, zorder = 1.
                )
            )
        arrow_patches.append(
            mpatches.Circle(
                ci, radius=.5, alpha=.1, zorder = .5
                )
            )
        colors.append((.1,.5,.1))
        colors.append((.5,.1,.1,.1))

    shapes = PatchCollection(
        arrow_patches,
        facecolors=colors
        )
    axes.add_collection(shapes)
    axes.set_title('t=%.3f $\\alpha$=%.3f U=%.3f $\\mu$=%.3f' % data)

def animation_runner():
    limit = 2
    axes.set_ylim(-limit, limit)
    axes.set_xlim(-limit, limit)
    axes.set_aspect(aspect='equal')
    magnets = system.System()
    states = gen_path(magnets)
    states = states[::4]
    movie = animation.FuncAnimation(
        fig,
        frame_gen,
        frames = states
        )
    movie.save('%d.mp4' % time.time())
    return

def gen_path(magnets):
    tolerance = 10.**-9.
    t0 = 0.
    path = [snap_shot(magnets)]
    del_t= 0
    print(magnets.total_delta_sqr())
    while magnets.alpha < 3.5:
        while magnets.total_delta_sqr() > tolerance or del_t< 1.0:
            del_t = magnets.elapsed-t0
            if del_t > 1.:
                magnets.gamma = magnets.gamma = 2*(magnets.alpha+1)*(del_t-1.)/del_t
            magnets.advance_in_time()
            path.append(snap_shot(magnets))
        magnets.alpha+=.05
        magnets.gamma=0.
        t0 = magnets.elapsed
        magnets.advance_in_time()
        del_t = magnets.elapsed-t0
        path.append(snap_shot(magnets))
    return path

def snap_shot(magnets):
    print(magnets.alpha)
    return (
        (
            magnets.elapsed,
            magnets.alpha,
            magnets.total_PE(),
            magnets.net_dipole_mag()
            ),
        tuple(magnets.state[:7])
        )


if __name__ == '__main__':
    fig, axes = pyplot.subplots()
    animation_runner()
