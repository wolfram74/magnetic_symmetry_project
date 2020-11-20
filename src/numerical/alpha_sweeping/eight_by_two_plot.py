import system
from matplotlib import pyplot
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import time

import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()

def vectors_strip(array):
    #array of w^2 vector pairs
    #want 7 arrays of w^2 coord pairs
    output = [[] for i in range(8)]
    for entry in array:
        output[0].append(entry[0])

        for i in range(7):
            output[i+1].append(entry[i])
    return output

def deltA_from_deltU(alpha, deltU):
    return -deltU*alpha**2/(1.+deltU*alpha)



def load_poly_values():
    source = 'poly'
    template = './saved_eigenmodes/'+source+'_%d_mode.txt'
    alphas = []
    e_vals = []
    eigen_drifting = [[] for i in range(7)]
    eigen_vals = [[] for i in range(7)]
    for i in range(7):
        address = template % (i+1)
        cache = open(address, 'r')
        for line in cache:
            split = line.rstrip().split(',')
            alpha = float(split[0])
            val = float(split[1])
            # print(split[2:])
            vector = [float(el) for el in split[2:-1]]
            if i==0:
                alphas.append(alpha)
            eigen_vals[i].append(val)
            eigen_drifting[i].append(vector)
        cache.close()
    for ind in range(len(eigen_vals[0])):
        e_vals.append([mode[ind] for mode in eigen_vals])

    magnets.load_state(2.47)
    magnets.shift_alpha_and_stablize(.00278)
    alphas.append(magnets.alpha)
    tidy_eigs = magnets.labeled_spectra()
    eigs = [tidy_eigs[i+1][0] for i in range(7)]
    e_vals.append(eigs)
    for i in range(7):
        eigen_drifting[i].append(tidy_eigs[i+1][1])
    # print(len(alphas),len(e_vals))
    print(len(alphas),len(eigen_drifting), len(eigen_drifting[0]))
    return alphas, e_vals, eigen_drifting

def load_mono_values():
    #x_values, 1-D list of alpha or 1/alpha
    #eigen values: same length as x_values, each entry 7 eigen values at that x-value
    # eigen_vecs: 7D list, one for each mode, x-value number of entries, of 7 vectors each
    inv_alphas = []
    e_vals = []
    eigen_drifting = [[] for i in range(7)]

    running = True
    last = False
    u_min = 10**-3
    deltU= -.01
    magnets.load_state(1.15, down=True)
    u = 1/magnets.alpha
    OMEGAS = [0.1815, 1.3229, 1.3229, 2.0000, 1.9126, 1.3229, 10.5203]


    while running:
        inv_alphas.append(1/magnets.alpha)
        tidy_eigs = magnets.labeled_spectra_mono()
        omegs = [tidy_eigs[i][0]/magnets.alpha for i in range(1,8)]
        e_vals.append(omegs)
        for i in range(7):
            eigen_drifting[i].append(tidy_eigs[i+1][1])
        if u < u_min and last:
            break
        if magnets.alpha < 5.:
            delta_alph = deltA_from_deltU(magnets.alpha, deltU)
        if magnets.alpha < 5. and magnets.alpha+delta_alph < 2.5:
            magnets.load_state(magnets.alpha+delta_alph, down=True)
        else:
            magnets.load_state(magnets.alpha*1.151)
        # magnets.shift_alpha_and_stablize(deltA_from_deltU(magnets.alpha, deltU))
        u = 1/magnets.alpha
        if u < u_min:
            last=True

    return inv_alphas, e_vals, eigen_drifting

def making_plot(poly=True):
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    sub_label = 'a b c d e f g h'.split(' ')
    figure, subplots = pyplot.subplots(4,2)
    figure.set_figheight(24)
    figure.set_figwidth(12)
    vector_yaxis_label = '$\\delta \\phi$'
    vector_title_template = '(%s) $\\omega_%d$'
    print(len(subplots))
    print(len(subplots[0]))

    if poly:
        print('poly mode')
        x_axis_label = '$\\alpha$'
        eig_val_yaxis = '$\\omega_i^2$'
        eig_val_title = '(a) $\\omega^2$ vs $\\alpha$'
        get_values = load_poly_values
        eig_limit = 10
        limit_sets = lambda x: (x[0], x[-1])
        annotate = poly_annotate
    else:
        print('mono mode')
        x_axis_label = '$\\alpha^{-1}$'
        eig_val_yaxis = '$\\omega^2 \\alpha^{-1}$'
        eig_val_title = '(a) $\\omega^2 \\alpha^{-1}$ vs $\\alpha^{-1}$'
        eig_limit = 4
        limit_sets = lambda x: (x[-1], x[0])
        get_values = load_mono_values
        annotate = mono_annotate

    x_values, eigen_vals, eigen_vecs = get_values()
    #x_values, 1-D list of alpha or 1/alpha
    #eigen values: same length as x_values, each entry 7 eigen values at that x-value
    # eigen_vecs: 7D list, one for each mode, x-value number of entries, of 7 vectors each
    # print(len(x_values))
    # print(len(eigen_vals), len(eigen_vals[0]))
    # print(len(eigen_vecs), len(eigen_vecs[0]))
    # print(eigen_vecs[0][0])
    x_min, x_max = limit_sets(x_values)
    subplots[0][0].set_xlabel(x_axis_label, fontsize=16)
    subplots[0][0].set_ylabel(eig_val_yaxis)
    subplots[0][0].set_ylim(0, eig_limit)
    subplots[0][0].set_title(
        label = eig_val_title,
        loc='right')

    for i in range(7):
        y_curve = [entry[i] for entry in eigen_vals]
        subplots[0][0].plot(x_values, y_curve)
        row = (i+1)/2
        col = (i+1)%2
        print(row, col)
        subplots[row][col].set_title(
            label=vector_title_template % (sub_label[i+1],(i+1)),
            loc='right')
        subplots[row][col].set_ylim(-1, 1)
        subplots[row][col].set_xlim(x_min, x_max)
        subplots[row][col].set_xlabel(x_axis_label, fontsize=16)
        subplots[row][col].set_ylabel(vector_yaxis_label, fontsize=16)
        subplots[row][col].xaxis.set_label_coords(.5, -0.035)
        subplots[row][col].yaxis.set_label_coords(-.0675, 0.5)

        subplots[row][col].axhline(y=0,color='k',linestyle='--')
        subplots[row][col].grid(which='both', axis='x')


        for j in range(7):
            line = [entry[j] for entry in eigen_vecs[i]]
            subplots[row][col].plot(
                x_values, line,
                color = colors[j],
                linestyle=styles[j],
                marker=markers[j],
                markevery=3, markersize=3
                )

    annotate(subplots)
    time_label = ("%d" % time.time())[-5:]
    pyplot.savefig(time_label+'-kit_and_kabootle.png', bbox_inches='tight')

def poly_annotate(subplots):
    subplots[0][0].annotate(s='$\\omega_1$', xy=(1.0,.5**2))
    subplots[0][0].annotate(s='$\\omega_2$', xy=(1.0,1.52**2))
    subplots[0][0].annotate(s='$\\omega_3$', xy=(1.0,1.72**2))
    subplots[0][0].annotate(s='$\\omega_4$', xy=(1.0,1.87**2))
    subplots[0][0].annotate(s='$\\omega_5$', xy=(1.0,2.2**2))
    subplots[0][0].annotate(s='$\\omega_6$', xy=(1.0,2.55**2))
    subplots[0][0].annotate(s='$\\omega_7$', xy=(1.0,2.92**2))


    subplots[0][1].annotate(s='$\\phi_0$', xy=(.05,.9))
    subplots[0][1].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,.25))
    subplots[0][1].annotate(s='$\\phi_1, \\phi_4$', xy=(1.0,0.1))
    subplots[0][1].annotate(s='$\\phi_2, \\phi_3$', xy=(.5,-.3))

    subplots[1][0].annotate(s='$\\phi_0$', xy=(.1,.05))
    subplots[1][0].annotate(s='$\\phi_5$', xy=(.5,.6))
    subplots[1][0].annotate(s='$\\phi_4$', xy=(1.0,.40))
    subplots[1][0].annotate(s='$\\phi_2$', xy=(.5,.2))
    subplots[1][0].annotate(s='$\\phi_3$', xy=(.5,-.2))
    subplots[1][0].annotate(s='$\\phi_1$', xy=(1.0,-.40))
    subplots[1][0].annotate(s='$\\phi_6$', xy=(.5,-.65))

    subplots[1][1].annotate(s='$\\phi_0$', xy=(.9,.05))
    subplots[1][1].annotate(s='$\\phi_2$', xy=(.25,.60))
    subplots[1][1].annotate(s='$\\phi_6$', xy=(.25,.35))
    subplots[1][1].annotate(s='$\\phi_4$', xy=(.25,.15))
    subplots[1][1].annotate(s='$\\phi_1$', xy=(.25,-.20))
    subplots[1][1].annotate(s='$\\phi_5$', xy=(.25,-.40))
    subplots[1][1].annotate(s='$\\phi_3$', xy=(.25,-.7))

    subplots[2][0].annotate(s='$\\phi_0$', xy=(.05,.075))
    subplots[2][0].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,.4))
    subplots[2][0].annotate(s='$\\phi_2, \\phi_3$', xy=(.25,.15))
    subplots[2][0].annotate(s='$\\phi_1, \\phi_4$', xy=(.25,-0.5))

    subplots[2][1].annotate(s='$\\phi_0$', xy=(.1,.1))
    subplots[2][1].annotate(s='$\\phi_1$', xy=(.25,.65))
    subplots[2][1].annotate(s='$\\phi_2$', xy=(.35,.38))
    subplots[2][1].annotate(s='$\\phi_5$', xy=(.25,.13))
    subplots[2][1].annotate(s='$\\phi_6$', xy=(.25,-.18))
    subplots[2][1].annotate(s='$\\phi_3$', xy=(.35,-.43))
    subplots[2][1].annotate(s='$\\phi_4$', xy=(.25,-.75))

    subplots[3][0].annotate(s='$\\phi_0$', xy=(.8,.5))
    subplots[3][0].annotate(s='$\\phi_2, \\phi_3$', xy=(.25,.5))
    subplots[3][0].annotate(s='$\\phi_1, \\phi_4$', xy=(.25,-0.22))
    subplots[3][0].annotate(s='$\\phi_5,\\phi_6$', xy=(.25,-.65))

    subplots[3][1].annotate(s='$\\phi_0$', xy=(.05,.075))
    subplots[3][1].annotate(s='$\\phi_2, \\phi_3$', xy=(.75,.65))
    subplots[3][1].annotate(s='$\\phi_1, \\phi_4$', xy=(1.2,0.3))
    subplots[3][1].annotate(s='$\\phi_5,\\phi_6$', xy=(1.2,-.3))

def mono_annotate(subplots):
    subplots[0][0].annotate(s='$\\omega_1$', xy=(.05,.25))
    subplots[0][0].annotate(s='$\\omega_2$', xy=(.1,.85))
    subplots[0][0].annotate(s='$\\omega_3$', xy=(.1,1.176))
    subplots[0][0].annotate(s='$\\omega_4$', xy=(.05,1.55))
    subplots[0][0].annotate(s='$\\omega_5$', xy=(.05,1.725))
    subplots[0][0].annotate(s='$\\omega_6$', xy=(.05,1.95))
    subplots[0][0].annotate(s='$\\omega_7 \Uparrow$', xy=(.10,3.5))

    subplots[0][1].annotate(s='$\\phi_0$', xy=(.1,.5))
    subplots[0][1].annotate(s='$\\phi_1, \\phi_4$', xy=(.1,-0.2))
    subplots[0][1].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.1,-.6))

    subplots[1][0].annotate(s='$\\phi_0$', xy=(.4,.05))
    subplots[1][0].annotate(s='$\\phi_2,\\phi_5$', xy=(.2,.6))
    subplots[1][0].annotate(s='$\\phi_1$', xy=(.05,.25))
    subplots[1][0].annotate(s='$\\phi_4$', xy=(.05,-.30))
    subplots[1][0].annotate(s='$\\phi_3,\\phi_6$', xy=(.2,-.63))

    subplots[1][1].annotate(s='$\\phi_0, \\phi_1, \\phi_4$', xy=(.1,.1))
    subplots[1][1].annotate(s='$\\phi_2, \\phi_6$', xy=(.1,0.6))
    subplots[1][1].annotate(s='$\\phi_3, \\phi_5$', xy=(.1,-.63))

    subplots[2][0].annotate(s='$\\phi_0, \\phi_1, \\phi_4$', xy=(.1,.1))
    subplots[2][0].annotate(s='$\\phi_2, \\phi_3$', xy=(.1,0.6))
    subplots[2][0].annotate(s='$\\phi_5, \\phi_6$', xy=(.1,-.63))

    subplots[2][1].annotate(s='$\\phi_0$', xy=(.1,.15))
    subplots[2][1].annotate(s='$\\phi_1, \\phi_4$', xy=(.1,0.8))
    subplots[2][1].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.3,-.4))

    subplots[3][0].annotate(s='$\\phi_0$', xy=(.3,.05))
    subplots[3][0].annotate(s='$\\phi_1$', xy=(.1,.80))
    subplots[3][0].annotate(s='$\\phi_3,\\phi_6$', xy=(.1,.3))
    subplots[3][0].annotate(s='$\\phi_2,\\phi_5$', xy=(.1,-.35))
    subplots[3][0].annotate(s='$\\phi_4$', xy=(.1,-.85))

    subplots[3][1].annotate(s='$\\phi_0$', xy=(.1,.85))
    subplots[3][1].annotate(s='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.1,.3))
    subplots[3][1].annotate(s='$\\phi_1, \\phi_4$', xy=(.1,-0.1))

if __name__ == '__main__':
    # making_plot()
    making_plot(poly=False)
