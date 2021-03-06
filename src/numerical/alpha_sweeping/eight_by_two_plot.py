import system
from matplotlib import pyplot
from matplotlib import backend_bases
from matplotlib import transforms
from matplotlib import rcParams
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

def gen_fig_and_subplots(square=False):
    if not square:
        figure, subplots = pyplot.subplots(4,2)
        figure.set_figheight(24)
        figure.set_figwidth(12)
        return figure,subplots

    figure = pyplot.figure()
    figure.set_figheight(18)
    figure.set_figwidth(18)
    dimen = (6,6)
    subplots = []
    for row in range(3):
        row_arr = []
        for col in range(3):
            if row==2 and col==2:
                break
            ri, cj = 2*row, col*2
            if row==2:
                cj+=1
            print(row,col,ri,cj)
            row_arr.append(
                pyplot.subplot2grid(
                    dimen, (ri, cj),
                    rowspan=2,colspan=2
                    )
                )
        subplots.append(row_arr)
    return figure, subplots

def making_plot(poly=True):
    markers = [' ',' ',' ',' ', 'x', 'x', 'x']
    styles = ['-','--', '--','--', ' ',' ', ' ']
    colors = ['k','r', 'g', 'b', 'r', 'b','g']
    sub_label = 'a b c d a b c d'.split(' ')
    pyplot.rc('xtick', labelsize=14)
    pyplot.rc('ytick', labelsize=14)
    # figure, subplots = pyplot.subplots(4,2)
    # figure.set_figheight(24)
    # figure.set_figwidth(12)
    # tick marks
    #https://matplotlib.org/3.1.1/gallery/ticks_and_spines/tick-locators.html#sphx-glr-gallery-ticks-and-spines-tick-locators-py
    square_layout = False
    figure, subplots = gen_fig_and_subplots(square=square_layout)

    vector_yaxis_label = '$\\delta \\phi$'
    vector_title_template = '(%s) $\\omega_%d$'
    print(len(subplots))
    print(len(subplots[0]))

    if poly:
        print('poly mode')
        prefix = 'poly'
        x_axis_label = '$\\alpha$'
        y_axis_off = -.08 #how shifted is the y axis label
        y_axis_0 = 0.03
        x_axis_off =  -0.035 #how shifted is the x axis label
        eig_val_yaxis = '$\\omega_i^2$'
        eig_val_title = '(a) $\\omega^2$ vs. $\\alpha$'
        get_values = load_poly_values
        eig_limit = 10
        limit_sets = lambda x: (0, x[-1])
        annotate = poly_annotate
    else:
        print('mono mode')
        prefix = 'mono'
        x_axis_label = '$\\alpha^{-1}$'
        y_axis_off = -.09 #how shifted is the y axis label
        y_axis_0 = -0.00
        x_axis_off = -.055 #how shifted is the x axis label
        eig_val_yaxis = '$\\omega^2 \\alpha^{-1}$'
        eig_val_title = '(a) $\\omega^2 \\alpha^{-1}$ vs. $\\alpha^{-1}$'
        eig_limit = 4
        limit_sets = lambda x: (x[-1], x[0])
        get_values = load_mono_values
        annotate = mono_annotate

    x_values, eigen_vals, eigen_vecs = get_values()
    #x_values, 1-D list of alpha or 1/alpha
    #eigen values: same length as x_values, each entry 7 eigen values at that x-value
    # eigen_vecs: 7D list, one for each mode, x-value number of entries, of 7 vectors each

    x_min, x_max = limit_sets(x_values)
    yticks = (-1., -.5, 0, .5, 1.)
    fontsize=20
    subplots[0][0].set_xlabel(x_axis_label, fontsize=fontsize)
    subplots[0][0].set_ylabel(eig_val_yaxis, fontsize=fontsize)
    subplots[0][0].set_ylim(0, eig_limit)
    subplots[0][0].set_xlim(x_min, x_max)
    subplots[0][0].xaxis.set_label_coords(.5, x_axis_off)
    subplots[0][0].yaxis.set_label_coords(y_axis_off+y_axis_0, .5)
    subplots[0][0].tick_params(axis='both', labelsize=14)
    subplots[0][0].set_title(
        label = eig_val_title,
        loc='right', fontsize=fontsize)

    for i in range(7):
        y_curve = [entry[i] for entry in eigen_vals]
        subplots[0][0].plot(x_values, y_curve)
        if square_layout:
            row = int((i+1)/3)
            col = (i+1)%3
        else:
            row = int((i+1)/2)
            col = (i+1)%2
        print(i, row, col, len(subplots[row]))
        subplots[row][col].set_title(
            label=vector_title_template % (sub_label[i+1],(i+1)),
            loc='right', fontsize=fontsize)
        subplots[row][col].set_ylim(-1, 1)
        subplots[row][col].set_xlim(x_min, x_max)
        subplots[row][col].set_xlabel(x_axis_label, fontsize=fontsize)
        subplots[row][col].set_ylabel(vector_yaxis_label, fontsize=fontsize)
        subplots[row][col].xaxis.set_label_coords(.5, x_axis_off)
        subplots[row][col].yaxis.set_label_coords(y_axis_off, .5)
        subplots[row][col].set_yticks(yticks)

        subplots[row][col].axhline(y=0,color='k',linestyle='--')
        # subplots[row][col].grid(which='both', axis='x')


        for j in range(7):
            line = [entry[j] for entry in eigen_vecs[i]]
            subplots[row][col].plot(
                x_values, line,
                color = colors[j],
                linestyle=styles[j],
                marker=markers[j],
                markevery=3, markersize=3
                    )

    annotate(subplots, square=square_layout)
    time_label = ("%d" % time.time())[-5:]
    # pyplot.subplots_adjust(wspace=.45,hspace=.45)
    pyplot.tight_layout()
    fig_dims = figure.get_size_inches()
    epsi=0.
    topHalf = transforms.Bbox([[0,(1+epsi)*fig_dims[1]/2],[fig_dims[0],fig_dims[1]]])
    botHalf = transforms.Bbox([[0,0],[fig_dims[0],(1-epsi)*fig_dims[1]/2]])
    #https://matplotlib.org/3.3.3/api/transformations.html#matplotlib.transforms.Bbox
    print(topHalf, botHalf)
    # print(rcParams)
    pyplot.savefig(time_label+('-%s-kit_and_kabootleB.png' % prefix), bbox_inches=botHalf)
    for i in range(4,8):
        print(int(i/2), i%2)
        newt = subplots[int(i/2)][i%2].set_title(label=' ', loc='right')
        print(newt)
    pyplot.savefig(time_label+('-%s-kit_and_kabootleA.png' % prefix), bbox_inches=topHalf)
    # pyplot.savefig(time_label+('-%s-kit_and_kabootle_all.png' % prefix))

def poly_annotate(subplots, square=False):
    fsize = 18
    if square:
        e=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)]
    else:
        e=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1)]
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_1$', xy=(1.0,.6**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_2$', xy=(.75,1.52**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_3$', xy=(1.25,1.75**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_4$', xy=(.95,1.87**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_5$', xy=(1.0,2.2**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_6$', xy=(1.0,2.55**2),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_7$', xy=(.90,2.92**2),fontsize=fsize)


    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_0$', xy=(.05,.9),fontsize=fsize)
    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_5,\\phi_6$', xy=(.25,.25),fontsize=fsize)
    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(1.5,-0.2),fontsize=fsize)
    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_2, \\phi_3$', xy=(.5,-.3),fontsize=fsize)

    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_0$', xy=(.1,.05),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_5$', xy=(.5,.6),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_4$', xy=(1.0,.40),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_2$', xy=(.5,.2),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_3$', xy=(.5,-.2),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_1$', xy=(1.0,-.40),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_6$', xy=(.5,-.65),fontsize=fsize)

    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_0$', xy=(.9,.05),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_2$', xy=(.25,.60),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_6$', xy=(.25,.35),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_4$', xy=(.25,.15),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_1$', xy=(.25,-.20),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_5$', xy=(.25,-.40),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_3$', xy=(.25,-.7),fontsize=fsize)

    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_0$', xy=(.05,.075),fontsize=fsize)
    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_5,\\phi_6$', xy=(.25,.4),fontsize=fsize)
    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_2, \\phi_3$', xy=(.25,.15),fontsize=fsize)
    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(.25,-0.5),fontsize=fsize)

    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_0$', xy=(.1,.075),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_1$', xy=(.25,.65),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_2$', xy=(.35,.38),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_5$', xy=(.5,.11),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_6$', xy=(.5,-.16),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_3$', xy=(.35,-.43),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_4$', xy=(.25,-.75),fontsize=fsize)

    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_0$', xy=(.8,.5),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_2, \\phi_3$', xy=(.25,.5),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(.25,-0.22),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_5,\\phi_6$', xy=(.25,-.65),fontsize=fsize)

    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_0$', xy=(.05,.075),fontsize=fsize)
    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_2, \\phi_3$', xy=(.75,.65),fontsize=fsize)
    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(1.2,0.3),fontsize=fsize)
    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_5,\\phi_6$', xy=(1.2,-.3),fontsize=fsize)

def mono_annotate(subplots, square=False):
    fsize = 20
    if square:
        e=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1)]
    else:
        e=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1)]

    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_1$', xy=(.05,.26),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_2$', xy=(.15,.57),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_3$', xy=(.15,1.156),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_4$', xy=(.275,2.1),fontsize=fsize)
    # subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_5$', xy=(.025,1.67),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_5$', xy=(.5,1.67),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_6$', xy=(.025,2.1),fontsize=fsize)
    subplots[e[0][0]][e[0][1]].annotate(text='$\\omega_7 \\Uparrow$', xy=(.10,3.5),fontsize=fsize)

    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_0$', xy=(.1,.5),fontsize=fsize)
    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(.1,-0.2),fontsize=fsize)
    subplots[e[1][0]][e[1][1]].annotate(text='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.1,-.6),fontsize=fsize)

    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_0$', xy=(.4,.05),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_2,\\phi_5$', xy=(.2,.6),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_1$', xy=(.05,.25),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_4$', xy=(.05,-.30),fontsize=fsize)
    subplots[e[2][0]][e[2][1]].annotate(text='$\\phi_3,\\phi_6$', xy=(.2,-.63),fontsize=fsize)

    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_0, \\phi_1, \\phi_4$', xy=(.1,.1),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_2, \\phi_6$', xy=(.1,0.6),fontsize=fsize)
    subplots[e[3][0]][e[3][1]].annotate(text='$\\phi_3, \\phi_5$', xy=(.1,-.63),fontsize=fsize)

    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_0, \\phi_1, \\phi_4$', xy=(.1,.1),fontsize=fsize)
    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_2, \\phi_3$', xy=(.1,0.6),fontsize=fsize)
    subplots[e[4][0]][e[4][1]].annotate(text='$\\phi_5, \\phi_6$', xy=(.1,-.63),fontsize=fsize)

    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_0$', xy=(.1,.15),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(.1,0.8),fontsize=fsize)
    subplots[e[5][0]][e[5][1]].annotate(text='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.3,-.4),fontsize=fsize)

    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_0$', xy=(.3,.05),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_1$', xy=(.1,.80),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_3,\\phi_6$', xy=(.1,.3),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_2,\\phi_5$', xy=(.1,-.35),fontsize=fsize)
    subplots[e[6][0]][e[6][1]].annotate(text='$\\phi_4$', xy=(.1,-.85),fontsize=fsize)

    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_0$', xy=(.1,.85),fontsize=fsize)
    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_2, \\phi_3,\\phi_5,\\phi_6$', xy=(.1,.3),fontsize=fsize)
    subplots[e[7][0]][e[7][1]].annotate(text='$\\phi_1, \\phi_4$', xy=(.1,-0.1),fontsize=fsize)

if __name__ == '__main__':
    # making_plot()
    making_plot(poly=False)
