"""
interested in plotting vector field of our set of magnets.
salient facts:
    physical extent of the magnets only goes out to 1.5, plot /atleast/ to 2, probably 3.
    going through the plane will lead to a lot of singularities and experiments take place z ~= 1
    
aesthetic choices:
    quantities are measured relative to the field strength of a unit strength dipole at it's north pole IE
    vec B = [3(vec mu dot vec rhat)vec rhat - vec mu] *(2 r**3)**-1 
    don't have to care about permetivity, or other such dimensions
    
process
    0) generate vector fields of a single dipole to get a feel for what would be demonstrative
    1) calculate B fields of whole system at 100/unit length from -3 to 3, 3600 total data points
    2) cache data for furhter use to avoid redundant calculations
    3) start refining plotting decisions  
"""

cached_data = './cached_plots_and_data/B_field_vectors'
image_destination = './cached_plots_and_data/B_field_plots'

import system
from matplotlib import pyplot
from matplotlib import backend_bases
from matplotlib import transforms
from matplotlib import rcParams
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

import time
import os

import numpy
pi, cos, sin = numpy.pi, numpy.cos, numpy.sin

magnets = system.System()

class FieldPlotter():
    def __init__(self):
        self.magnets = system.System()
        alpha = 2.45
        mono = False
        self.magnets.load_state(desired_alpha=alpha, down=mono)
        
        figure, subplots = pyplot.subplots(1,1)
        self.figure = figure
        self.subplots = subplots

        self.determine_prefix()
        image_name = self.prefix+'simple_fields.png'

        self.gen_field_data()


    def determine_prefix(self):
        old_files = os.listdir(image_destination)
        next_num = len(old_files)+1
        self.prefix = "%04d_" % next_num

    def gen_field_data(self, N=30, z=1., size=2.):
        Xs = numpy.linspace(-size, size, num=points_per_line, endpoint=True)
        Ys = numpy.linspace(-size, size, num=points_per_line, endpoint=True)
        z = 1.
        self.plotX = []
        self.plotY = []
        self.BX = []
        self.BY = []
        self.B_vals = []
        self.colors = []
        self.largest_B = 0


def plot_prefix():
    old_files = os.listdir(image_destination)
    next_num = len(old_files)+1
    return "%04d_" % next_num

def test_runs():
    magnets.load_state(2.04, down=True)
    loc = [0.,0.,1.]
    B_vec = magnets.mag_field_of_i_at(i=2, xyz=loc)
    print(B_vec)
    print(magnets.total_field_at(xyz=loc))

def low_res():
    figure, subplots = pyplot.subplots(1,1)
    alpha = 2.45
    mono=True
    magnets.load_state(desired_alpha=alpha, down=mono)
    points_per_line=30
    size = 2
    # peak_B=.38
    Xs = numpy.linspace(-size, size, num=points_per_line, endpoint=True)
    Ys = numpy.linspace(-size, size, num=points_per_line, endpoint=True)
    z = 1.
    plotX = []
    plotY = []
    BX = []
    BY = []
    B_vals = []
    colors = []
    largest_B = 0
    for i in range(points_per_line):
        for j in range(points_per_line):
            x = Xs[i]
            y = Ys[j]
            loc = [x,y,z]
            B_hat,B_mag = magnets.total_field_at(xyz=loc)
            if B_mag> largest_B:
                largest_B = B_mag
            plotX.append(x)
            plotY.append(y)
            BX.append(B_hat[0])
            BY.append(B_hat[1])
            # blue = B_mag/peak_B
            # if blue > 1:
            #     blue=1

            B_vals.append(B_mag)
    for B_val in B_vals:
        blue = B_val/largest_B
        colors.append((1-blue, 1-blue, 1.))


    print(largest_B)



    add_circles(subplots)
    subplots.set_aspect('equal')
    subplots.set_title('%.03f,mono=%s' % (alpha, mono))

    subplots.quiver(
        plotX, plotY,
        BX, BY,
        color=colors
        # X=plotX, Y=plotY,
        # U=BX, V=BY
        )
    image_name = plot_prefix()+'simple_fields.png'
    destination = os.path.join(image_destination, image_name)
    pyplot.savefig(destination)

def add_circles(subplot):
    colors = []
    circle_patches = []
    for ind in range(7):
        tht0i = 2*pi*(ind-1)/6.
        if ind ==0:
            ci = [0.,0.]
        else:
            ci = [cos(tht0i), sin(tht0i)]
        circle_patches.append(
            mpatches.Circle(
                ci, radius=.12, alpha=.1, zorder=.5
                )
            )
    shapes = PatchCollection(circle_patches)
    subplot.add_collection(shapes)
    return

def void_check():
    magnets.load_state(2.04, down=True)
    epsi = .05
    for i in range(1,7):
        magnets.state[i]+=epsi
    print('total dipole', magnets.net_dipole_mag())
    # loc = [0.,0.,0.]
    # loc = [-0.1,0.,0.]
    # loc = [0.,-0.1,0.]
    spots = [
        [0.,0.,0.],
        # [0.1,0.,0.],
        # [0.,0.1,0.],
        # [-0.1,0.,0.],
        # [0.,-0.1,0.],
    ]
    for spot in spots:
        loc=spot
        tote_B = numpy.zeros(3)
        for j in range(1,7):
            B_vec = magnets.mag_field_of_i_at(i=j, xyz=loc)
            # print(B_vec)
            tote_B+= B_vec
        print(tote_B)
    # neighborhood of 4.5=B :/

if __name__ == "__main__":
    # test_runs()
    # print(plot_prefix())
    # low_res()
    void_check()