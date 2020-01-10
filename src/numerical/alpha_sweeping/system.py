'''
state variables: 7 phi values (0, 2pi), 7 dot phi values (real)
parameters: alpha (>= 0), gamma (0, `)
calculated geometry: r_ij matrix, theta_ij matrix
calculate forces: examine https://wolfram74.github.io/worked_problems/spring_20/week_2019_12_30/view.html
set up flow:
    calculate the geometry values
        delta x/y -> rij+theta_ij
        store those
    generate force equations, include drag term
    set initial phi values (0, and equi potentials)
'''

import numpy
cos = numpy.cos
sin = numpy.sin
atan = numpy.arctan2
pi = numpy.pi

class System():
    def __init__(self):
        self.r_vals = numpy.zeros((7,7))
        self.theta_vals = numpy.zeros((7,7))
        self.state = numpy.zeros(7*2)
        self.torques = [[0 for i in range(7)] for j in range(7)]
        self.calc_geometry()

    def calc_geometry(self):
        for i in range(7):
            for j in range(7):
                self.set_geo_vals(i,j)

    def set_geo_vals(self, i ,j):
        if i==j:
            return
        xi = cos(2*pi*(i-1)/6)
        yi = sin(2*pi*(i-1)/6)
        xj = cos(2*pi*(j-1)/6)
        yj = sin(2*pi*(j-1)/6)
        if i == 0:
            xi, yi = 0,0
        if j == 0:
            xj, yj = 0,0
        del_x = xi-xj
        del_y = yi-yj
        r = (del_x**2+del_y**2)**.5
        tht = atan(del_y, del_x)
        tht = atan(del_y, del_x)
        self.r_vals[i][j]=r
        self.theta_vals[i][j]=tht
