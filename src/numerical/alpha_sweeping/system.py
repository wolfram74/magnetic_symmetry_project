'''
state variables: 7 phi values (0, 2pi), 7 dot phi values (real)
parameters: alpha (>= 0), gamma (0, 1)
calculated geometry: r_ij matrix, theta_ij matrix
calculate forces: examine https://wolfram74.github.io/worked_problems/spring_20/week_2019_12_30/view.html
set up flow:
        calculate the geometry values
            delta x/y -> rij+theta_ij
            store those
        generate force equations, include drag term
            force equations will take in full state vector and return their contribution to change of state
                2*7 in, 2*7 out
    set initial phi values (0, and equi potentials)
    time evolution outline
        have state
            calculate kernel 1 by summing forces
                force equation could have as input self, state, and i,j, determine output from that, 1 function for all terms?
                    i==j -> drag, else torque
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
        self.alpha = 0.
        self.gamma = 0.
        self.state = numpy.zeros(7*2)
        self.kernels = numpy.zeros((6,7*2))
        self.kernel_step = 0
        self.delta_4 = numpy.zeros(2*7)
        self.delta_5 = numpy.zeros(2*7)
        self.step_size = 2**-3
        self.tolerance = 2**-20 #2**-10 ~= 1E-3
        self.elapsed=0
        # self.torques = [[0 for i in range(7)] for j in range(7)]
        self.calc_geometry()
        # self.calc_torques_eqs()
        self.set_init_state()

    def calc_geometry(self):
        for i in range(7):
            for j in range(7):
                self.set_geo_vals(i,j)
        return

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
        return

    # def calc_torques_eqs(self):
    #     for i in range(7):
    #         for j in range(7):
    #             self.set_torque_eq(i, j)
    #     return

    # def set_torque_eq(self, i,j):
    #     if i==j:
    #         self.torques[i][j] = self.gen_drag(i)
    #         return
    #     self.torques[i][j] = self.gen_torque(i,j)
    #     return

    # def gen_drag(self, i):
    #     def dragi(self):
    #         delta = numpy.zeros(2*7)
    #         delta[i+7] = -self.gamma*self.state[i+7]
    #         return delta
    #     return dragi

    # def gen_torque(self, i, j):
    #     deltas = numpy.zeros(2*7)
    #     strength = 1
    #     if 0 in (i,j):
    #         strength = self.alpha
    #     t1 = sin(self.state[i] - self.state[j])
    #     t2 = 3.*sin(self.state[i]+self.state[j]-2*self.theta_vals[i][j])
    #     r3 = self.r_vals[i][j]**(-3)
    #     deltas[i] =streng*r3*(t1+t2)
    #     return deltas

    def force_ij(self, state, i,j):
        deltas = numpy.zeros(2*7)
        if i==j:
            delta = numpy.zeros(2*7)
            delta[i+7] = -self.gamma*state[i+7]
            return delta
        strength = 1
        if 0 in (i,j):
            strength = self.alpha
        t1 = sin(state[i] - state[j])
        t2 = 3.*sin(state[i]+state[j]-2*self.theta_vals[i][j])
        r3 = self.r_vals[i][j]**(-3)
        deltas[i] =strength*r3*(t1+t2)
        return deltas

    def set_init_state(self):
        for i in range(1,7):
            self.state[i]=pi/2+2*pi*(i-1)/6
        return


    def total_force(self):
        tote_delta = numpy.zeros(2*7)
        temp_state = self.temp_state()
        for i in range(7):
            for j in range(7):
                tote_delta+= self.force_ij(temp_state, i, j)
        return tote_delta

    def temp_state(self):
        self.kernel_step += 1
        return self.state

    def rk45_step(self):
        for k in range(6)
            self.kernels[k] = self.step_size*self.total_force()
