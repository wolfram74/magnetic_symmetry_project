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
        # self.state[0]+=pi/120
        self.kernels = numpy.zeros((6,7*2))
        self.kernel_step = 0
        self.delta_4 = numpy.zeros(2*7)
        self.delta_5 = numpy.zeros(2*7)
        self.step_size = 2**-6
        self.tolerance = 2**-30 #2**-10 ~= 1E-3
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

    def force_ij(self, state, i,j):
        deltas = numpy.zeros(2*7)
        if i==j:
            delta = numpy.zeros(2*7)
            delta[i+7] = -self.gamma*state[i+7]
            delta[i] = state[i+7]
            return delta
        strength = self.dipole_strength(i)*self.dipole_strength(j)
        t1 = sin(state[i] - state[j])
        t2 = 3.*sin(state[i]+state[j]-2*self.theta_vals[i][j])
        r3 = self.r_vals[i][j]**(-3)
        # if i == 0:
        #     print(-strength*r3*(t1+t2)/2)
        deltas[i+7] =-strength*r3*(t1+t2)/2
        return deltas

    def PE_ij(self, i,j):
        if i==j:
            return 0.
        strength = self.dipole_strength(i)*self.dipole_strength(j)
        t1 = cos(self.state[i] - self.state[j])
        t2 = 3.*cos(self.state[i]+self.state[j]-2*self.theta_vals[i][j])
        r3 = self.r_vals[i][j]**(-3)
        return -strength*r3*(t1+t2)/2.

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
            # if i== 0:
            #     print('deep', tote_delta)
        return tote_delta

    def temp_state(self):
        if self.kernel_step==0:
            temp_vals = (
                self.state
                +self.kernels[0]*0.
                )
        if self.kernel_step==1:
            temp_vals = (
                self.state
                +self.kernels[0]/4.
                )
        if self.kernel_step==2:
            temp_vals = (
                self.state
                +3.*self.kernels[0]/32.
                +9.*self.kernels[1]/32.
                )
        if self.kernel_step==3:
            temp_vals = (
                self.state
                +1932.*self.kernels[0]/2197.
                -7200.*self.kernels[1]/2197.
                +7296.*self.kernels[2]/2197.
                )
        if self.kernel_step==4:
            temp_vals = (
                self.state
                +439.*self.kernels[0]/216.
                -8.*self.kernels[1]/1.
                +3680.*self.kernels[2]/513.
                -845.*self.kernels[3]/4104.
                )
        if self.kernel_step==5:
            temp_vals = (
                self.state
                -8.*self.kernels[0]/27.
                +2.*self.kernels[1]/1.
                -3544.*self.kernels[2]/2565.
                +1859.*self.kernels[3]/4104.
                -11.*self.kernels[4]/40.
                )
        return temp_vals

    def rk45_step(self):
        for k in range(6):
            self.kernels[k] = self.step_size*self.total_force()
            self.kernel_step += 1
        self.kernel_step = 0
        self.delta_4 = (
            25.*self.kernels[0]/216.
            +1408.*self.kernels[2]/2565.
            +2197.*self.kernels[3]/4104.
            -1.*self.kernels[4]/5.
            )
        self.delta_5 = (
            16.*self.kernels[0]/135.
            +6656.*self.kernels[2]/12825.
            +28561.*self.kernels[3]/56430.
            -9.*self.kernels[4]/50.
            +2.*self.kernels[5]/55.
            ) #O
        return

    def step_rescale(self):
        rescale = 10.**6
        numer = self.tolerance*self.step_size
        for ind in range(len(self.delta_4)):
            diff = abs(self.delta_4[ind]-self.delta_5[ind])
            if diff ==0:
                proposed = 1
                continue
            proposed = (
                numer/(2.*diff)
            )**.25
            if proposed < rescale:
                rescale = proposed
        return rescale

    def advance_in_time(self):
        self.rk45_step()
        adjust = self.step_rescale()
        # print(adjust, self.step_size, self.elapsed)
        if adjust < .9:
            self.step_size *= adjust
            self.advance_in_time()
        self.state+=self.delta_5
        self.elapsed += self.step_size
        if adjust > 1.25:
            adjust = 1.25
        self.step_size *= adjust
        return self.elapsed

    def total_KE(self):
        KE = sum((self.state**2)[7:])
        return KE

    def total_PE(self):
        PE = 0.
        for i in range(7):
            for j in range(7):
                PE += self.PE_ij(i,j)
        return PE

    def total_delta_sqr(self):
        changes = self.total_force()
        return sum( (changes**2) )


    def net_dipole_moment(self):
        moment = numpy.zeros(2)
        for i in range(7):
            moment+=self.mag_moment_i(i)
        return moment

    def net_dipole_mag(self):
        return numpy.linalg.norm(self.net_dipole_moment())

    def dipole_strength(self, i):
        if i== 0:
            return self.alpha
        return 1.

    def mag_moment_i(self, i):
        direction = numpy.zeros(2)
        direction[0] = cos(self.state[i])
        direction[1] = sin(self.state[i])
        return self.dipole_strength(i)*direction

    def lim_U(self):
        limit = -13.8440856265381*self.alpha
        return self.total_PE()/limit

    def lim_moment(self):
        limit = self.alpha + 1.24407105398155
        return self.net_dipole_mag()/limit

    def load_state(self, desired_alpha):
        '''
        checks if stored_states.txt exists
        lines are formatted alpha, phi_0, ..., phi_6
        looks for largest alpha lower than desired
        sets alpha and phi's to those values
        '''
        if desired_alpha==0:
            self.set_init_state()
            return
        try:
            saves = open('stored_states.txt', 'r')
            last_alpha = 0.
            last_state = None
            for line in saves:
                txt_list = line.split(',')
                txt_list.pop()
                num_list = [float(el) for el in txt_list]
                if num_list[0]>desired_alpha:
                    self.alpha = last_alpha
                    self.state = last_state
                    return
                last_alpha = num_list[0]
                # print(num_list[1:])
                last_state = numpy.concatenate((numpy.array(num_list[1:]),numpy.zeros(7)))
            saves.close()
        except:
            print('failed to load saved state, alpha unchanged. Please run sweeper.py')
            return
        return

    def shift_alpha_and_stablize(self, del_alpha):
        T0 = self.elapsed
        T_damp = 2.
        running = True
        self.alpha += del_alpha
        self.gamma = 0.
        while running:
            self.advance_in_time()
            del_T = self.elapsed-T0
            if del_T > T_damp:
                scale = (del_T-T_damp)/del_T
                self.gamma = 2.*(self.alpha+1)*scale
            force_crit = self.total_delta_sqr() < 10**-8
            time_crit = del_T > 5.
            if not force_crit or not time_crit:
                continue
            running = False
        print(del_T,)
        return


    def central_diff(self, ind, j, step):
        '''
        from rubin's comp physics, eq:5.7
        '''
        self.state[j]+=step/2.
        forward = self.total_force()
        self.state[j]-=step
        backward = self.total_force()
        self.state[j]+=step/2.
        return (forward[ind+7]-back[ind+7])/step

    def extended_diff(self, ind, j, step):
        '''
        from rubin's comp physics, eq:5.11
        '''
        half = self.central_diff(j, step/2.)
        full = self.central_diff(j, step)
        return (4.*half-full)/3.

