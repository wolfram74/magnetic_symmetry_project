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
        self.L_mat = numpy.zeros((7,7))
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
        self.state*=0
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

    def load_state(self, desired_alpha, down=False):
        '''
        checks if stored_states.txt exists
        lines are formatted alpha, phi_0, ..., phi_6
        looks for largest alpha lower than desired
        sets alpha and phi's to those values
        alternative values available for .85 through 2.5
        '''
        if float(desired_alpha)==0.:
            self.set_init_state()
            return
        try:
            if down and desired_alpha> .85 and desired_alpha <2.5:
                saves = open('stored_states_down 2.txt', 'r')
            else:
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
        self.shift_alpha_and_stablize(0.)
        return

    def shift_alpha_and_stablize(self, del_alpha):
        T0 = self.elapsed
        T_damp = 2.
        running = True
        old_alpha = self.alpha
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
        if del_T>20.:
            print("going from %f to %f took %f units" % (old_alpha, self.alpha, del_T))
        # print(del_T,)
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
        return (forward[ind+7]-backward[ind+7])/step

    def extended_diff(self, ind, j, step):
        '''
        from rubin's comp physics, eq:5.11
        https://en.wikipedia.org/wiki/Richardson_extrapolation
        '''
        half = self.central_diff(ind, j, step/2.)
        full = self.central_diff(ind, j, step)
        return (4.*half-full)/3.

    def calc_L_mat(self, step):
        for i in range(7):
            for j in range(7):
                self.L_mat[i][j] = self.extended_diff(i, j, step)

    def spectrum_finder(self):
        if self.total_delta_sqr()>10**-8:
            print('%f calculated while delta^2 at %fE-6' %(self.alpha, self.total_delta_sqr()*10**6))
        self.calc_L_mat(.01)
        eigs = numpy.linalg.eig(self.L_mat)
        output = []
        for ind in range(7):
            output.append(
                [-eigs[0][ind], eigs[1][:, ind]]
                )
        # for ind in range(7):
        #     if output[ind][1][2]<0:
        #         flipped = numpy.array(output[ind][1])
        #         output[ind][1]= tuple(-1*flipped)
        return output

    def labeled_spectra(self):
        spectra = self.spectrum_finder()
        naive_sort = sorted(spectra, key=lambda el : el[0])
        #spectra is list of w^2, vector pairs
        labeled_vectors = {1:[],2:[],3:[],4:[],5:[],6:[],7:[]}
        for index in range(7):
            mode = naive_sort[index]
            val = mode[0]
            vec = mode[1]
            p14 = vec[1]/vec[4]
            p23 = vec[2]/vec[3]
            p56 = vec[5]/vec[6]
            if index >2:
                mode_id = index+1
            elif p14>0:
                mode_id = 1

            else:
                mode_id = self.mode_checker235(vec)
            new_vec = self.consistent_direction(mode_id, vec)
            labeled_vectors[mode_id].append(val)
            labeled_vectors[mode_id].append(new_vec)
        return labeled_vectors

    def mode_checker235(self, vector):
        #mode 2 has |p5|=|p6|>others
        #mode 3 has |p2|=|p3|>others
        #mode 5 has |p1|=|p4|>others
        abs_vals = numpy.abs(vector)
        if abs_vals[5]>abs_vals[2] and abs_vals[5]>abs_vals[1]:
            return 2
        if abs_vals[2]>abs_vals[1] and abs_vals[2]>abs_vals[5]:
            return 3
        return 5

    def mode_checker1467(self, vector):
        a01,a02,a05,a12,a15,a25 = self.sign_compare(vector)
        negatives = (a01,a02,a05,a12,a15,a25)
        pre1 = self.alpha <= 1.0
        pre2p01 = self.alpha <= 2.01
        if a01 and a02 and not a05 and not a12 and a15 and a25:
            return 1
        if pre1:
            if not any(negatives):
                return 7
        else:
            if all((a05, a15, a25)) and not any((a01, a02, a12)):
                return 7
        if pre1:
            if a01 and a05 and not a15:
                return 6
        elif not pre2p01:
            if all(negatives[1:5]) and not (a01 or a25):
                return 6
        else:
            if all(negatives[:3]) and not any(negatives[3:]):
                return 6
        # if a12:
        return 4


    def consistent_direction(self, label, vector):
        flipped = numpy.array(vector)
        flip = False
        new_vec = tuple(flipped)
        if label in [1,7]:
            if vector[0]<0:
                flip =True
        if label == 3:
            #3 -> el2 >0
            if vector[2]< 0:
                flip =True
        if label == 4:
            #4 -> el1 <0
            if vector[1] > 0:
                flip =True
        if label in [2,6]:
            #2 -> el5 >0
            #6 -> el5 >0
            if vector[5] >0:
                flip=True
        if label == 5:
            #5 -> el1 >0
            if vector[1] <0:
                flip=True
        if flip:
            new_vec = tuple(-1*flipped)
        # print(label, new_vec)
        return new_vec

    def sign_compare(self, vector):
        #a01, a02, a05, a12, a15, a25
        a01 = vector[0]*vector[1]<0
        a02 = vector[0]*vector[2]<0
        a05 = vector[0]*vector[5]<0
        a12 = vector[1]*vector[2]<0
        a15 = vector[1]*vector[5]<0
        a25 = vector[2]*vector[5]<0
        negatives = (a01,a02,a05,a12,a15,a25)
        return negatives



    def labeled_spectra_mono(self):
        spectra = self.spectrum_finder()
        naive_sort = sorted(spectra, key=lambda el : el[0])
        #spectra is list of w^2, vector pairs
        labeled_vectors = {1:[],2:[],3:[],4:[],5:[],6:[],7:[]}
        for index in range(7):
            mode = naive_sort[index]
            val = mode[0]
            vec = mode[1]
            p14 = vec[1]/vec[4]
            p23 = vec[2]/vec[3]
            p56 = vec[5]/vec[6]
            p25 =vec[2]/vec[5]
            sim2356 = self.nearly_equal(
                [vec[2],vec[3],vec[5],vec[6]]
                )
            sim014 = self.nearly_equal(
                [vec[0],vec[1],vec[4]]
                )
            #at high alpha  modes 4,5,6 (index 3,4,5) need differentiation
            if index == 6:
                mode_id = index+1
            if sim2356:
                #mode 1 or 5
                mode_id = index+1
            elif index in [0,1]:
                if p14 >0:
                    mode_id=1
                else:
                    mode_id=2
            else:
                #mode 4 and 5 (index 3 and 4)
                # if abs(vec[0])< 0.01 and abs(vec[1]) < .01:
                if p14>0:
                    mode_id=5
                else:
                    mode_id=4
            #     mode_id = self.mode_checker235(vec)
            new_vec = self.consistent_direction(mode_id, vec)
            labeled_vectors[mode_id].append(val)
            labeled_vectors[mode_id].append(new_vec)
        return labeled_vectors

    def nearly_equal(self, nums, tolerance=10**-6):
        #check all numbers are nearly equal
        if len(nums) in [0,1]:
            return True
        total = sum(nums)
        average = total/len(nums)
        for num in nums:
            if abs(num-average)>tolerance:
                return False
        return True
