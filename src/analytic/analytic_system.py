import sympy
from sympy.matrices import Matrix, zeros
cos = sympy.cos
sin = sympy.sin
atan = sympy.atan2
pi = sympy.pi
import numpy

class System():
    def __init__(self,):
        self.mu_hats = self.init_phi_vals()
        self.positions = self.init_x_y_vals()
        self.limit_mu_hats = self.high_alpha_limits()
        self.velocities = self.init_velo_vals()
        self.rad_matrix = self.make_r_vals()
        self.orientation_matrix = self.make_theta_vals()
        self.magnet_strengths = self.init_strengths()

    def init_phi_vals(self):
        return sympy.symbols('phi0 phi1 phi2 phi3 phi4 phi5 phi6', real=True)

    def init_velo_vals(self):
        return sympy.symbols('v0 v1 v2 v3 v4 v5 v6', real=True)

    def init_x_y_vals(self):
        outs = [(0.,0.), ]
        for i in range(6):
            ang = 2*pi*i/6
            xi = cos(ang)
            yi = sin(ang)
            outs.append( (xi, yi))
        return outs

    def make_r_vals(self):
        rads = []
        for i in range(7):
            row = []
            for j in range(7):
                if i==j:
                    row.append(0)
                    continue
                if 0 in [i,j]:
                    row.append(1)
                    continue
                values = (1,1,sympy.sqrt(3),2,sympy.sqrt(3),1)
                row.append( values[(i-j)%6])
            rads.append(row)
        return rads

    def make_theta_vals(self):
        thetas = []
        for i in range(7):
            row = []
            for j in range(7):
                sep_vec = self.sep_vec_ij(i, j)
                row.append(atan(sep_vec[1],sep_vec[0]))
            thetas.append(row)
        return thetas

    def sep_vec_ij(self,i,j):
        del_x = self.positions[j][0]-self.positions[i][0]
        del_y = self.positions[j][1]-self.positions[i][1]
        return (del_x,del_y)


    def init_strengths(self):
        a = sympy.symbols('alpha',positive=True)
        return [a,1,1,1,1,1,1]

    def b_hat_at_xy(self, xy_vec):
        #assumes B field originates from magnet 0
        theta_xy = atan(xy_vec[1],xy_vec[0])
        b_dir_x = 2*cos(theta_xy)**2 - sin(theta_xy)**2
        b_dir_y = 3*cos(theta_xy)* sin(theta_xy)
        magnitude = sympy.sqrt(b_dir_x**2+b_dir_y**2)
        b_hat = (b_dir_x/magnitude, b_dir_y/magnitude)
        return b_hat

    def high_alpha_limits(self):
        lim_phis = [0]
        for i in range(1,7):
            b_hat = self.b_hat_at_xy(self.positions[i])
            b_theta = atan(b_hat[1], b_hat[0])
            lim_phis.append(b_theta)
        return lim_phis


    def PE_ij(self, i,j):
        if i==j:
            return 0.
        strength = self.magnet_strengths[i]*self.magnet_strengths[j]
        t1 = cos(self.mu_hats[i] - self.mu_hats[j])
        t2 = 3.*cos(self.mu_hats[i]+self.mu_hats[j]-2*self.orientation_matrix[i][j])
        r3 = self.rad_matrix[i][j]**(-3)
        return -strength*r3*(t1+t2)/2.

    def PE_i(self,i):
        expression = 0
        for j in range(7):
            expression += self.PE_ij(i, j)
        return expression

    def PE_total(self):
        output = 0
        for i in range(7):
            output += self.PE_i(i)
        return output

    def acceleration_i(self, i):
        PE_term = self.PE_i(i)
        acceleration_term = -sympy.diff(PE_term, self.mu_hats[i])
        return acceleration_term

    def linearization(self, values):
        # high alpha limit
        array = []
        subs_array = self.make_subs_array(values)
        for i in range(7):
            row = []
            force_i = self.acceleration_i(i)
            for j in range(7):
                M_ij = sympy.diff(force_i, self.mu_hats[j])
                M_ij = M_ij.subs(subs_array)
                M_ij /= self.magnet_strengths[0]
                M_ij = sympy.limit(M_ij, self.magnet_strengths[0], sympy.oo)
                M_ij = M_ij.subs(subs_array).evalf()
                row.append(M_ij)
            array.append(row)
        matrix = sympy.Matrix(array)
        return matrix


    def make_subs_array(self, angle_vals):
        swap_vals = []
        for i in range(7):
            swap_vals.append((self.mu_hats[i], angle_vals[i]))
        return swap_vals

    def eval_func_specifically(self, func_name, angle_vals, i=None,j=None):
        concrete_vals = self.make_subs_array(angle_vals)
        if j!=None:
            expression = getattr(self, func_name)(i, j)
            return expression.subs(concrete_vals)
        if i!=None:
            expression = getattr(self, func_name)(i)
            return expression.subs(concrete_vals)
        expression = getattr(self, func_name)()
        return expression.subs(concrete_vals)

