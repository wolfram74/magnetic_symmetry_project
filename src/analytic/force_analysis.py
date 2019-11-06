import sympy
from generate_U import full_potential, phi_gen
import cartesian_U_gen as c_U_gen
import pickle
import random
pi = sympy.pi
cos = sympy.cos
sin = sympy.sin



def parse_U_bin():
    Ubin = open('U_expr.pkl', 'rb')
    return pickle.load(Ubin)

def gen_phi_vals(symbols):
    # soln4 =(-0.523598754451836, -0.4168911504489822, 5.75958670639739, 5.65287959564158, 5.86629439213736, 5.75958670690163, 5.65287934383298)
    soln46 =(1.047185375472269, 1.280586331557914, 2.220244342353899, 3.482771775423115, 4.894816587679744, 6.15733787379985, 7.09698436863931)
    outputs = []
    for i in range(len(symbols)):
        outputs.append((symbols[i], soln46[i]))
    print(outputs)
    return outputs

def force_equations():
    phis = phi_gen()
    gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi

    total_U = parse_U_bin()
    # sympy.pprint(total_U)
    # total_U = sympy.simplify(total_U)
    # sympy.pprint(total_U)
    phi_dots = []
    phi_vals = gen_phi_vals(phis)
    for phi in phis:
        sympy.pprint(phi)
        phi_dots.append(-sympy.diff(total_U, phi))
        # sympy.pprint(phi_dots[-1])
        sympy.pprint(phi_dots[-1].subs(phi_vals).subs(gamma,1).evalf())

def bulk_evals():
    phis = phi_gen()
    gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi
    total_U = parse_U_bin()
    sympy.pprint(total_U*72)
    # for j in range(4):
    #     tht = random.random()
    #     print(tht)
    #     for i in range(1):
    #         # phi_vals = [(phi, random.random()) for phi in phis]
    #         phi_vals = [(phi, tht+i*pi/6.) for phi in phis]
    #         # print(phi_vals)
    #         sympy.pprint(total_U.subs(phi_vals).evalf())

def cart_force_equations():
    all_terms, coords = c_U_gen.full_potential()
    total_U = 0
    # output = open('U_expr.pkl', 'wb')
    for row in all_terms:
        for term in row:
            total_U+= term/2
    equations = []
    for i in range(7):
        for dim in coords:
            # sympy.pprint(dim[i])
            force = total_U.diff(dim[i]).simplify()
            # sympy.pprint(force)
            equations.append(force)
    # lam_vals = sympy.solve(equations)
    # sympy.pprint(lam_vals)
    # just_xy_forces = []
    # for i in range(len(equations)):
    #     if i%3==2:
    #         continue
    #     just_xy_forces.append(equations[i])
    # x_y_solutions = sympy.solve(just_xy_forces)
    # sympy.pprint(x_y_solutions)
    sympy.pprint(equations[:3])
    lam0_solnx = sympy.solve(equations[0], coords[2][0])
    lam0_solny = sympy.solve(equations[1], coords[2][0])
    sympy.pprint(lam0_solnx)
    sympy.pprint(lam0_solnx)
    sympy.pprint(lam0_solny)

# force_equations()

if __name__=='__main__':
    # bulk_evals()
    cart_force_equations()
