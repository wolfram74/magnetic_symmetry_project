import sympy
from generate_U import full_potential, phi_gen
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

# force_equations()

bulk_evals()
