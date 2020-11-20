import analytic_system
import sympy
import copy

magnets = analytic_system.System()

def simple_checks():
    sympy.pprint(magnets.mu_hats)
    sympy.pprint(magnets.velocities)
    sympy.pprint(magnets.positions)
    print('sep vec from 1 to 3')
    sympy.pprint(magnets.sep_vec_ij(1,3))
    [sympy.pprint(row) for row in magnets.rad_matrix]
    print('angles for sep vector from i to j')
    [sympy.pprint(row) for row in magnets.orientation_matrix]

def b_field_orientation():
    # for i in range(1,7):
    #     print(i)
    #     sympy.pprint(magnets.b_hat_at_xy(magnets.positions[i]))
    sympy.pprint(magnets.limit_mu_hats)
    print('sanity check pairs sum to 0')
    pairs = [(1,4),(2,6),(3,5),]
    for pair in pairs:
        phi_i = magnets.limit_mu_hats[pair[0]]
        phi_j = magnets.limit_mu_hats[pair[1]]
        sympy.pprint(phi_j+phi_i)

def term_calculation():
    sympy.pprint(magnets.PE_ij(0,1))
    sympy.pprint(magnets.PE_ij(2,5))
    # sympy.pprint(
    #     magnets.eval_func_specifically(
    #         'PE_ij', magnets.limit_mu_hats, 2,5
    #     )
    # )
    # total_PE_at_limit = magnets.eval_func_specifically(
    #         'PE_total', magnets.limit_mu_hats
    #     )

    # sympy.pprint(total_PE_at_limit)
    # sympy.pprint(total_PE_at_limit.expand().simplify())
    # sympy.pprint(total_PE_at_limit.expand().simplify().evalf())
    magnets.magnet_strengths[0] = 1
    total_PE_at_flat = magnets.eval_func_specifically(
            'PE_total', [0,0,0,0,0,0,0]
        )
    sympy.pprint(total_PE_at_flat.evalf())

def forces_in_alpha():
    F_2 = magnets.acceleration_i(2)
    sympy.pprint(F_2)
    off_bal = copy.copy(magnets.limit_mu_hats)
    off_bal[2]+=.01
    for k in range(7):
        F_2_equilibrium = magnets.eval_func_specifically(
            'acceleration_i', magnets.limit_mu_hats, i=k #no alpha term
            # 'acceleration_i', off_bal, i=k #possesses alpha term
            )
        # sympy.pprint(F_2_equilibrium)
        sympy.pprint(F_2_equilibrium.evalf())

def linearization():
    linear_mat = magnets.linearization(magnets.limit_mu_hats)
    x_var = sympy.symbols('x')
    sympy.pprint(linear_mat)
    print(sympy.latex(linear_mat))
    # characteristic = linear_mat.charpoly(x=x_var)
    # char_expr = characteristic.as_expr().collect(x_var)
    # for term in char_expr.args:
    #     print(term)
    #     sympy.pprint(term)
    # sympy.pprint(characteristic)
    eigen_vals = linear_mat.eigenvals()
    for key in eigen_vals.keys():
        # print(key, eigen_vals[key])
        # sympy.pprint(key.evalf())
        degen_eval_pair = (eigen_vals[key],float(sympy.re(key.evalf())))
        print("%d, %.04f" % degen_eval_pair)


def charpoly_play():
    x_var, a = sympy.symbols('x a')
    test = sympy.Matrix([
        [a,2,0],
        [2,2,0],
        [0,0,3],
        ])
    characteristic = test.charpoly(x=x_var)
    sympy.pprint(characteristic)
    char_expr = characteristic.as_expr().collect(x_var)
    sympy.pprint(char_expr)
    sympy.pprint(char_expr.simplify())
    sympy.pprint(sympy.limit(char_expr, a, sympy.oo))
    # sympy.pprint(char_expr.simplify().args)
    # sympy.pprint(char_expr.simplify().expand())
    # sympy.pprint(char_expr.simplify().expand().args)
    # eigen_vals = test.eigenvals()
    # for key in eigen_vals.keys():
    #     print(key, eigen_vals[key])
    #     sympy.pprint(key)
    #     sympy.pprint(eigen_vals[key])


if __name__=='__main__':
    # simple_checks()
    # b_field_orientation()
    # term_calculation()
    # forces_in_alpha()
    linearization()
    # charpoly_play()
'''
1, -0.1815 (mode 1)
3, -1.3229 (modes 2, 3, 6)
1, -1.9126 (mode 4)
1, -2.0000 (mode 5)
1, -10.5203 (mode 7)
(modes 2, 3, and 6 have same e-val)
3 and 6 have qualitative similaries,
2 odd one out, more like 4
1,5 and 7 have mostly synchronized movement
OMEGAS = [0.1815, 1.3229, 1.3229, 2.0000, 1.9126, 1.3229, 10.5203]
'''
