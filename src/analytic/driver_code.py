import analytic_system
import sympy

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
    sympy.pprint(
        magnets.eval_func_specifically(
            'PE_ij', magnets.limit_mu_hats, 2,5
        )
    )
    total_PE_at_limit = magnets.eval_func_specifically(
            'PE_total', magnets.limit_mu_hats
        )

    sympy.pprint(total_PE_at_limit)
    sympy.pprint(total_PE_at_limit.expand().simplify())
    sympy.pprint(total_PE_at_limit.expand().simplify().evalf())


if __name__=='__main__':
    # simple_checks()
    # b_field_orientation()
    term_calculation()
