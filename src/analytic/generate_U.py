'''
based on the work here
https://wolfram74.github.io/worked_problems/fall_19/week_2019_10_14/view.html
assumes l = 1
'''

import sympy
import pickle
pi = sympy.pi
cos = sympy.cos
sin = sympy.sin

def r_ij(i,j):
    if i == 0 or j==0:
        return 1
    values = (1,1,sympy.sqrt(3),2,sympy.sqrt(3),1)
    return values[(i-j)%6]

def theta_ij(i,j):
    if i ==0:
        return (j-1)*2*pi/6
    if j ==0:
        return (i-1)*2*pi/6+pi
    return ((1+2*i+(j-i)%6)%12)*pi/6

def phi_gen():
    return sympy.symbols('phi0 phi1 phi2 phi3 phi4 phi5 phi6')

def full_potential():
    phis = phi_gen()
    gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi
    sympy.pprint(phis)
    all_terms = []
    for i in range(7):
        row = []
        for j in range(7):
            if i==j:
                row.append(0)
                continue
            row.append(gamma*(
                cos(phis[i]-phis[j])
                +3*cos(
                    phis[i]+phis[j]-2*theta_ij(i,j)
                    )
                )/(r_ij(i,j)**3))
            # sympy.pprint([i,j, r_ij(i,j)**3])
        all_terms.append(row)
    return (all_terms, phis)

def mathjax_out(expr):
    return( sympy.latex(expr).replace('\\', '\\\\') )

def mathjax_table(all_terms):
    for ind in range(len(all_terms)):
        row = 'U_{%d,j}'%ind
        for ele_ind in range(len(all_terms[ind])):
            row+=", " + mathjax_out(all_terms[ind][ele_ind])
            if ele_ind%2==0:
                row+='\n \\\\\\\\ \\\\qquad \n'
        print(row)
        print('\\\\\\\\')

def for_mathjax():
    mathjax_table(full_potential()[0])

def pickle_potential():
    U_terms, phis = full_potential()
    gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi

    total_U = 0
    output = open('U_expr.pkl', 'wb')
    for row in U_terms:
        for term in row:
            total_U-= term/2
    sympy.pprint(total_U)
    total_U = sympy.collect(total_U, gamma)
    # sympy.pprint(total_U)
    # total_U = sympy.simplify(total_U)
    # sympy.pprint(total_U)
    # phi_vals = [(phi, 0.) for phi in phis]
    # U_tote = total_U.subs(phi_vals).evalf()
    # sympy.pprint(U_tote)
    # sympy.pprint(U_tote/7)

    pickle.dump(total_U, output)
    output.close()

if __name__=='__main__':
    # full_potential()
    pickle_potential()

'''
r + '^-5*(
    -3*(
        '+R_hat[0] +'*'+A[0]+'+' +R_hat[1] +'*'+ A[1] +'
        )*(
        '+R_hat[0] +'*'+B[0]+'+' +R_hat[1] +'*'+ B[1] +'
        )
) +'+ r + '^-3*('+ A[0]+'*'+B[0]+'+'+A[1]+'*'+B[1]+')'
question: is R_hat[0]**2+R_hat[1]**2=1 or =|R|**2

'''

