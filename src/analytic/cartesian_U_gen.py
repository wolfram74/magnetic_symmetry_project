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
    return ((1+2*i+(j-i)%6)%12)*pi/6

def U_ij(i,j,xs, ys):
    rneg3 = r_ij(i,j)**-3
    tht = theta_ij(i,j)
    co_tht = cos(tht)
    sn_tht = sin(tht)
    return (
        xs[i]*xs[j]
        +ys[i]*ys[j]
        +3*(
            xs[i]*xs[j]
            -ys[i]*ys[j]
            )*co_tht
        +3*(
            xs[i]*ys[j]
            +ys[i]*xs[j]
            )*sn_tht
        )

def symb_gen(char):
    string = ' '.join(['%s%d' % (char, i) for i in range(7)])
    return sympy.symbols(string, real=True)


def full_potential():
    x_s = symb_gen('x')
    y_s = symb_gen('y')
    lam_s = symb_gen('lambda')
    # gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi
    # sympy.pprint(phis)
    all_terms = []
    for i in range(7):
        row = []
        for j in range(7):
            if i==j:
                row.append(0)
                continue
            row.append(U_ij(i,j,x_s,y_s))
        row.append(lam_s[i]*(x_s[i]**2+y_s[i]**2-1))
        all_terms.append(row)
    return (all_terms, (x_s, y_s, lam_s))

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
    U_terms, xs_ys_lams = full_potential()
    # gamma = sympy.symbols('gamma') #mu0*mu1*mu2/8pi

    total_U = 0
    # output = open('U_expr.pkl', 'wb')
    for row in U_terms:
        for term in row:
            total_U+= term/2
    sympy.pprint(total_U)
    total_U = sympy.simplify(total_U)
    sympy.pprint(total_U)
    # pickle.dump(total_U, output)
    # output.close()

if __name__=='__main__':
    # full_potential()
    pickle_potential()
