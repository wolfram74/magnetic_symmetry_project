'''
based on the work here
https://wolfram74.github.io/worked_problems/fall_19/week_2019_10_14/view.html
assumes l = 1
'''
import sympy
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

def full_potential():
    phis = sympy.symbols('phi0 phi1 phi2 phi3 phi4 phi5 phi6')
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
                )/r_ij(i,j)**3)
        all_terms.append(row)
    mathjax_table(all_terms)
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


full_potential()
