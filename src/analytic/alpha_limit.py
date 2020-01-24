import sympy



def para_bhat():
    tht = sympy.symbols('theta', real=True)
    xhat = 2*sympy.cos(tht)**2 - sympy.sin(tht)**2
    yhat = 3*sympy.cos(tht)* sympy.sin(tht)
    sympy.pprint(xhat)
    sympy.pprint(xhat.simplify())
    sympy.pprint(yhat)
    sympy.pprint(yhat.simplify())


if __name__ == '__main__':
    para_bhat()
