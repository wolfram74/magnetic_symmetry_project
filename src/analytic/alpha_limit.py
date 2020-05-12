import sympy



def para_bhat():
    tht = sympy.symbols('theta', real=True)
    xhat = 2*sympy.cos(tht)**2 - sympy.sin(tht)**2
    yhat = 3*sympy.cos(tht)* sympy.sin(tht)
    sympy.pprint(xhat)
    xhat = xhat.simplify()
    yhat = yhat.simplify()
    print('simpler xhat, yhat')
    sympy.pprint([xhat, yhat])
    r2 = xhat**2+yhat**2
    sympy.pprint(r2)
    r2 = r2.simplify()
    print('simplified rhat')
    sympy.pprint(r2)
    # sympy.pprint(sympy.sqrt(r2).simplify())
    # r2g = (2+sympy.sqrt(3)*sympy.sin(tht))*(2-sympy.sqrt(3)*sympy.sin(tht))
    # sympy.pprint(r2 - r2g.expand())
    # arg = yhat/sympy.sqrt(r2)
    # sympy.pprint(sympy.asin(arg))
    # sympy.pprint(sympy.asin(arg.expand()))
    # sympy.pprint(sympy.asin(arg.expand().simplify()))
    # norming
    phix = sympy.acos(xhat/sympy.sqrt(r2))
    phiy = sympy.asin(yhat/sympy.sqrt(r2))
    #angle phi_of_theta
    phit = sympy.atan2(yhat, xhat)
    # sympy.pprint(phix.subs(tht, sympy.pi/3).evalf())
    # sympy.pprint(phiy.subs(tht, sympy.pi/3).evalf())
    print('phi of theta')
    # sympy.pprint(phit.expand().simplify())
    sympy.pprint(phit.subs(tht, sympy.pi/3).evalf())

    sympy.pprint((phit.subs(tht, sympy.pi/3)/sympy.pi).evalf())
    ratio = (phit.subs(tht, sympy.pi/3)/sympy.pi).evalf()
    # print('pi guess')
    # sympy.pprint(10*ratio**2)
    # sympy.pprint((10*ratio**2-sympy.pi).evalf())
    # sympy.pprint(phit.subs(tht, 2*sympy.pi/3).evalf())
    # sympy.pprint(phit.subs(tht, 3*sympy.pi/3))
    tht2rat = phit.subs(tht, sympy.pi/3)/sympy.pi
    print('phi_2 in units of pi')
    sympy.pprint(tht2rat)
    sympy.pprint(tht2rat.expand().simplify())
    sympy.pprint(tht2rat.expand().simplify().evalf())
    # sympy.pprint(10*(tht2rat**2).expand().simplify())
    # sympy.pprint(10*(tht2rat**2).expand().simplify().evalf())
    sympy.pprint(((10*tht2rat)**2).expand().simplify().evalf())
    sympy.pprint((
        phit.subs(tht, sympy.pi/3).evalf() - (.1*sympy.pi)**2
        ).evalf())

def limit_phis():
    tht = sympy.symbols('theta', real=True)
    xhat = 2*sympy.cos(tht)**2 - sympy.sin(tht)**2
    yhat = 3*sympy.cos(tht)* sympy.sin(tht)
    xhat = xhat.simplify()
    yhat = yhat.simplify()
    phit = sympy.atan2(yhat, xhat)
    return [0, 0,
        phit.subs(tht, sympy.pi/3),phit.subs(tht, 2*sympy.pi/3),
        0,
        -phit.subs(tht, 2*sympy.pi/3), -phit.subs(tht, sympy.pi/3)
        ]

def dipole_limit():
    alpha = sympy.symbols('alpha', real=True)
    phis = limit_phis()
    mu_x = alpha
    mu_y = 0
    for phi in phis[1:]:
        mu_x += sympy.cos(phi)
        mu_y += sympy.sin(phi)
    sympy.pprint([mu_x, mu_y])
    mag2 = mu_x**2 + mu_y**2
    sympy.pprint(mag2)
    sympy.pprint(mag2.expand())
    sympy.pprint(mag2.expand().simplify())
    sympy.pprint(mag2.expand().evalf())
    mag = sympy.sqrt(mag2.evalf())
    sympy.pprint(mag)
    sympy.pprint(mag.evalf())

if __name__ == '__main__':
    para_bhat()
    # sympy.pprint(limit_phis())
    # dipole_limit()

