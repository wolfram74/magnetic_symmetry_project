import sympy

cos = sympy.cos
sin = sympy.sin
atan = sympy.atan2
acos = sympy.acos
asin = sympy.asin
pi = sympy.pi
sqrt = sympy.sqrt


tht = sympy.symbols('theta', real=True)
# xhat = 2*cos(tht)**2-sin(tht)**2
# yhat = 3*cos(tht)*sin(tht)
xhat = (2*cos(tht)**2-sin(tht)**2).expand().simplify()
yhat = (3*cos(tht)*sin(tht)).expand().simplify()
rsqr = xhat**2+yhat**2
# rsqr = xhat.expand().simplify()**2+yhat.expand().simplify()**2

sympy.pprint(xhat)
sympy.pprint(yhat)
# sympy.pprint(rsqr)
tidy_r2 = rsqr.expand().simplify()
sympy.pprint(tidy_r2)

# sympy.pprint(sqrt(rsqr))
sympy.pprint(sqrt(tidy_r2).expand().simplify())
tidy_r1 = sqrt(tidy_r2)



phi_tan = atan(yhat,xhat)
phi_sin = asin(yhat/tidy_r1)
phi_cos = acos(xhat/tidy_r1)

term1 = cos(phi_cos)
term2 = 3*cos(phi_cos-2*tht)
term2_alt = 3*(cos(phi_cos)*cos(2*tht)+sin(phi_sin)*sin(2*tht))

print('term1')
sympy.pprint(term1)
sympy.pprint(term2)
sympy.pprint(term2.expand())
sympy.pprint(term2_alt)
print('term2')
sympy.pprint(term2_alt.expand().simplify())

