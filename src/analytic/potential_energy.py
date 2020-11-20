import sympy

sin = sympy.sin

x, y = sympy.symbols('x y', real=True)

expr = sin(x-y)+3*sin(x+y)
expr2 = 2sin(x-y)+1*sin(x+y)

sympy.pprint(sympy.solve([expr, expr2], x))
