import sympy
import numpy as np
import scipy
from sympy import sin, cos, sqrt
from sympy.parsing.sympy_parser import parse_expr
from sympy.utilities.iterables import flatten





def derive_functions(potential_expression_string, num_unknowns, var_set=None):

	if var_set == None:
		var_set = ['a'+str(i) for i in range(num_unknowns)]

	potential_expression = parse_expr(potential_expression_string) 

	# build the sum of squares of torques function, f
	taus = sum([potential_expression.diff(var)*potential_expression.diff(var) for var in var_set])

	#get hessian
	hessian_matrix = [[potential_expression.diff(var1).diff(var2) for var1 in var_set] for var2 in var_set]


	# convert expression into function that can be accepted by scipy.optimize
	_temp_func = sympy.lambdify(var_set, taus)
	def torque_function(x):
	    #return _temp_func(x[0], x[1], x[2] , x[3], x[4], x[5], x[6])
	    return _temp_func(*flatten(x))
	_temp_func2 = sympy.lambdify(var_set, potential_expression)
	def potential_function(x):
	    return _temp_func2(*flatten(x))

	_temp_func3 = sympy.lambdify(var_set, hessian_matrix)
	def hessian_function(x):
		return _temp_func3(*flatten(x))

	return torque_function, potential_function, hessian_function