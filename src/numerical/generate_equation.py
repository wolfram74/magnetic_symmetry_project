import numpy as np
import math


def approx_equal(a, b, epsilon=0.00001):
	return abs(a-b) < epsilon

# modify to guess different things as required...
def guess_length(x):
	if approx_equal(x, math.sqrt(3)):
		return 'sqrt(3)'
	if approx_equal(x, 2):
		return '2'
	if approx_equal(x, 1):
		return '1'
	if approx_equal(x, 0):
		return '0'
	for num in range(-6, 7):
		for den in range(1, 7):
			if approx_equal(x, num/den):
				return str(num)+'/'+str(den)
			if approx_equal(x, num*math.sqrt(3)/den):
				return str(num)+'*sqrt(3)/' + str(den)
	raise ValueError('couldn\'t guess length: ', str(x))

def generate_tables(p):

	lengths_table = [[guess_length(np.linalg.norm(p[i]-p[j])) for j in range(len(p))] for i in range(len(p))]
	R_table  = [[[guess_length(p[j][0]-p[i][0]), guess_length(p[j][1]-p[i][1])] for j in range(len(p))] for i in range(len(p))]
	return lengths_table, R_table


# naming_mode: can either be "trig" or "cart"
# if substitute_root_3 is true, then every instance of "sqrt(3)" is replaced with "e"
def get_equations(positions, naming_mode='trig', substitute_root_3= False):
	parameter_names = []
	if naming_mode == "cart":
		for i in range(len(positions)):
			parameter_names.append(['X'+str(i), 'Y'+str(i)])
	elif naming_mode == 'trig':
		for i in range(len(positions)):
			parameter_names.append(['cos(a'+str(i)+')', 'sin(a'+str(i)+')'])
	else:
		raise ValueError(str(naming_mode) + ' is not a valid parameter')
	lengths_table, R_table = generate_tables(positions)
	terms = []
	for i in range(len(positions)):
		for j in range(i+1, len(positions)):
			terms.append(get_term(lengths_table[i][j], R_table[i][j], parameter_names[i], parameter_names[j]))
			
	string =  '+'.join(terms)
	string = string.replace('sqrt(3)^-3', 'sqrt(3)/9')
	string = string.replace('sqrt(3)^-5', 'sqrt(3)/27')
	if substitute_root_3:
		string = string.replace('sqrt(3)', 'e')
	return string

'''
apply U12 = - (1/r^3)*(m1.m2) +(3/r^5) * (m1.r)*(m2.r) 
'''
def get_term(r, R_hat, A, B):
	return r + '^-5*(-3*('+R_hat[0] +'*'+A[0]+'+' +R_hat[1] +'*'+ A[1] +')*('+R_hat[0] +'*'+B[0]+'+' +R_hat[1] +'*'+ B[1] +')) +'+ r + '^-3*('+ A[0]+'*'+B[0]+'+'+A[1]+'*'+B[1]+')'


if __name__ == '__main__':
	positions = [np.array([0,0])] + [np.array([np.cos(i*np.pi/3), np.sin(i*np.pi/3)]) for i in range(6)]
	#print(get_equations(positions, naming_mode='cart')) # <- for trying to solve it as a system of simultaneous equations
	print(get_equations(positions, naming_mode='trig'))
	input()
