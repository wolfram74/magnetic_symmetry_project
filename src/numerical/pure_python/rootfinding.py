from scipy import optimize
import scipy
import symbolic_manipulations
import generate_equation
import numpy as np
import random


def find_solutions(positions, filename='Collage.png', iters=100):
    # generate important functions

    torque_function, potential_function, hessian_function = symbolic_manipulations.derive_functions(generate_equation.get_equations(positions), len(positions))


    '''
    Determine nature of a solution using the 2nd derivative (Hessian Matrix)

    '''
    def determine_nature(x):
        hessian = hessian_function(x)
        eigen_values = np.linalg.eig(hessian)[0]
        positive = np.sum([eigen_values > 0])
        negative = np.sum([eigen_values < 0])
        print(eigen_values)
        print(positive, negative)
        if (positive) > 0 and (negative) == 0:
            return 'local minimum'

        if (positive) ==0  and (negative) > 0:
            return 'local maximum'

        if (positive) > 0  and (negative) > 0:
            return 'saddle point'

        return 'indeterminate' # highly unlikely to happen


    # test near equality of two floats
    def has_seen(seen_u, u):
        for x in seen_u:
            if abs(x-u) < 0.00001:
                return True
            
        return False

    # seen_u: array of the potential energies of the solutions that have been found
    # seen: array of the corresponding rotations for each solution
    seen_u = []
    seen = []

    for i in range(1, iters):
        if not i % 1000:
            print('   '+str(i))
        output = scipy.optimize.minimize(torque_function, [random.uniform(0, 2*np.pi) for i in range(len(positions))])
        sol = output['x']
        #print(sol)
        u = potential_function(sol)
        print(u)
        if not has_seen(seen_u, u):    # don't double count identical or degenerate solutions
            tau = torque_function(sol)
            seen_u.append(u)
            seen.append(sol)
            print('candidate solution no.'+str(len(seen))+' found on iter.'+str(int(i)))
            print('     u='+str(float(u)) + '; tau^2=' + str(float(tau)))
            
            print('    ' + str(sol))

    sorted_u = sorted(seen_u)
    indeces = [seen_u.index(x) for x in sorted_u]
    sorted_sols = [seen[i] for i in indeces]




    '''
    now we draw a collage of the solutions we've found:

    '''
    import graphics

    torque_cutoff = 0.0000000001
    # ^^ we are quite confident that geniune solutions will have converged their sum of torques-squared to within 10^-10 of zero 
    candidate_sols = sorted_sols
    n=1
    cb = graphics.collage_builder()
    for sol in candidate_sols:
        tau = torque_function(sol)
        if tau < torque_cutoff:
            print('solution no.:'+str(n)+':'+str(sol))
            n = n + 1
            rotations = [np.array([np.cos(theta), np.sin(theta)]) for theta in sol]
            cb.add_solution(positions,rotations, potential_function(sol), tau, determine_nature(sol))
    cb.create_collage(filename)





