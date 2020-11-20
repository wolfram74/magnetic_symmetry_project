from rootfinding import find_solutions
import numpy as np

if __name__ == '__main__':

	# setup the positions of a spheres in a triangle
	positions = [np.array([0,0])] + [np.array([np.cos(i*np.pi/3), np.sin(i*np.pi/3)]) for i in range(2)]
	find_solutions(positions, filename='Collage_triangle.png', iters=10)

	# setup the positions of a spheres in a rhombus
	positions = [np.array([0,0])] + [np.array([np.cos(i*np.pi/3), np.sin(i*np.pi/3)]) for i in range(3)]
	find_solutions(positions, filename='Collage_rhombus.png', iters=20)

	# setup the positions of a spheres in a square
	positions = [np.array([np.cos(i*np.pi/2), np.sin(i*np.pi/2)])/np.sqrt(2) for i in range(4)]
	find_solutions(positions, filename='Collage_square.png', iters=30)

	# setup the positions of a spheres in a filled hexagon
	positions = [np.array([0,0])] + [np.array([np.cos(i*np.pi/3), np.sin(i*np.pi/3)]) for i in range(6)]
	find_solutions(positions, filename='Collage_filled_hexagon.png', iters=1000)
