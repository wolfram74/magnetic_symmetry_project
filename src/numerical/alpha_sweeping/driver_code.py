import system

magnets = system.System()
# print(magnets.r_vals)
# print(magnets.theta_vals)
print(magnets.state)
print(magnets.kernels)
# print(magnets.total_force())
magnets.rk45_step()
print(magnets.delta_4)
print(magnets.delta_5)
# print(magnets.forces[0][0])
