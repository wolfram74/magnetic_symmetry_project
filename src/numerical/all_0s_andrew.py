import numpy as np
from math import cos, sin, pi

positions = [np.array([0, 0])] + [np.array([cos(i*pi/3),sin(i*pi/3)]) for i in range(6)]

def u(i, j, positions):
    R_ij = positions[i]-positions[j]
    r = np.linalg.norm(R_ij)
    # print(r)
    M = np.array([1, 0])
    return (-3 * np.dot(R_ij,M) * np.dot(R_ij,M) / r**5
        + np.dot(M , M) / r**3)

total = 0
for i in range(0, 7):
    print(i)
    for j in range(0, 7):
        if i==j:
            continue
        print(i, j, u(i, j, positions))
        total+=u(i, j, positions)/2

print(total)

