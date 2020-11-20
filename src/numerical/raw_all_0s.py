import sympy
import numpy
import pickle
pi = numpy.pi
cos = numpy.cos
sin = numpy.sin


def r_ij(i,j):
    if i == 0 or j==0:
        return 1
    values = (1,1,numpy.sqrt(3),2,numpy.sqrt(3),1)
    return values[(i-j)%6]

def theta_ij(i,j):
    if i ==0:
        return (j-1)*2*pi/6.
    if j ==0:
        return (i-1)*2*pi/6.+pi
    return ((1+2*i+(j-i)%6)%12)*pi/6.

def U_ij(i,j):
    return -1*(
        cos(0+0)+3*cos(0+0-2*theta_ij(i,j))
        )/(r_ij(i,j)**3)

def total_U_at_0():
    sum_U = 0
    for i in range(7):
        print(i)
        for j in range(7):
            if i==j:
                continue
            if i==0 or j ==0:
                # print(i,j, U_ij(i,j), r_ij(i,j), theta_ij(i,j))
                pass
            print(i,j, U_ij(i,j))
            sum_U += U_ij(i,j)/2
    print(sum_U)
    print(sum_U/7.)

if __name__=='__main__':
    total_U_at_0()
