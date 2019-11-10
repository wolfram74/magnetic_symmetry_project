#!/usr/bin/env python
# coding: utf-8

# In[1]:


#define variables
a0, a1, a2, a3, a4, a5, a6 = var('a0 a1 a2 a3 a4 a5 a6')


# In[2]:


# the potential energy function, generated by generate_equation.py
pot(a0, a1, a2, a3, a4, a5, a6) = 1^-5*(-3*(1*cos(a0)+0*sin(a0))*(1*cos(a1)+0*sin(a1))) +1^-3*(cos(a0)*cos(a1)+sin(a0)*sin(a1))+1^-5*(-3*(1/2*cos(a0)+1*sqrt(3)/2*sin(a0))*(1/2*cos(a2)+1*sqrt(3)/2*sin(a2))) +1^-3*(cos(a0)*cos(a2)+sin(a0)*sin(a2))+1^-5*(-3*(-3/6*cos(a0)+1*sqrt(3)/2*sin(a0))*(-3/6*cos(a3)+1*sqrt(3)/2*sin(a3))) +1^-3*(cos(a0)*cos(a3)+sin(a0)*sin(a3))+1^-5*(-3*(-6/6*cos(a0)+0*sin(a0))*(-6/6*cos(a4)+0*sin(a4))) +1^-3*(cos(a0)*cos(a4)+sin(a0)*sin(a4))+1^-5*(-3*(-3/6*cos(a0)+-3*sqrt(3)/6*sin(a0))*(-3/6*cos(a5)+-3*sqrt(3)/6*sin(a5))) +1^-3*(cos(a0)*cos(a5)+sin(a0)*sin(a5))+1^-5*(-3*(1/2*cos(a0)+-3*sqrt(3)/6*sin(a0))*(1/2*cos(a6)+-3*sqrt(3)/6*sin(a6))) +1^-3*(cos(a0)*cos(a6)+sin(a0)*sin(a6))+1^-5*(-3*(-3/6*cos(a1)+1*sqrt(3)/2*sin(a1))*(-3/6*cos(a2)+1*sqrt(3)/2*sin(a2))) +1^-3*(cos(a1)*cos(a2)+sin(a1)*sin(a2))+sqrt(3)/27*(-3*(-6/4*cos(a1)+1*sqrt(3)/2*sin(a1))*(-6/4*cos(a3)+1*sqrt(3)/2*sin(a3))) +sqrt(3)/9*(cos(a1)*cos(a3)+sin(a1)*sin(a3))+2^-5*(-3*(-6/3*cos(a1)+0*sin(a1))*(-6/3*cos(a4)+0*sin(a4))) +2^-3*(cos(a1)*cos(a4)+sin(a1)*sin(a4))+sqrt(3)/27*(-3*(-6/4*cos(a1)+-3*sqrt(3)/6*sin(a1))*(-6/4*cos(a5)+-3*sqrt(3)/6*sin(a5))) +sqrt(3)/9*(cos(a1)*cos(a5)+sin(a1)*sin(a5))+1^-5*(-3*(-3/6*cos(a1)+-3*sqrt(3)/6*sin(a1))*(-3/6*cos(a6)+-3*sqrt(3)/6*sin(a6))) +1^-3*(cos(a1)*cos(a6)+sin(a1)*sin(a6))+1^-5*(-3*(-6/6*cos(a2)+0*sin(a2))*(-6/6*cos(a3)+0*sin(a3))) +1^-3*(cos(a2)*cos(a3)+sin(a2)*sin(a3))+sqrt(3)/27*(-3*(-6/4*cos(a2)+-3*sqrt(3)/6*sin(a2))*(-6/4*cos(a4)+-3*sqrt(3)/6*sin(a4))) +sqrt(3)/9*(cos(a2)*cos(a4)+sin(a2)*sin(a4))+2^-5*(-3*(-6/6*cos(a2)+-6*sqrt(3)/6*sin(a2))*(-6/6*cos(a5)+-6*sqrt(3)/6*sin(a5))) +2^-3*(cos(a2)*cos(a5)+sin(a2)*sin(a5))+sqrt(3)/27*(-3*(0*cos(a2)+-6*sqrt(3)/6*sin(a2))*(0*cos(a6)+-6*sqrt(3)/6*sin(a6))) +sqrt(3)/9*(cos(a2)*cos(a6)+sin(a2)*sin(a6))+1^-5*(-3*(-3/6*cos(a3)+-3*sqrt(3)/6*sin(a3))*(-3/6*cos(a4)+-3*sqrt(3)/6*sin(a4))) +1^-3*(cos(a3)*cos(a4)+sin(a3)*sin(a4))+sqrt(3)/27*(-3*(0*cos(a3)+-6*sqrt(3)/6*sin(a3))*(0*cos(a5)+-6*sqrt(3)/6*sin(a5))) +sqrt(3)/9*(cos(a3)*cos(a5)+sin(a3)*sin(a5))+2^-5*(-3*(1*cos(a3)+-6*sqrt(3)/6*sin(a3))*(1*cos(a6)+-6*sqrt(3)/6*sin(a6))) +2^-3*(cos(a3)*cos(a6)+sin(a3)*sin(a6))+1^-5*(-3*(1/2*cos(a4)+-3*sqrt(3)/6*sin(a4))*(1/2*cos(a5)+-3*sqrt(3)/6*sin(a5))) +1^-3*(cos(a4)*cos(a5)+sin(a4)*sin(a5))+sqrt(3)/27*(-3*(3/2*cos(a4)+-3*sqrt(3)/6*sin(a4))*(3/2*cos(a6)+-3*sqrt(3)/6*sin(a6))) +sqrt(3)/9*(cos(a4)*cos(a6)+sin(a4)*sin(a6))+1^-5*(-3*(1*cos(a5)+0*sin(a5))*(1*cos(a6)+0*sin(a6))) +1^-3*(cos(a5)*cos(a6)+sin(a5)*sin(a6))
pot = pot.expand()


# In[3]:


# build the sum of squares of torques function, f
f = sum(k*k for k in pot.diff())


# In[4]:


# matrix derivative - for second derivative test
H = pot.diff(2)


# In[5]:


def has_seen(seen_u, u):
    for x in seen_u:
        if abs(x-u) < 0.00001:
            return True
        
    return False

# seen_u: array of the potential energies of the solutions that have been found
# seen: array of the corresponding rotations for each solution
seen_u = []
seen = []

import random
import numpy as np
for i in range(1, 100):
    if not i % 1000:
        print('   '+str(i))
    sol = minimize(f, [random.uniform(0, np.pi/6), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi)])
    u = pot(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n()
    if not has_seen(seen_u, u):
        tau = f(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n()
        seen_u.append(u)
        seen.append(sol)
        print('candidate solution no.'+str(len(seen))+' found on iter.'+str(int(i)))
        print('     u='+str(float(u)) + '; tau^2=' + str(float(tau)))
        
        print('    ' + str(sol))
print(seen_u)
print(seen)


# In[6]:


print(len(seen_u))
sorted_u = sorted(seen_u)
print(sorted_u)
indeces = [seen_u.index(x) for x in sorted_u]
sorted_sols = [seen[i] for i in indeces]


# In[7]:


import matplotlib.pyplot as plt
def display_balls(positions, rotations, u, tau):
    n = len(positions)
    circles = [plt.Circle((positions[i][0], positions[i][1]), 0.5, color='gray', ec = 'k') for i in range(len(positions))]

    fig, ax = plt.subplots(figsize=(5, 5)) # note we must use plt.subplots, not plt.subplot

    ax.set_aspect('equal')
    ax.set_ylim(-5, 5)
    ax.set_xlim(-5, 5)

    for c in circles:
        ax.add_artist(c)

    X = np.array([float(p[0]) for p in positions])
    Y = np.array([float(p[1]) for p in positions])
    U = np.array([float(r[0]) for r in rotations])
    V = np.array([float(r[1]) for r in rotations])
    X -= U/2
    Y -= V/2
    ax.quiver(X, Y, U, V, units='xy', scale=1, color='r', zorder=10)
    
    dx, dy = 0,0
    for x in rotations:
        dx+=x[0]
        dy+=x[1]
    
    plt.text(int(-4), int(-4), 'U='+str(u)+'\nm='+str(sqrt(dx*dx+dy*dy))+'\ntau^2='+str(tau), color='k')
    plt.show()


# In[8]:


def interpret_eigen(x):
    p = 0
    n = 0
    for k in x:
        if k < 0:
            n+=1
        if k > 0:
            p+=1
    if n == 0:
        return 'local minimum'
    if p == 0:
        return 'local maximum'
    return 'saddle point'


# In[9]:


def d_hex(sol):
    
    positions = [[0, 0]]+[[np.cos(x*np.pi/3), np.sin(x*np.pi/3)] for x in range(6)]
    rotations = [[np.cos(x), np.sin(x)] for x in sol]
    u= (pot(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n())
    tau= float(f(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n())
    eigens = H(a0=sol[0], a1=sol[1], a2=sol[2], a3=sol[3], a4=sol[4], a5=sol[5], a6=sol[6]).n().eigenvalues()
    #print('eigens='+str(eigens))
    print('solution:'+ str(sol)+ '\nnature = '+interpret_eigen(eigens))
    display_balls(positions, rotations, u, tau)


# In[10]:



torque_cutoff = 0.0000000001
# ^^ we are quite confident that geniune solutions will have converged to within 10^-10 of zero 
candidate_sols = sorted_sols
n=1
for sol in candidate_sols:
    tau = f(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n()
    if tau < torque_cutoff:
        print('solution no.:'+str(n))
        n = n + 1
        d_hex(sol)




# In[11]:


# display solutions which were rejected so we can check them
for sol in candidate_sols:
    tau = f(sol[0], sol[1], sol[2], sol[3], sol[4], sol[5], sol[6]).n()
    if tau >= torque_cutoff:
        print('solution no.:'+str(n))
        n = n + int(1)
        d_hex(sol)


# In[12]:


# code for testing

def torque_array(Phi):
    a0, a1, a2, a3, a4, a5, a6 = Phi
    e= sqrt(3)
    torque_eqs = [1^-5*(3*(1*cos(a1)+0*sin(a1))*(cos(a0)*0-sin(a0)*1)) +1^-3*(cos(a1)*sin(a0)-sin(a1)*cos(a0))+1^-5*(3*(1/2*cos(a2)+1*e/2*sin(a2))*(cos(a0)*1*e/2-sin(a0)*1/2)) +1^-3*(cos(a2)*sin(a0)-sin(a2)*cos(a0))+1^-5*(3*(-3/6*cos(a3)+1*e/2*sin(a3))*(cos(a0)*1*e/2-sin(a0)*-3/6)) +1^-3*(cos(a3)*sin(a0)-sin(a3)*cos(a0))+1^-5*(3*(-6/6*cos(a4)+0*sin(a4))*(cos(a0)*0-sin(a0)*-6/6)) +1^-3*(cos(a4)*sin(a0)-sin(a4)*cos(a0))+1^-5*(3*(-3/6*cos(a5)+-3*e/6*sin(a5))*(cos(a0)*-3*e/6-sin(a0)*-3/6)) +1^-3*(cos(a5)*sin(a0)-sin(a5)*cos(a0))+1^-5*(3*(1/2*cos(a6)+-3*e/6*sin(a6))*(cos(a0)*-3*e/6-sin(a0)*1/2)) +1^-3*(cos(a6)*sin(a0)-sin(a6)*cos(a0)),    1^-5*(3*(-6/6*cos(a0)+0*sin(a0))*(cos(a1)*0-sin(a1)*-6/6)) +1^-3*(cos(a0)*sin(a1)-sin(a0)*cos(a1))+1^-5*(3*(-3/6*cos(a2)+1*e/2*sin(a2))*(cos(a1)*1*e/2-sin(a1)*-3/6)) +1^-3*(cos(a2)*sin(a1)-sin(a2)*cos(a1))+e/27*(3*(-6/4*cos(a3)+1*e/2*sin(a3))*(cos(a1)*1*e/2-sin(a1)*-6/4)) +e/9*(cos(a3)*sin(a1)-sin(a3)*cos(a1))+2^-5*(3*(-6/3*cos(a4)+0*sin(a4))*(cos(a1)*0-sin(a1)*-6/3)) +2^-3*(cos(a4)*sin(a1)-sin(a4)*cos(a1))+e/27*(3*(-6/4*cos(a5)+-3*e/6*sin(a5))*(cos(a1)*-3*e/6-sin(a1)*-6/4)) +e/9*(cos(a5)*sin(a1)-sin(a5)*cos(a1))+1^-5*(3*(-3/6*cos(a6)+-3*e/6*sin(a6))*(cos(a1)*-3*e/6-sin(a1)*-3/6)) +1^-3*(cos(a6)*sin(a1)-sin(a6)*cos(a1)),    1^-5*(3*(-3/6*cos(a0)+-3*e/6*sin(a0))*(cos(a2)*-3*e/6-sin(a2)*-3/6)) +1^-3*(cos(a0)*sin(a2)-sin(a0)*cos(a2))+1^-5*(3*(1/2*cos(a1)+-3*e/6*sin(a1))*(cos(a2)*-3*e/6-sin(a2)*1/2)) +1^-3*(cos(a1)*sin(a2)-sin(a1)*cos(a2))+1^-5*(3*(-6/6*cos(a3)+0*sin(a3))*(cos(a2)*0-sin(a2)*-6/6)) +1^-3*(cos(a3)*sin(a2)-sin(a3)*cos(a2))+e/27*(3*(-6/4*cos(a4)+-3*e/6*sin(a4))*(cos(a2)*-3*e/6-sin(a2)*-6/4)) +e/9*(cos(a4)*sin(a2)-sin(a4)*cos(a2))+2^-5*(3*(-6/6*cos(a5)+-6*e/6*sin(a5))*(cos(a2)*-6*e/6-sin(a2)*-6/6)) +2^-3*(cos(a5)*sin(a2)-sin(a5)*cos(a2))+e/27*(3*(0*cos(a6)+-6*e/6*sin(a6))*(cos(a2)*-6*e/6-sin(a2)*0)) +e/9*(cos(a6)*sin(a2)-sin(a6)*cos(a2)),    1^-5*(3*(1/2*cos(a0)+-3*e/6*sin(a0))*(cos(a3)*-3*e/6-sin(a3)*1/2)) +1^-3*(cos(a0)*sin(a3)-sin(a0)*cos(a3))+e/27*(3*(3/2*cos(a1)+-3*e/6*sin(a1))*(cos(a3)*-3*e/6-sin(a3)*3/2)) +e/9*(cos(a1)*sin(a3)-sin(a1)*cos(a3))+1^-5*(3*(1*cos(a2)+0*sin(a2))*(cos(a3)*0-sin(a3)*1)) +1^-3*(cos(a2)*sin(a3)-sin(a2)*cos(a3))+1^-5*(3*(-3/6*cos(a4)+-3*e/6*sin(a4))*(cos(a3)*-3*e/6-sin(a3)*-3/6)) +1^-3*(cos(a4)*sin(a3)-sin(a4)*cos(a3))+e/27*(3*(0*cos(a5)+-6*e/6*sin(a5))*(cos(a3)*-6*e/6-sin(a3)*0)) +e/9*(cos(a5)*sin(a3)-sin(a5)*cos(a3))+2^-5*(3*(1*cos(a6)+-6*e/6*sin(a6))*(cos(a3)*-6*e/6-sin(a3)*1)) +2^-3*(cos(a6)*sin(a3)-sin(a6)*cos(a3)),    1^-5*(3*(1*cos(a0)+0*sin(a0))*(cos(a4)*0-sin(a4)*1)) +1^-3*(cos(a0)*sin(a4)-sin(a0)*cos(a4))+2^-5*(3*(2*cos(a1)+0*sin(a1))*(cos(a4)*0-sin(a4)*2)) +2^-3*(cos(a1)*sin(a4)-sin(a1)*cos(a4))+e/27*(3*(3/2*cos(a2)+1*e/2*sin(a2))*(cos(a4)*1*e/2-sin(a4)*3/2)) +e/9*(cos(a2)*sin(a4)-sin(a2)*cos(a4))+1^-5*(3*(1/2*cos(a3)+1*e/2*sin(a3))*(cos(a4)*1*e/2-sin(a4)*1/2)) +1^-3*(cos(a3)*sin(a4)-sin(a3)*cos(a4))+1^-5*(3*(1/2*cos(a5)+-3*e/6*sin(a5))*(cos(a4)*-3*e/6-sin(a4)*1/2)) +1^-3*(cos(a5)*sin(a4)-sin(a5)*cos(a4))+e/27*(3*(3/2*cos(a6)+-3*e/6*sin(a6))*(cos(a4)*-3*e/6-sin(a4)*3/2)) +e/9*(cos(a6)*sin(a4)-sin(a6)*cos(a4)),    1^-5*(3*(1/2*cos(a0)+1*e/2*sin(a0))*(cos(a5)*1*e/2-sin(a5)*1/2)) +1^-3*(cos(a0)*sin(a5)-sin(a0)*cos(a5))+e/27*(3*(3/2*cos(a1)+1*e/2*sin(a1))*(cos(a5)*1*e/2-sin(a5)*3/2)) +e/9*(cos(a1)*sin(a5)-sin(a1)*cos(a5))+2^-5*(3*(1*cos(a2)+e*sin(a2))*(cos(a5)*e-sin(a5)*1)) +2^-3*(cos(a2)*sin(a5)-sin(a2)*cos(a5))+e/27*(3*(0*cos(a3)+e*sin(a3))*(cos(a5)*e-sin(a5)*0)) +e/9*(cos(a3)*sin(a5)-sin(a3)*cos(a5))+1^-5*(3*(-3/6*cos(a4)+1*e/2*sin(a4))*(cos(a5)*1*e/2-sin(a5)*-3/6)) +1^-3*(cos(a4)*sin(a5)-sin(a4)*cos(a5))+1^-5*(3*(1*cos(a6)+0*sin(a6))*(cos(a5)*0-sin(a5)*1)) +1^-3*(cos(a6)*sin(a5)-sin(a6)*cos(a5)),    1^-5*(3*(-3/6*cos(a0)+1*e/2*sin(a0))*(cos(a6)*1*e/2-sin(a6)*-3/6)) +1^-3*(cos(a0)*sin(a6)-sin(a0)*cos(a6))+1^-5*(3*(1/2*cos(a1)+1*e/2*sin(a1))*(cos(a6)*1*e/2-sin(a6)*1/2)) +1^-3*(cos(a1)*sin(a6)-sin(a1)*cos(a6))+e/27*(3*(0*cos(a2)+e*sin(a2))*(cos(a6)*e-sin(a6)*0)) +e/9*(cos(a2)*sin(a6)-sin(a2)*cos(a6))+2^-5*(3*(-6/6*cos(a3)+e*sin(a3))*(cos(a6)*e-sin(a6)*-6/6)) +2^-3*(cos(a3)*sin(a6)-sin(a3)*cos(a6))+e/27*(3*(-6/4*cos(a4)+1*e/2*sin(a4))*(cos(a6)*1*e/2-sin(a6)*-6/4)) +e/9*(cos(a4)*sin(a6)-sin(a4)*cos(a6))+1^-5*(3*(-6/6*cos(a5)+0*sin(a5))*(cos(a6)*0-sin(a6)*-6/6)) +1^-3*(cos(a5)*sin(a6)-sin(a5)*cos(a6))]
    return [x.n() for x in torque_eqs] # evaluate numerically

test = (1.973378658980664, 0.998083693363702, 6.01230155193543, 1.145904482094263, 6.72292920526922, 3.890451766786223, 5.16145340322297)
print(torque_array(test))
print(sum([x*x for x in torque_array(test)]))                                # sum of squares of torques computed from torque array
print(f(test[0], test[1], test[2], test[3], test[4], test[5], test[6]).n())  # same but from function used earlier
print(pot(test[0], test[1], test[2], test[3], test[4], test[5], test[6]).n())# potential energy
