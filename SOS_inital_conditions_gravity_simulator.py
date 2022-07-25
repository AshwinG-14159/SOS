#%%

import numpy as np
import matplotlib.pyplot as plt
import math as m

#%%
#single body motion

class particle():
    def __init__(self, m, r, v, a, f):
        self.r = np.array(r)
        self.v = np.array(v)
        self.a = np.array(a)
        self.m = m
        self.f = f
    
    def set_acc(self, m, r, v, a):
        return(self.f(m, r, v, a))
    
    
    def update(self,dt):
        self.v[0] += self.a[0]*dt/2
        self.v[1] += self.a[1]*dt/2
        self.v[2] += self.a[2]*dt/2
        
        self.r[0] += self.v[0]*dt
        self.r[1] += self.v[1]*dt
        self.r[2] += self.v[2]*dt
        
        self.v[0] += self.a[0]*dt/2
        self.v[1] += self.a[1]*dt/2
        self.v[2] += self.a[2]*dt/2

        self.a = self.set_acc(self.m, self.r, self.v, self.a)
        
    def simulate(self,duration,dt):
        r_arr1 = [self.r[0]]
        
        r_arr2 = [self.r[1]]
        r_arr3 = [self.r[2]]
        v_arr1 = [self.v[0]]
        v_arr2 = [self.v[1]]
        v_arr3 = [self.v[2]]
        t_arr = [0]

        for i in range(1, int(duration/dt)+1):
            self.update(dt)
           # print(self.r[0])
            r_arr1.append(self.r[0])
            r_arr2.append(self.r[1])
            r_arr3.append(self.r[2])
            v_arr1.append(self.v[0])
            v_arr2.append(self.v[1])
            v_arr3.append(self.v[2])
            t_arr.append(i)

        return(r_arr1,r_arr2,r_arr3,v_arr1,v_arr2,v_arr3,t_arr)
#%%

#forces

G = 6.67 * 10**-11
m1 = 10
m2 = 6 * 10**24
mu = G*(m2)


def gravity_earth(m,r,v,a):
    #  print(-(G*m2/np.dot(r,r))*r)
    return -(G*m2/(np.dot(r,r))**1.5)*r

#%%

#Calculations

def orbital_velocity(r, m2):
    return m.sqrt(G*m2/r)

def angular_momentum(r,v):
    return np.cross(r,v)

def eccentricity_vector(r,v,h,mu):
    e = (np.cross(v,h))/mu - np.array(r)/m.sqrt(np.dot(r,r))
    return e

def mag(r):
    return m.sqrt(np.dot(r,r))

def average_of_arrays(arr):
    l = []
    for i in range(len(arr[0])):
        sum = 0
        for j in range(len(arr)):
            sum += arr[j,i]
        l.append(sum/len(arr))
    
    return l

def angle_bw(r1,r2):
    return np.dot(r1,r2)/(mag(r1) * mag(r2))


#%%

#simulate a circular orbit around earth

r = 360*10**3.


p1 = particle(10, [r,0.,0.], [0.3*orbital_velocity(r,6 * 10**24),1.2*orbital_velocity(r,6 * 10**24),0.], [0.,0.,0.], gravity_earth)

r_x,r_y,r_z,v_x,v_y,v_z,t = p1.simulate(200,0.001)
plt.figure(figsize = (7,7))
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-1.1*10**6, 1*10**6)
plt.ylim(-1.1*10**6, 1*10**6)
plt.scatter(r_x,r_y, s = 0.01)


#%%
#calculating angular momentum per unit mass


l_h = np.zeros((len(r_x),3))

for i in range(len(r_x)):
    l_h[i] = angular_momentum([r_x[i], r_y[i], r_z[i]],[v_x[i], v_y[i], v_z[i]])

plt.figure(figsize = (7,7))
plt.xlabel('time')
plt.ylabel('angular momentum')
plt.plot(t, l_h[:,2])

 #variation of angular momentum with time is very small
#smaller the timestep, smaller the variation in the value of h over the motion. Hence it can be used to determine how accurate a simulation is
#have proved that angular momentum stays constant over the motion


#%%
#calculating the eccentricity vector
#It is also found to be approximately constant
l_e =  np.zeros((len(r_x),3))

for i in range(len(r_x)):
    l_e[i] = eccentricity_vector([r_x[i],r_y[i],r_z[i]], [v_x[i],v_y[i],v_z[i]], l_h[i], mu)

plt.figure(figsize = (7,7))
plt.xlabel('time')
plt.ylabel('eccentricity vector (y component)')
plt.scatter(t, l_e[:,1], s = 0.01)


#%%

#finding actual values of h, apse line vector, eccentricity for the motion

h = average_of_arrays(l_h)
apse = average_of_arrays(l_e)
e = mag(apse)
print(h,apse,e)


#%%


#%%




#%%






