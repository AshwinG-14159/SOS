#%%

import numpy as np
import matplotlib.pyplot as plt
import math as m



#%%
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

def acc_due_to_g(m,r,v,a):
    return a

def gravity_planet(m,r,v,a):
    #  print(-(G*m2/np.dot(r,r))*r)
    return -(G*m2/(np.dot(r,r))**1.5)*r


#%%

m2 = 6*10**24
mu = G*(m2)
min_dist = 36000*10**3
e = 1.5

r = [min_dist,0.,0.]
v = [0.,m.sqrt((mu)*(1+e)/min_dist),0.]


p1 = particle(m1,r,v,[0.,0.,0.], gravity_planet)

r_x,r_y,r_z,v_x,v_y,v_z,t = p1.simulate(100000,0.1)
#print(r_y[:1000])
plt.figure(figsize = (7,7))
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(-300000*10**3, 60000*10**3)
plt.ylim(-180000*10**3, 180000*10**3)

plt.scatter(r_x,r_y, s = 0.01)




#%%