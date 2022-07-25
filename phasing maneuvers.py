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
    
    def delta_v(self,d_v):
        self.v += d_v
    
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

def orbital_velocity(r):
    return m.sqrt(G*m2/r)

def v_perigee(a,e):
    return m.sqrt(abs((mu/a)*(1+e)/(1-e)))

def v_appogee(a,e):
    return m.sqrt(abs((mu/a)*(1-e)/(1+e)))

def a_val(r1,r2):
    return (r1+r2)/2

def e_val(r1,r2,a):
    return abs(r2-r1)/(2*a)

def time_period(a):
    return m.pi*2* m.sqrt((a**3)/mu)

def a_given_t(T):
    return (T*m.sqrt(mu)/(2*m.pi))**(2/3)

#%%

#p1 moves and cathes p2

m2 = 6*10**24
mu = G*(m2)

r = 36000*10**3
theta = m.radians(20)
dt = 1

v = orbital_velocity(r)

r1 = [r,0.,0.]
r2 = [r*m.cos(theta), r*m.sin(theta), 0.]
v1 = [0., v, 0.]
v2 = [-v*m.sin(theta), v*m.cos(theta), 0.]

p1 = particle(m1,r1,v1,[0.,0.,0.], gravity_planet)
p2 = particle(m1,r2,v2,[0.,0.,0.], gravity_planet)

r_x1a,r_y1a,r_z1a,v_x1a,v_y1a,v_z1a,t1a = p1.simulate(time_period(r),dt)
r_x2a,r_y2a,r_z2a,v_x2a,v_y2a,v_z2a,t2a = p2.simulate(time_period(r),dt)

T = time_period(r)
delta_t = T*theta/(2*m.pi)
T_req = T-delta_t
a_req = a_given_t(T_req)


r_other = 2*a_req - r
e = e_val(r,r_other,a_req)

v_req = v_perigee(a_req, e)
dv = np.array([0.,v_req,0.]) - np.array(v1)

p1.delta_v(-dv)
r_x1b,r_y1b,r_z1b,v_x1b,v_y1b,v_z1b,t1b = p1.simulate(time_period(a_req),dt)
r_x2b,r_y2b,r_z2b,v_x2b,v_y2b,v_z2b,t2b = p2.simulate(time_period(r),dt)

p1.delta_v(dv)

r_x1c,r_y1c,r_z1c,v_x1c,v_y1c,v_z1c,t1c = p1.simulate(time_period(r),dt)
r_x2c,r_y2c,r_z2c,v_x2c,v_y2c,v_z2c,t2c = p2.simulate(time_period(r),dt)

plt.scatter(r_x1b,r_y1b, c = 'blue', s = 0.1)
plt.scatter(r_x2b,r_y2b, c = 'red', s = 0.1)
plt.show()




#%%

