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


def t_hohmann(r1,r2):
    return time_period(a_val(r1,r2))/2

def t_bi_elliptical(r1,r2,r3):
    return time_period(a_val(r1,r2))/2  +  time_period(a_val(r2,r3))/2

def max_r_bi(r1,r3,t0):
    r2 = max(r1,r3)
    dr = min(r1,r3)/1000
    while t_bi_elliptical(r1, r2, r3) < t0:
        r2+=dr
    
    return r2-dr


def dv_hohmann(alpha, r1):
    return (1/m.sqrt(alpha) - m.sqrt(2)*(1-alpha)/m.sqrt(alpha*(1+alpha)) - 1) * m.sqrt(mu/r1)
    
    

def dv_bi(alpha, beta, r1):
    return ( m.sqrt(2*(alpha+beta)/(alpha*beta))  -  1 - 1/m.sqrt(alpha)  - m.sqrt(2)*(1-beta)/m.sqrt(beta*(1+beta)) ) * m.sqrt(mu/r1)



def efficient_path(r1,r2,t_max):
    
    if t_max < t_hohmann(r1,r2):
        return -1
    elif t_max <= 2*t_hohmann(r1,r2):
        return 'hohmann'
    
    else:
        r_bi_max = max_r_bi(r1,r2,t_max)
        if r_bi_max < r2:
            return 'hohmann'
    alpha = r2/r1
    beta = r_bi_max/r1
    
    if dv_hohmann(alpha, r1) > dv_bi(alpha, beta, r1):
        return 'bi-elliptic', r_bi_max
    
    return 'hohmann'


#%%

def hohmann_transfer(r1,r2,dt):
    v1 = [0.,orbital_velocity(r1),0.]
    v4 = [0.,-orbital_velocity(r2),0.]
    
    
    r = [r1,0.,0.]
    #v = [0.,orbital_velocity(r1),0.]
    
    
    p1 = particle(m1,r,v1,[0.,0.,0.], gravity_planet)
    
    r_x,r_y,r_z,v_x,v_y,v_z,t = p1.simulate(time_period(r1),dt)
    
    a_ellipse = a_val(r1,r2)
    e_ellipse = e_val(r1,r2,a_ellipse)
    
    if r2>r1:
        v2 = [0.,v_perigee(a_ellipse, e_ellipse),0.]
        v3 = [0.,-v_appogee(a_ellipse, e_ellipse),0.]
    
    else:
        v3 = [0.,-v_perigee(a_ellipse, e_ellipse),0.]
        v2 = [0.,v_appogee(a_ellipse, e_ellipse),0.]
    
    delta_v1 = np.array(v2) - np.array(v1)
    delta_v2 = np.array(v4) - np.array(v3)
    
    p1.delta_v(delta_v1)
    
    r_x2, r_y2, r_z2, v_x2, v_y2, v_z2, t2 = p1.simulate(time_period(a_ellipse)/2, dt)
    
    p1.delta_v(delta_v2)
    
    r_x3, r_y3, r_z3, v_x3, v_y3, v_z3, t3 = p1.simulate(time_period(r2), dt)
    
    #print(r_x[:1000])
    plt.figure(figsize = (7,7))
    
    
    plt.scatter(r_x,r_y, s = 0.01)
    plt.scatter(r_x2,r_y2, s = 0.01)
    plt.scatter(r_x3,r_y3, s = 0.01)
    plt.show()



#%%

def bi_elliptic_transfer(ra,rb,rc,dt):
    va1 = [0.,orbital_velocity(ra),0.]
    vc2 = [0.,orbital_velocity(rc),0.]
    
    
    r = [ra,0.,0.]
    
    
    p1 = particle(m1,r,va1,[0.,0.,0.], gravity_planet)
    
    r_x,r_y,r_z,v_x,v_y,v_z,t = p1.simulate(time_period(ra),dt)
    
    a_ellipse1 = a_val(ra,rb)
    e_ellipse1 = e_val(ra,rb,a_ellipse1)
    
    a_ellipse2 = a_val(rb,rc)
    e_ellipse2 = e_val(rb,rc,a_ellipse2)
    
    
    va2 = [0.,v_perigee(a_ellipse1, e_ellipse1),0.]
    vb1 = [0.,-v_appogee(a_ellipse1, e_ellipse1),0.]
    vb2 = [0.,-v_appogee(a_ellipse2, e_ellipse2),0.]
    vc1 = [0.,v_perigee(a_ellipse2, e_ellipse2),0.]
    
    
    delta_v1 = np.array(va2) - np.array(va1)
    delta_v2 = np.array(vb2) - np.array(vb1)
    delta_v3 = np.array(vc2) - np.array(vc1)
    
    
    
    p1.delta_v(delta_v1)
    
    r_x2, r_y2, r_z2, v_x2, v_y2, v_z2, t2 = p1.simulate(time_period(a_ellipse1)/2, dt)
    
    p1.delta_v(delta_v2)
    
    r_x3, r_y3, r_z3, v_x3, v_y3, v_z3, t3 = p1.simulate(time_period(a_ellipse2)/2, dt)
    
    p1.delta_v(delta_v3)
    
    r_x4, r_y4, r_z4, v_x4, v_y4, v_z4, t4 = p1.simulate(time_period(rc), dt)
    
    
    #print(r_x[:1000])
    plt.figure(figsize = (7,7))
    
    
    plt.scatter(r_x,r_y, s = 0.01)
    plt.scatter(r_x2,r_y2, s = 0.01)
    plt.scatter(r_x3,r_y3, s = 0.01)
    plt.scatter(r_x4,r_y4, s = 0.01)
    plt.show()


#%%

m2 = 6*10**24
mu = G*(m2)
r1 = 18000*10**3
r2 = 10 * 40000*10**3
dt = 5


t_max = 2.5*time_period(a_val(r1,r2))


if r2>r1:
    path = efficient_path(r1,r2,t_max)
else:
    path = efficient_path(r2,r1,t_max)

print(path)

if path == 'hohmann':
    hohmann_transfer(r1, r2, dt)
    
elif path != -1:
    bi_elliptic_transfer(r1, path[1], r2, dt)



#%%



#%%
