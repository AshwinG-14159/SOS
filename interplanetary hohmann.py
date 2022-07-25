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




#%%

#assumed r2>r1
plt.figure(figsize = (7,7))
m2 = 6*10**24
mu = G*(m2)
r1 = 18000*10**3
r2 = 25000*10**3
dt = 0.1
theta1 = m.radians(20)
theta2 = m.radians(80)

v1 = orbital_velocity(r1)
v2 = orbital_velocity(r2)



v1_vec = [-v1*m.sin(theta1),v1*m.cos(theta1),0.]
v2_vec = [-v2*m.sin(theta2),v2*m.cos(theta2),0.]


r1_vec = [r1*m.cos(theta1), r1*m.sin(theta1), 0.]
r2_vec = [r2*m.cos(theta2), r2*m.sin(theta2), 0.]

a_ellipse = a_val(r1,r2)
e_ellipse = e_val(r1,r2,a_ellipse)


p1 = particle(m1,r1_vec,v1_vec,[0.,0.,0.], gravity_planet)
p2 = particle(m1,r2_vec,v2_vec,[0.,0.,0.], gravity_planet)


omega1 = 2*m.pi/time_period(r1)
omega2 = 2*m.pi/time_period(r2)

omega_rel = omega2-omega1

t_transfer = time_period(a_ellipse)/2

delta_theta_req = m.pi - omega2*t_transfer

delta_theta_current = theta2-theta1

t_wait = (delta_theta_req - delta_theta_current)/omega_rel

#print(t_wait)
if t_wait < 0:
    t_wait += 2*m.pi/omega_rel


r_x1, r_y1, r_z1, v_x1, v_y1, v_z1,t_1 =  p1.simulate(t_wait, dt)
r_x2, r_y2, r_z2, v_x2, v_y2, v_z2,t_2 =  p2.simulate(t_wait, dt)

plt.scatter(r_x1, r_y1, s = 0.01)
plt.scatter(r_x2, r_y2, s = 0.01)
plt.show()

theta1 = theta1 + t_wait*omega1
theta2 = theta2 + t_wait*omega2

v1_vec = [-v1*m.sin(theta1),v1*m.cos(theta1),0.]
v4_vec = [-v2*m.sin(theta2),v2*m.cos(theta2),0.]

r_probe = [r1*m.cos(theta1), r1*m.sin(theta1), 0.]
v2 = v_perigee(a_ellipse, e_ellipse)
v3 = v_appogee(a_ellipse, e_ellipse)

v2_vec = [-v2*m.sin(theta1), v2*m.cos(theta1), 0.]
delta_v1 = np.array(v2_vec)-np.array(v1_vec)



r1_vec = [r1*m.cos(theta1), r1*m.sin(theta1), 0.]
r2_vec = [r2*m.cos(theta2), r2*m.sin(theta2), 0.]

p1 = particle(m1,r1_vec,v1_vec,[0.,0.,0.], gravity_planet)
p2 = particle(m1,r2_vec,v4_vec,[0.,0.,0.], gravity_planet)



p3 = particle(m1,r_probe, v1_vec, [0., 0., 0.], gravity_planet)
p3.delta_v(delta_v1)



r_x1b, r_y1b, r_z1b, v_x1b, v_y1b, v_z1b,t_1b =  p1.simulate(t_transfer, dt)
r_x2b, r_y2b, r_z2b, v_x2b, v_y2b, v_z2b,t_2b =  p2.simulate(t_transfer, dt)
r_x3, r_y3, r_z3, v_x3, v_y3, v_z3,t_3 =  p3.simulate(t_transfer, dt)

plt.figure(figsize = (7,7))

plt.scatter(r_x3, r_y3, s = 0.01)
plt.scatter(r_x1b, r_y1b, s = 0.01)
plt.scatter(r_x2b, r_y2b, s = 0.01)


'''
if r2>r1:
    v2 = [0.,v_perigee(a_ellipse, e_ellipse),0.]
    v3 = [0.,-v_appogee(a_ellipse, e_ellipse),0.]

else:
    v3 = [0.,-v_perigee(a_ellipse, e_ellipse),0.]
    v2 = [0.,v_appogee(a_ellipse, e_ellipse),0.]
'''


'''
delta_v1 = np.array(v2) - np.array(v1)
delta_v2 = np.array(v4) - np.array(v3)

p1.delta_v(delta_v1)

r_x2, r_y2, r_z2, v_x2, v_y2, v_z2, t2 = p1.simulate(time_period(a_ellipse)/2, 1)

p1.delta_v(delta_v2)

r_x3, r_y3, r_z3, v_x3, v_y3, v_z3, t3 = p1.simulate(time_period(r2), 1)

#print(r_x[:1000])
plt.figure(figsize = (7,7))



plt.scatter(r_x2,r_y2, s = 0.01)
plt.scatter(r_x3,r_y3, s = 0.01)
'''
plt.show()

#%%

print(len(r_x1))
print(np.array(r_x1[90:110])/1000000)
print(np.array(r_x1[-110:-90])/1000000)
print(np.array(r_x1b[:20])/1000000)


#%%
