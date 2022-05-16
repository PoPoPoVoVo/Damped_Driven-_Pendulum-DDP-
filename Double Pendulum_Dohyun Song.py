#!/usr/bin/env python3
 # -*- coding: utf-8 -*-
"""
Created on Wed May  4 17:51:44 2022

@author: dohyun
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

def rk4(func, y0, x):
    y0 = np.array(y0)
    n = len(x)
    f = lambda xi,yi: np.array(func(yi,xi))
    y = np.zeros( (n,) + y0.shape)
    h = (x[n-1] - x[0])/(n-1) 
    xi = x[0]
    y[0] = yi = y0[:]
    for i in range(1, n):
        k1 = h * f(xi, yi)
        k2 = h * f(xi + 0.5 * h, yi + 0.5 * k1)
        k3 = h * f(xi + 0.5 * h, yi + 0.5 * k2)
        k4 = h * f(xi + h, yi + k3)
        xi = x[i]
        y[i] = yi = yi + (k1 + 2*(k2 + k3) + k4) / 6
    return y

def derivs(y0, x):
	dydx = np.zeros_like(y0)
	delta = y0[0] - y0[2]
	a1 = L2/L1 * (M2/(M1+M2)) * np.cos(delta)
	a2 = L1/L2 * np.cos(delta)
	f1 = -(L2/L1) * (M2/(M1+M2)) * y0[3]*y0[3]*np.sin(delta) - (g/L1)*np.sin(y0[0])
	f2 = (L1/L2) *  y0[1] *y0[1] * np.sin(delta)-(g/L2)*np.sin(y0[2])
	
	dydx[0]=y0[1]
	dydx[1]=(f1-a1*f2)/(1-a1*a2)
	dydx[2]=y0[3]
	dydx[3]=(-a2*f1+f2)/(1-a1*a2)
	return dydx

L1 = 1
L2 = 1
M1 =1
M2 =1
g = 9.8


y0 = [60,0,30,0]

n_steps = 20000
x0 = 0
xn = 400
x = np.linspace(x0, xn, n_steps+1)
dx = (xn - x0)/n_steps



y = integrate.odeint(derivs, y0, x)
y_rk4 = rk4(derivs, y0, x)
xp1 = L1* np.sin(y[:,0])
yp1 = -L2*np.cos(y[:,0])
xp2 = xp1+ L2*np.sin(y[:,2])
yp2 = yp1+ -L2*np.cos(y[:,2])


fig=plt.figure(figsize=(10,10))
plt.suptitle("Double Pendulum     Dohyun Song")
ax=plt.axes()
plt.xlim(-3,3)
plt.ylim(-3,3)
line1, = ax.plot([],[],'o-',lw=2)
line2, = ax.plot([],[],'o-',lw=2)
ax.grid()

def init():
	line1.set_data([],[])
	line2.set_data([],[])
	
	return line1, line2

def animate(i):
	line1.set_data([0,xp1[i]],[0,yp1[i]])
	line2.set_data([xp1[i],xp2[i]],[yp1[i],yp2[i]])
	
	return line1, line2


ani = animation.FuncAnimation(fig, animate,range(1,len(y)),
							  interval=dx, blit=True, init_func = init, repeat=False)

plt.show()
print(xp1)

