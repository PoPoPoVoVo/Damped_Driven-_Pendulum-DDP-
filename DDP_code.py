#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  8 03:48:21 2022

@author: dohyun
"""

"""
proj3-2_ani.py                   Joonha Kim

Animation of damped driven pendulum, 
and plot phase diagram and  time series graph
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
    """ differential equations"""
    dydx = np.zeros_like(y0)
    dydx[0] = y0[1]
    dydx[1] = -  c*y0[1] - np.sin(y0[0]) + F* np.cos(omega_D*x)

    return dydx


def frequency(y):
    freq = np.array([-1.]*(n_steps+1))
    for i in range(n_steps):
        if y[i,1]*y[i+1,1]< 0.:
            a = (y[i+1,1] - y[i,1])/(x[i+1]-x[i])
            b = y[i,1] - a*x[i]
            freq[i]=-b/a
    j_temp = 0
    for j in range(n_steps+1):
        if freq[j]>0:
            freq[j_temp:j] = np.pi/(freq[j] - freq[j_temp])
            j_temp = j
    freq[j_temp:] = freq[j_temp-1]
    return freq

def poincare_plot(freq):
	return int(15)

A=0
n_steps =50000
initial_n_phase_diagram = int(1)
x0 = 0
xn = 1000
c = 0.05
F = 0.7
omega_D = 0.7

x = np.linspace(x0, xn, n_steps+1)
dx = (xn - x0)/n_steps

y0 = [A,0]

x = np.linspace(x0, xn, n_steps+1)
dx = (xn - x0)/n_steps


# initial state
y0 = np.array([A, 0])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, y0, x)
y_rk4 = rk4(derivs, y0, x)
xp = np.sin(y[:,0])
yp = -np.cos(y[:,0])
freq = frequency(y)
period = omega_D/freq
x_poin=[]
y_poin=[]


phi = y.copy()
for j in range(n_steps):
    while phi[j,0]>np.pi or phi[j,0] < -np.pi:
        if(phi[j,0]>np.pi): phi[j,0]= phi[j,0]-2*np.pi
        if(phi[j,0]<-np.pi):phi[j,0]= phi[j,0]+2*np.pi

fig= plt.figure(figsize=(10,10))
plt.suptitle("Damped Driven Pendulum    Dohyun Song")
ax = plt.subplot(212)#ax = plt.gca()
plt.xlim(0, xn)
plt.ylim(1.5*y[:,0].min(), 1.5*y[:,0].max())
param_text = ax.text(0.1,0.95, r"c=%.3f, F=%.4f, $\omega_D$=%.4f" 
                     % (c,F,omega_D), transform=ax.transAxes)

line,= plt.plot([], [],'r-', lw=0.5)
line2, = plt.plot([], [], 'b-', lw=0.5)
ax.grid()
ax.set_ylabel(r"$\phi$",rotation='horizontal')
ax.set_xlabel(r"$\tau$")

#진자 운동 구현
ax3 = plt.subplot(221)
line4, = ax3.plot([],[], 'o-', lw=2)
plt.xlim(-2,2)
plt.ylim(-2,2)
ax3.text(.1,0.95, "pendulum animation",transform=ax3.transAxes)
ax3.grid()

#위상 도표 구현
ax2 = plt.subplot(222)
line3, = ax2.plot([],[], 'r', marker='.',linestyle='',lw=0.5)
line5, = ax2.plot([],[], 'b', marker='.', linestyle ='',lw=2) #poincare plot
plt.xlim(-np.pi, np.pi)
plt.ylim(-5, 5)
ax2.set_xlabel(r"$\phi$",verticalalignment='top')
ax2.set_ylabel(r"$\dot{\phi}$", rotation='horizontal')
ax2.text(0.1,0.95, 'phase diagram',transform=ax2.transAxes)
ax2.grid()
time_template = 'time = %.1f,'
freq_template = r' frequency = %.5f, period = %.5f in $T_D$'
time_text = ax.text(0.43, 0.95, '', transform=ax.transAxes)
freq_text  =ax.text(0.55, 0.95, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([],[])
    line5.set_data([],[])
    time_text.set_text('')
    freq_text.set_text('')
    return line, line2, line3, line4,line5,time_text, freq_text


def animate(i):
     global x_poin, y_poin
     #print(omega_D*i*dx)%(2*np.pi)
     line.set_data(x[0:i], y[0:i,0])
     if i>initial_n_phase_diagram:
          line3.set_data(phi[initial_n_phase_diagram:i,0], phi[initial_n_phase_diagram:i,1])
     if i>initial_n_phase_diagram and (omega_D*i*dx)%(6.28)<0.02:
          x_poin.append(phi[i,0])
          y_poin.append(phi[i,1])
          line5.set_data(x_poin[:],y_poin[:])
     line2.set_data(x[0:i], y_rk4[0:i,0])
     line4.set_data([0,xp[i]],[0,yp[i]])
     time_text.set_text(time_template % (i*dx))
     freq_text.set_text(freq_template % (freq[i], period[i]))
     return line, line2, line3, line4,line5, time_text, freq_text
    


ani = animation.FuncAnimation(fig, animate, range(1, len(y)),
                              interval=dx, blit=True, init_func=init, repeat=False)

plt.show()
