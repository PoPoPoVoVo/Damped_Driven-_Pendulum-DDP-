"""
proj3-2_bifurcation.py                 Joonha Kim

Animation of physical pendulum gratphs
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


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

A=0
n_steps = 1000
c = 0.5
omega_D = 2./3
#initial_n_phase_diagram = int(n_steps/2)
x0 = 0
xn = 351*np.pi/omega_D #400
dx = (xn - x0)/n_steps
#ni = int(300*3*np.pi/2/dx + 0.5)

F0 = 1.42
Fn = 1.50
nF = 1000
Fx = np.linspace(F0,Fn,nF+1)

x = np.linspace(x0, xn, n_steps+1)
dx = (xn - x0)/n_steps

y0 = [A,0]

x = np.linspace(x0, xn, n_steps+1)

phi_dot_bi = np.zeros((nF+1,50))

# initial state
y0 = np.array([A, 0])
for iF in range(nF+1):
    F = Fx[iF]
    # integrate your ODE using scipy.integrate.
    y = rk4(derivs, y0, x) #integrate.odeint(derivs, y0, x)
    jj = 0
    for ii in range(n_steps+1):
        if x[ii] > (301+jj)*np.pi/omega_D:
            a = (y[ii,1]-y[ii-1,1])/(x[ii]-x[ii-1])
            b = y[ii,1] - a*x[ii]
            phi_dot_bi[iF,jj] = ((-1)**(jj+1))*(a*(301+jj)*np.pi/omega_D+b)
            jj += 1
        if jj>=50: break
            

#y_rk4 = rk4(derivs, y0, x)

fig= plt.figure(figsize=(10,8))
plt.title('Bifurcation diagram - Damped driven pendulum          Joonha Kim')
ax = plt.gca()
plt.xlim(F0, Fn)
plt.ylim(-1, 4)
plt.grid()
ax.set_xlabel('F')
ax.set_ylabel(r'$\dot{\phi}$', rotation='horizontal')
param_text = ax.text(0.2, 0.95, r"c=%.3f, $\omega_D$=%.4f" % (c,omega_D), transform=ax.transAxes)
for iii in range(50):
    plt.plot(Fx,phi_dot_bi[:,iii],'r,')
plt.savefig('proj3-2_bifurcation%.2f_%.2f.png' % (F0, Fn))
plt.show()
