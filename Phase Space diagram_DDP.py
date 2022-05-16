
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


A=0
n_steps =500000
initial_n_phase_diagram = int(n_steps/2)
x0 = 0
xn = 10000
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

x_poin=[]
y_poin=[]


phi = y.copy()
for j in range(n_steps):
    while phi[j,0]>np.pi or phi[j,0] < -np.pi:
        if(phi[j,0]>np.pi): phi[j,0]= phi[j,0]-2*np.pi
        if(phi[j,0]<-np.pi):phi[j,0]= phi[j,0]+2*np.pi

fig= plt.figure(figsize=(10,8))
plt.title('Phase Space diagram - Damped driven pendulum          Dohyun Song')
ax = plt.gca()
plt.xlim(-np.pi, np.pi)
plt.ylim(-5, 5)
plt.grid()
ax.set_xlabel(r'$\phi$',verticalalignment='top')
ax.set_ylabel(r'$\dot{\phi}$', rotation='horizontal')

'''
for i in range(1,len(y)):
    if (omega_D*i*dx)%(6.28)<0.001 and i>initial_n_phase_diagram:
        x_poin.append(phi[i,0])
        y_poin.append(phi[i,1])
'''        
        
#print(len(x_poin))
plt.plot(phi[initial_n_phase_diagram:,0], phi[initial_n_phase_diagram:,1],'r', marker='.',linestyle='',lw=0.5)
#plt.plot(x_poin[:],y_poin[:],'b', marker='.', linestyle ='',lw=2)

plt.savefig('Phase Space diagram_DDP_%.2f.png' %F)
plt.show()

