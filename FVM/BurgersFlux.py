from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np

dx = 0.005
dt = 0.5*dx
x0 = 0
x1 = 1
Nx = int(1+(x1-x0)/dx)
print(Nx)
Nt = 1500
psi = np.zeros((Nt,Nx))

#psi[0,0:int(Nx/2)] = 1 #To the left we have 1, to the right we have 0

psi[:,0] = 1
psi2 = np.zeros((Nt,Nx))
#psi2[0,0:int(Nx/2)] = 1 #To the left we have 1, to the right we have 0
psi2[:,0] = 1

psi3 = np.zeros((Nt,Nx))
#psi2[0,0:int(Nx/2)] = 1 #To the left we have 1, to the right we have 0
psi3[:,0] = 1

def UpstreamAdvection():
	for nt in range(Nt-1):
		for i in range(1,Nx):
			psi[nt+1,i] = psi[nt,i] - dt*((psi[nt,i]+psi[nt,i-1])/2)**2*(psi[nt,i]-psi[nt,i-1])/dx

def UpstreamAdvection2():
	for nt in range(Nt-1):
		for i in range(1,Nx):
			psi3[nt+1,i] = psi3[nt,i] - dt*psi3[nt,i]**2*(psi3[nt,i]-psi3[nt,i-1])/dx	
			
def UpstreamFlux():
	for nt in range(Nt-1):
		for i in range(1,Nx):
			psi2[nt+1,i] = psi2[nt,i] - dt*(psi2[nt,i]**3-psi2[nt,i-1]**3)/(3*dx)
			
UpstreamAdvection()
UpstreamAdvection2()
UpstreamFlux()

xs = np.linspace(x0,x1,Nx)
time = int(2.4/dt)
plt.title("time = {}".format(time*dt))
plt.plot(xs,psi[time,:],color="black")
plt.plot(xs,psi2[time,:],color="red",linestyle="--")
plt.plot(xs,psi3[time,:],color="blue",linestyle="--")
plt.legend(["Up,Advec","Up,Advec2","Up,Flux"])
plt.axes([0.5,1,-0.1,1.1])
plt.show()



