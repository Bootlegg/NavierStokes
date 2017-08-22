import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d

def SolvePoisson(psi,zeta):
	"""
	We set BCs at y = top and y = bottom, but not x, i think
	Following lec12.pdf here
	"""
	#Compute streamfunction from vorticity
	#ζ_(i,j)^1=(((ψ_(i+1,j)-2ψ_(i,j)+ψ_(i-1,j) ))/dx^2 +((ψ_(i,j+1)-2ψ_(i,j)+ψ_(i,j-1) ))/dy^2 )

	#solve psi from lec12.pdf
	psin = np.zeros((ny,nx),dtype=float64)
	psin = psi.copy()
	
	#Skal måske lave en special case for enten x = 0 eller x = L, hvor jeg copy over?
	
	dtau = 0.5*0.5*(0.5*dx**2+0.5*dy**2)
	for r in range(500): #pseudo-time
		psin = psi.copy()
		psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(
				+(psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2
				+(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2
				-zeta[1:ny-1,1:nx-1])

		psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
				+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-1])/dx**2
				+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
				-zeta[1:ny-1,0])
		#psi[1:ny-1,-1] = psi[1:ny-1,0] 
		
		psi[1:ny-1,-1] = psin[1:ny-1,-1]+dtau*(
		+(psin[1:ny-1,0]-2*psin[1:ny-1,-1]+psin[1:ny-1,-2])/dx**2
		+(psin[2:ny,-1]-2*psin[1:ny-1,-1]+psin[0:ny-2,-1])/dy**2
		-zeta[1:ny-1,-1])
		
		#boundary at x = L, i'm gonna set equal to x = 0...
		#Boundary for psi, maybe i should remove these!
		#psi[:,-1] = 0 #right boundary
		#psi[:,0] = 0	#left boundary
		psi[-1,:] = 0
		psi[0,:] = 0 
		
		

Lx = 6000
Ly = 6000
nx = 65
ny = 65

x = np.linspace(0,Lx,nx)
y = np.linspace(0,Ly,ny)
x,y = np.meshgrid(x*1000,y*1000)
k = 2*pi/(Lx*1000)
m = pi/(Ly*1000)
dx = 1000*Lx/(nx-1)
dy = 1000*Ly/(ny-1)
U0 = 20				#zonal wind
beta = 1.62*10**(-11)	#he set 0 infront?
Av4 = 10**(-6)
A = 10**(-4)
#initial vorticity and streamfunction

zeta = np.array(A*np.exp(-2*(k**2*x**2+m**2*y**2)),dtype=float64)
zetan = zeta
#time integration parameters					#hours
time_end = 3*3600 				#second
dt = 100

psi = np.zeros((ny,nx),dtype=float64)
dypsi = np.zeros((ny,nx),dtype=float64)
dxpsi = np.zeros((ny,nx),dtype=float64)
u = np.zeros((ny,nx),dtype=float64)
v = np.zeros((ny,nx),dtype=float64)
dfly = np.zeros((ny,nx),dtype=float64)
dflx = np.zeros((ny,nx),dtype=float64)

SolvePoisson(psi,zeta)
#Compute streamfunction from vorticity
#ζ_(i,j)^1=(((ψ_(i+1,j)-2ψ_(i,j)+ψ_(i-1,j) ))/dx^2 +((ψ_(i,j+1)-2ψ_(i,j)+ψ_(i,j-1) ))/dy^2 )

# #solve psi from lec12.pdf
# psin = np.zeros((ny,nx),dtype=float64)
# psin = psi.copy()
# dtau = 0.5*0.5*(0.5*dx**2+0.5*dy**2)
# for r in range(500): #pseudo-time
	# psin = psi.copy()
	# psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(\
			# (psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2\
			# +(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2\
			# +zeta[1:ny-1,1:nx-1])


	# #Boundary for psi, maybe i should remove these!
	# #psi[:,-1] = 0 #right boundary
	# #psi[:,0] = 0	#left boundary
	# psi[-1,:] = 0
	# psi[0,:] = 0 

psi+=U0*(Ly/2*1000 - y)

#Calculate gradients for air velocity
#forward on bottom? Backward on top?
dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference


dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
dxpsi[:,-1]= (psi[:,0]-psi[:,-2])/(2*dx)
dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx);  #centered difference

u = -dypsi
v = dxpsi


#Forward time difference
zetan[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1]-dt*(\
	(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
	+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
	+beta*v[1:ny-1,1:nx-1])
zeta = zetan







#plt.hold(True)
fig = plt.figure()
ax = plt.gca()
ax.contour(x,y,zeta,colors='black')
ax.quiver(x,y,u,v)
ax.set_title('Barotropic Vorticity Equation')
ax.set_xlabel('x')
ax.set_ylabel('y')

def animate(i):

	global zetan,u,v,zeta,psi,psin,dxpsi,dypsi
	#leapfrog
	zetan[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1]-2*dt*(\
		(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
		+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
		+beta*v[1:ny-1,1:nx-1])
	zeta = zetan
	
	SolvePoisson(psi,zeta)
	# #Calculate psi streamfunction
	# #solve psi from lec12.pdf
	# psin = np.zeros((ny,nx),dtype=float64)
	# psin = psi.copy()
	# dtau = 0.5*0.5*(0.5*dx**2+0.5*dy**2)
	# for r in range(500): #pseudo-time
		# psin = psi.copy()
		# psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(\
				# (psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2\
				# +(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2\
				# +zeta[1:ny-1,1:nx-1])

		# #psi[:,-1] = 0 #right boundary
		# #psi[:,0] = 0	#left boundary
		# psi[-1,:] = 0
		# psi[0,:] = 0 

	psi+=U0*(Ly/2*1000-y)
 	
	#Calculate gradients for air velocity
	#forward on bottom? Backward on top?
	dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
	dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
	dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference


	dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
	dxpsi[:,-1]= (psi[:,0]-psi[:,-2])/(2*dx)
	dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx);  #centered difference

	u = -dypsi
	v = dxpsi

	ax.clear()
	ax.contour(x,y,zeta,colors='black')
	ax.quiver(x,y,u,v)
	ax.set_title('Barotropic Vorticity Equation')
	ax.set_xlabel('x')
	ax.set_ylabel('y')

	if i == 4:
		fig.savefig('BVE.png', bbox_inches='tight')



anim = animation.FuncAnimation(fig, animate, frames=5, interval=0.5)#,blit=False)
plt.show()