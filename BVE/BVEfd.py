"""
Dan Krog

Numerical Solution of Barotropic Vorticity Equation

vorticity, Holton page 101,
Absolute vorticity eta,

eta = zeta + f

f is coriolis parameter, zeta is relative vorticity, zeta = nabla x (u,v)'

dzeta/dt + J(psi,zeta+f) = 0

where J is Jacobian operator


Barotropic fluid, pressure depends only only density.
Hence, we can take the shallow water equations, and make the height constant, h and z are constant,
then, vertical velocity w = 0, and the fluid is horizontally nondivergent,
du/dx + dv/dy = 0


Boundary conditions, for a periodic in x, wall in ybottom and ytop
psi = 0
nablapsi.n = 0 (no-slip on side walls)


DIfferent methods to solve the Poisson equation
Should add some terms to the BVE, like vorticity drag, forcing, viscosity



====================================================

Change functions to return


Perhaps instead of setting psi = 0, see if it changes if it set psi = psi_old on northern and southern boundary
for Poisson solver. and hold psi_old constant throughout pseudo-time iteration on nothern and southern boundaries


lige nu holder jeg zeta = 0 på boundaries, maybe
i hvert fald gøre den periodic etc?


maybe er min velocity fra stream function psi derivatives også forkerte


lige nu er beta constant, bør den vidst ikke være, right?


#Mangler denne her ikke noget? Måske?Passer det her med units? Ikke særlig godt tror jeg?
psi+=U0*(Ly/2*1000-y)
jo det passer vel med units faktisk, men, hvorfor er den der?


I may have translated matlab code wrong, forgot python numpy slicing doesn't include last point, so x[1:nx-1] e.g


Add in Arawakian grid?


"""


import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d

#0 = zeta0
#1 -> zeta0
#2 -> 0, 0 = zeta0, 1 = zeta
#3 -> 1, 1 = zeta0, 2 = zeta
#4 -> 2, 2 = zeta0, 3 = zeta
#5 -> 3, 3 = zeta0, 4 = zeta

def SolvePoisson(psi,zeta):
	"""
	Solve Poisson Equation
	Holton page 465
	
	nabla^2psi = zeta
	
	where we are solving for psi the streamfunction, and zeta is the relative vorticity
	
	
	We set BCs at y = top and y = bottom, but not x, i think
	Following lec12.pdf here
	
	
	From lec12.pdf,
	Solving poisson equation means putting some boundary condition on psi, right?
	Because otherwise you could add any constant to psi, and it would still be able to give some curvature
	nabla^2psi = f...
	so it makes sense, to give some unique solution, we impose boundary conditions, 
	but i think, for the sake of calculating velocity, it also wouldn't matter if we added some constant to psi
	
	
	right now i'm doing psi[1:ny-1,1:nx-1]
	#So that means i'm not setting psi[bottom,x] but I AM setting psi[top,x]
	#So i think it should be psi[1:ny-2,1:nx-1]
	#Or maybe, x should also be changed...
	But, remember in Python, 1:ny-1, then 1 is included, but ny-1 is excluded, so in effect, we're doing 1:ny-2
	
	psi is streamfunction
	zeta is vorticity
	"""
	#Compute streamfunction from vorticity
	#ζ_(i,j)^1=(((ψ_(i+1,j)-2ψ_(i,j)+ψ_(i-1,j) ))/dx^2 +((ψ_(i,j+1)-2ψ_(i,j)+ψ_(i,j-1) ))/dy^2 )

	#solve psi from lec12.pdf
	psin = np.zeros((ny,nx),dtype=np.float64)
	psin = psi.copy()
	
	
	dtau = 0.5*0.5*(0.5*dx**2+0.5*dy**2)
	for r in range(500): #pseudo-time
		psin = psi.copy()
		
		
		#Interior points
		psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(
				+(psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2
				+(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2
				-zeta[1:ny-1,1:nx-1])
		
		
		#x = 0 boundary
		#Jeg tror, vi skal tage et note ud af Holtons bog, og her, skipper vi psin[1:ny-1,-1] og bruger psin[1:ny-1,-2]
		psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
				+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-2])/dx**2
				+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
				-zeta[1:ny-1,0])
		#Old
		#psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
		#		+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-1])/dx**2
		#		+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
		#		-zeta[1:ny-1,0])		
		
		#x = L boundary
		#Enten så gør den periodic sådan her
		#psi[1:ny-1,-1] = psi[1:ny-1,0]
		#eller sådan her, Holton
		psi[:,-1] = psi[:,0]
		
		#Eller gør den periodic sådan her
		#psi[1:ny-1,-1] = psin[1:ny-1,-1]+dtau*(
		#+(psin[1:ny-1,0]-2*psin[1:ny-1,-1]+psin[1:ny-1,-2])/dx**2
		#+(psin[2:ny,-1]-2*psin[1:ny-1,-1]+psin[0:ny-2,-1])/dy**2
		#-zeta[1:ny-1,-1])
		
		#boundary at x = L, i'm gonna set equal to x = 0...
		#Boundary for psi, maybe i should remove these!
		#psi[:,-1] = 0 #right boundary
		#psi[:,0] = 0	#left boundary
		
		
		#What if i don't set it to 0? afaik it should be constant along the edges, but maybe not 0...
		#Det har i hvert fald bestemt en effekt om man sætter det til 0 eller ej..
		#psi[0,:] = 0 
		#psi[-1,:] = 0
	
		
def FirstStepZeta(zeta,zetan,u,v,beta):
	"""
	Euler-forward probably?
	
	Denne her kan jeg indføre nogle [:] til, i stedet for[1:ny-1] fx.... tror godt vi kan gå HELT UD til boundaries
	Men lad os starte et sted dog
	
	When we calculate zetan for x = 0, we use x -> -2 instead of x -> -1
	"""
	
	#Interior points
	zetan[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1]-dt*(\
		(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
		+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
		+beta*v[1:ny-1,1:nx-1])
	
	#x = 0
	zetan[1:ny-1,0] = zeta[1:ny-1,0]-dt*(\
		(u[1:ny-1,1]*zeta[1:ny-1,1]-u[1:ny-1,-2]*zeta[1:ny-1,-2])/(2*dx)\
		+(v[2:ny,0]*zeta[2:ny,0]-v[0:ny-2,0]*zeta[0:ny-2,0])/(2*dy)\
		+beta*v[1:ny-1,0])
	
	#x = L borders
	#zetan[1:ny-1,-1] = zeta[1:ny-1,-1]-dt*(\
	#	(u[1:ny-1,0]*zeta[1:ny-1,0]-u[1:ny-1,-2]*zeta[1:ny-1,-2])/(2*dx)\
	#	+(v[2:ny,-1]*zeta[2:ny,-1]-v[0:ny-2,-1]*zeta[0:ny-2,-1])/(2*dy)\
	#	+beta*v[1:ny-1,-1])
	
	zetan[:,-1] = zetan[:,0]
	
	#zetan[1:ny-1,1] = zeta[1:ny-1,0]-dt*(\
	#	(u[1:ny-1,1]*zeta[1:ny-1,1]-u[1:ny-1,-1]*zeta[1:ny-1,-1])/(2*dx)\
	#	+(v[2:ny,0]*zeta[2:ny,0]-v[0:ny-2,0]*zeta[0:ny-2,0])/(2*dy)\
	#	+beta*v[1:ny-1,0])
	
	#zeta = zetan		

	
def UpdateZetaLeapFrog(zetan,zeta0,zeta,u,v,beta):

	"""
	Denne her kan jeg indføre nogle [:] til, i stedet for[1:ny-1] fx.... tror godt vi kan gå HELT UD til boundaries
	Men lad os starte et sted dog
	
	Hmm... jeg skal faktisk også have old velocities?:O
	"""
	zetan[1:ny-1,1:nx-1] = zeta0[1:ny-1,1:nx-1]-2*dt*(\
		(u[1:ny-1,2:nx]*zeta0[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta0[1:ny-1,0:nx-2])/(2*dx)\
		+(v[2:ny,1:nx-1]*zeta0[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta0[0:ny-2,1:nx-1])/(2*dy)\
		+beta*v[1:ny-1,1:nx-1])
	
	#x = 0, x = L borders
	zetan[1:ny-1,0] = zeta0[1:ny-1,0]-2*dt*(\
		(u[1:ny-1,1]*zeta0[1:ny-1,1]-u[1:ny-1,-1]*zeta0[1:ny-1,-1])/(2*dx)\
		+(v[2:ny,0]*zeta0[2:ny,0]-v[0:ny-2,0]*zeta0[0:ny-2,0])/(2*dy)\
		+beta*v[1:ny-1,0])
	
	
	zetan[1:ny-1,-1] = zeta0[1:ny-1,-1]-2*dt*(\
		(u[1:ny-1,0]*zeta0[1:ny-1,0]-u[1:ny-1,-2]*zeta0[1:ny-1,-2])/(2*dx)\
		+(v[2:ny,-1]*zeta0[2:ny,-1]-v[0:ny-2,-1]*zeta0[0:ny-2,-1])/(2*dy)\
		+beta*v[1:ny-1,-1])
	

	
	#Holton does this
	#zetan[:,nx-1] = zetan[:,0] #meget unstable
	
	
	
	#y = 0, y = L borders
	zetan[0,1:nx-1] = zeta0[0,1:nx-1]-2*dt*(\
		(u[0,2:nx]*zeta0[0,2:nx]-u[0,0:nx-2]*zeta0[0,0:nx-2])/(2*dx)\
		+(v[1,1:nx-1]*zeta0[1,1:nx-1]-v[-1,1:nx-1]*zeta0[-1,1:nx-1])/(2*dy)\
		+beta*v[0,1:nx-1])
	
	
	zetan[-1,1:nx-1] = zeta0[-1,1:nx-1]-2*dt*(\
		(u[-1,2:nx]*zeta0[-1,2:nx]-u[-1,0:nx-2]*zeta0[-1,0:nx-2])/(2*dx)\
		+(v[0,1:nx-1]*zeta0[0,1:nx-1]-v[-2,1:nx-1]*zeta0[-2,1:nx-1])/(2*dy)\
		+beta*v[-1,1:nx-1])
	
	
	#zeta0 = zetan
	
def UpdateZetaLeapFrogHolton(zetan,zeta0,zeta,u,v,beta,numdif):
	"""

	"""
	
	#zetan[1:ny-2,0:nx-2] = (zeta0[1:ny-2,0:nx-2] -beta*2*dt*v[1:ny-2,0:nx-2]	
	#						-2*dt*(dflx[1:ny-2,0:nx-2]+dfly[1:ny-2,0:nx-2])
     #       -2*dt*numdif[1:ny-2,0:nx-2])
			
	zetan[1:ny-1,0:nx-1] = (zeta0[1:ny-1,0:nx-1] -beta*2*dt*v[1:ny-1,0:nx-1]	
							-2*dt*(dflx[1:ny-1,0:nx-1]+dfly[1:ny-1,0:nx-1])
            -2*dt*numdif[1:ny-1,0:nx-1])
	
	zetan[:,nx-1]=zetan[:,0]

def divflux(P,u,v,dx,dy):
	dflx = np.zeros((ny,nx))
	dfly = np.zeros((ny,nx))
	#FØRST dfly
	#han har sat 0 foran y divflux, så vi har ingen y directioon divflux langs I bottom og top
	#dfly[0,:] = 0*(P[1,:]*v[1,:] - P[0,:]*v[0,:])/dy;
	#dfly[ny-1,:] = 0*(P[ny-1,:]*v[ny-1,:] - P[ny-2,:]*v[ny-2,:])/dy;
	#I center:
	#dfly[1:ny-2,:] = (P[2:ny-1,:]*v[2:ny-1,:]-P[0:ny-3,:]*v[0:ny-3,:])/(2*dy);

	#NU TAGER VI dflx
	#% Take cyclic differences on left and right boundaries

	#dflx[:,0]=(P[:,1]*u[:,1]-P[:,nx-2]*u[:,nx-2])/(2*dx);
	#dflx[:,nx-1]= dflx[:,0];

	#% take centered differences on interior points

	#dflx[:,1:nx-2]= (P[:,2:nx-1]*u[:,2:nx-1]-P[:,0:nx-3]*u[:,0:nx-3])/(2*dx);
	
	
	
	#FØRST dfly
	#han har sat 0 foran y divflux, så vi har ingen y directioon divflux langs I bottom og top
	dfly[0,:] = 0*(P[1,:]*v[1,:] - P[0,:]*v[0,:])/dy;
	dfly[ny-1,:] = 0*(P[ny-1,:]*v[ny-1,:] - P[ny-2,:]*v[ny-2,:])/dy;
	#I center:
	dfly[1:ny-1,:] = (P[2:ny,:]*v[2:ny,:]-P[0:ny-2,:]*v[0:ny-2,:])/(2*dy);

	#NU TAGER VI dflx
	#% Take cyclic differences on left and right boundaries

	dflx[:,0]=(P[:,1]*u[:,1]-P[:,nx-2]*u[:,nx-2])/(2*dx);
	dflx[:,nx-1]= dflx[:,0];

	#% take centered differences on interior points

	dflx[:,1:nx-1]= (P[:,2:nx]*u[:,2:nx]-P[:,0:nx-2]*u[:,0:nx-2])/(2*dx);
	
	return dflx,dfly


	#Interessant det der sker ved boundary periodic… vi springer en gridpoint over I centered differences, look it…

def Damping4(Dk4,nx,ny,U):
	numdif = np.zeros((ny,nx));          

	#Do smoothing in y space for 1st derivative zero at boundaries
	numdif[3:ny-4,:] = Dk4*(U[5:ny-2,:] -4*U[4:ny-3,:]+6*U[3:ny-4,:]
		-4*U[2:ny-5,:]+U[1:ny-6,:])
	numdif[2,:] = Dk4*(-3*U[1,:] +6*U[2,:]-4*U[3,:]+U[4,:])
	numdif[1,:] = Dk4*(2*U[1,:] -3*U[2,:] +U[3,:])
	numdif[ny-3,:] = Dk4*(-3*U[ny-2,:]+6*U[ny-3,:]-4*U[ny-3,:]+U[ny-4,:])
	numdif[ny-2,:] = Dk4*(2*U[ny-2,:] -3*U[ny-3,:] + U[ny-4,:])

	#%do smoothing in x space with periodicity
	numdif[:,2:nx-3] = numdif[:,2:nx-3]+Dk4*(U[:,4:nx-1] -4*U[:,3:nx-2]+6*U[:,2:nx-3]
		-4*U[:,1:nx-4]+U[:,0:nx-5])
	numdif[:,1] = numdif[:,1]+Dk4*(U[:,3] -4*U[:,2]+6*U[:,1]
		-4*U[:,0]+U[:,nx-1])
	numdif[:,0] = numdif[:,0]+Dk4*(U[:,2] -4*U[:,1]+6*U[:,0]
		-4*U[:,nx-1]+U[:,nx-2])
	numdif[:,nx-1] = numdif[:,0]
	numdif[:,nx-2] = numdif[:,nx-2]+Dk4*(U[:,1] -4*U[:,0]+6*U[:,nx-2]
		-4*U[:,nx-3]+U[:,nx-4])

	return numdif
	
def CalcVelocity(u,v,dxpsi,dypsi,psi,dx,dy):
	"""
	Calculate velocity (u,v) from stream function psi
	
	u = -dpsi/dy
	v = dpsi/dx
	
	Use finite difference
	
	More closely following Holton, he does something different on last x coordinate, dxpsi[:,-1]
	Also note, that we do
	dxpsi[:,0]=(psi[:,1]-psi[:,-2])/(2*dx)
	and NOT
	dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
	so, in the last psi, we jump to SECOND LAST coordinate, instead of last coordinate.
	so we skip the last column x = L, when calculating for x = 0
	
	"""
	#Calculate gradients for air velocity
	#forward on bottom? Backward on top?

	#dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
	#dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
	#dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference

	
	#x = 0, x = L
#	dxpsi[:,0]=(psi[:,1]-psi[:,-2])/(2*dx)
	#dxpsi[:,-1]= dxpsi[:,0]#(psi[:,0]-psi[:,-2])/(2*dx)
	#
	#Interior
	#dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx)  #centered difference

	
	dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
	dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
	dypsi[1:-1,:] = (psi[2:,:]-psi[0:-2,:])/(2*dy)  #Centered difference

	
	#x = 0, x = L
	dxpsi[:,0]=(psi[:,1]-psi[:,-2])/(2*dx)
	dxpsi[:,-1]= dxpsi[:,0]#(psi[:,0]-psi[:,-2])/(2*dx)
	
	#Interior
	dxpsi[:,1:-1]= (psi[:,2:]-psi[:,0:-2])/(2*dx)  #centered difference

	u = -dypsi
	v = dxpsi
	
	return u,v

def KineticEnergy(u,v):
	"""
	Calculate Domain integrated Kinetic Energy
	
	Should be conserved.
	Can also be done as |nablapsi|^2
	
	"""
	KE = np.sum(u**2+v**2)
	
	
	return KE
	
	
	
def Enstrophy(zeta):
	"""
	Calculate Domain integrated Enstrophy
	
	Domain Integrated Enstrophy should be conserved.
	calculated as,
	|zeta+f|^2
	"""
	Ens = np.sum(zeta*zeta)
	return Ens
	


def animate(i):

	global zetan,u,v,zeta,psi,psin,dxpsi,dypsi,zeta0,dx,dy,t,dt
	#leapfrog
	#zetan[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1]-2*dt*(\
	#	(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
	#	+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
	#	+beta*v[1:ny-1,1:nx-1])
	#zeta = zetan
	t+=dt
	zeta0 = zeta.copy()
	zeta = zetan.copy()
	
	numdif = Damping4(Av4,nx,ny,zeta0)
	
	#UpdateZetaLeapFrog(zetan,zeta0,zeta,u,v,beta)
	UpdateZetaLeapFrogHolton(zetan,zeta0,zeta,u,v,beta,numdif)
	zeta0 = zeta.copy()
	zeta = zetan.copy()
	
	
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
 	
	# #Calculate gradients for air velocity
	# #forward on bottom? Backward on top?
	#dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
	#dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
	#dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference


	#dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
	#dxpsi[:,-1]= dxpsi[:,0]#(psi[:,0]-psi[:,-2])/(2*dx)
	#dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx)  #centered difference

	#u = -dypsi
	#v = dxpsi
	
	u,v = CalcVelocity(u,v,dxpsi,dypsi,psi,dx,dy)
	
	dflx,dfly = divflux(zeta,u,v,dx,dy)
	
	KE = KineticEnergy(u,v)
	Ens = Enstrophy(zeta)
	print("KE = {}".format(KE))
	print("En = {}".format(Ens))

	ax.clear() #Hvorfor skal jeg bruge denne her? Den bruger de andre animations ikke
	C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
	#C = ax.contour(x/1000,y/1000,psi/100000,8,colors='black')
	ax.quiver(x/1000,y/1000,u,v)
	ax.set_title('Barotropic Vorticity Equation t = {}'.format(t))
	#ax.set_xlabel('x')
	#ax.set_ylabel('y')
	ax.set_xticks([-Lx*2/6,-Lx/6,0,Lx/6,Lx*2/6])
	ax.set_yticks([-Ly*2/6,-Ly/6,0,Ly/6,Ly*2/6])
	ax.set_xlabel("x/km")
	ax.set_ylabel("y/km")
	#plt.xticks([-Lx/6,Lx/6],[-Lx/6,Lx/6])
	#plt.xticks([-Lx/6,Lx/6])
	#ax.set_xticks([0],('hey'))
	plt.clabel(C,inline=1,fontsize=10,fmt="%1.1f")
	#print(zeta)

	if i == 4:
		fig.savefig('BVE.png', bbox_inches='tight')


if __name__ == "__main__":
	
	
	
	Lx = 6000 #km
	Ly = 6000 #km
	nx = 65
	ny = 65
	pi = 3.141592
	#x = np.linspace(0,Lx,nx)
	#y = np.linspace(0,Ly,ny)
	X = np.linspace(-Lx/2,Lx/2,nx)
	Y = np.linspace(-Ly/2,Ly/2,ny)
	x,y = np.meshgrid(X*1000,Y*1000)
	k = 2*pi/(Lx*1000)
	m = pi/(Ly*1000)
	dx = 1000*Lx/(nx-1)
	dy = 1000*Ly/(ny-1)
	U0 = 20				#zonal wind
	beta = 1.62*10**(-11)	#he set 0 infront?
	Av4 = 10**(-6)
	A = 10**(-4)
	#initial vorticity and streamfunction

	zeta0 = np.array(A*np.exp(-2*(k**2*x**2+m**2*y**2)),dtype=np.float64)
	zeta = zeta0
	zetan = zeta0
	#time integration parameters					#hours
	time_end = 3*3600 				#second
	dt = 100

	psi = np.zeros((ny,nx),dtype=np.float64)
	dypsi = np.zeros((ny,nx),dtype=np.float64)
	dxpsi = np.zeros((ny,nx),dtype=np.float64)
	u = np.zeros((ny,nx),dtype=np.float64)
	v = np.zeros((ny,nx),dtype=np.float64)
	dfly = np.zeros((ny,nx),dtype=np.float64)
	dflx = np.zeros((ny,nx),dtype=np.float64)

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
	# dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
	# dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
	# dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference


	# dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
	# dxpsi[:,-1]= dxpsi[:,0]#(psi[:,0]-psi[:,-2])/(2*dx)
	# dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx)  #centered difference

	# u = -dypsi
	# v = dxpsi

	u,v = CalcVelocity(u,v,dxpsi,dypsi,psi,dx,dy)

	dflx,dfly = divflux(zeta,u,v,dx,dy)

	#Forward time difference
	FirstStepZeta(zeta0,zeta,u,v,beta)
	# zetan[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1]-dt*(\
		# (u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
		# +(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
		# +beta*v[1:ny-1,1:nx-1])
	# zeta = zetan




	t=0
	t += dt
	#plt.hold(True)
	fig = plt.figure()
	ax = plt.gca()
	#ax.contour(x,y,zeta,colors='black')
	#ax.quiver(x,y,u,v)
	ax.set_title('Barotropic Vorticity Equation')
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	#ax.set_xticks([0])
	#ax.set_yticks([0])
	#ax.set_xticks([Lx/3,0],["hey","loL"])
	#ax.set_yticks([Ly/3,2*Ly/3],[1,1])

	
	
	
	#From animation function,
	for timestep in range(100):
		t+=dt
		zeta0 = zeta.copy()
		zeta = zetan.copy()
		
		numdif = Damping4(Av4,nx,ny,zeta0)
		
		#UpdateZetaLeapFrog(zetan,zeta0,zeta,u,v,beta)
		UpdateZetaLeapFrogHolton(zetan,zeta0,zeta,u,v,beta,numdif)
		zeta0 = zeta.copy()
		zeta = zetan.copy()
		
		
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
		
		# #Calculate gradients for air velocity
		# #forward on bottom? Backward on top?
		#dypsi[0,:] = (psi[1,:] - psi[0,:])/dy
		#dypsi[-1,:] = (psi[-1,:] - psi[-2,:])/dy
		#dypsi[1:-2,:] = (psi[2:-1,:]-psi[0:-3,:])/(2*dy)  #Centered difference


		#dxpsi[:,0]=(psi[:,1]-psi[:,-1])/(2*dx)
		#dxpsi[:,-1]= dxpsi[:,0]#(psi[:,0]-psi[:,-2])/(2*dx)
		#dxpsi[:,1:-2]= (psi[:,2:-1]-psi[:,0:-3])/(2*dx)  #centered difference

		#u = -dypsi
		#v = dxpsi
		
		u,v = CalcVelocity(u,v,dxpsi,dypsi,psi,dx,dy)
		
		dflx,dfly = divflux(zeta,u,v,dx,dy)
		
		KE = KineticEnergy(u,v)
		Ens = Enstrophy(zeta)
		print("KE = {}".format(KE))
		print("En = {}".format(Ens))

		ax.clear() #Hvorfor skal jeg bruge denne her? Den bruger de andre animations ikke
		C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
		#C = ax.contour(x/1000,y/1000,psi/100000,8,colors='black')
		ax.quiver(x/1000,y/1000,u,v)
		ax.set_title('Barotropic Vorticity Equation t = {}'.format(t))
		#ax.set_xlabel('x')
		#ax.set_ylabel('y')
		ax.set_xticks([-Lx*2/6,-Lx/6,0,Lx/6,Lx*2/6])
		ax.set_yticks([-Ly*2/6,-Ly/6,0,Ly/6,Ly*2/6])
		ax.set_xlabel("x/km")
		ax.set_ylabel("y/km")
		#plt.xticks([-Lx/6,Lx/6],[-Lx/6,Lx/6])
		#plt.xticks([-Lx/6,Lx/6])
		#ax.set_xticks([0],('hey'))
		plt.clabel(C,inline=1,fontsize=10,fmt="%1.1f")
		#print(zeta)

	#===========================
	#USe animation function
	# anim = animation.FuncAnimation(
	# fig, 
	# animate, 
	# #frames=5, 
	# interval=0.5,
	# blit=False #blit=False default, 
	# )
	# plt.show()