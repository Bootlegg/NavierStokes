
def PsiSolvePoissonJacobiChannel(psi,zeta,dx,dy,nx,ny,epstol,Ngc):
#Compute streamfunction from vorticity
	#ζ_(i,j)^1=(((ψ_(i+1,j)-2ψ_(i,j)+ψ_(i-1,j) ))/dx^2 +((ψ_(i,j+1)-2ψ_(i,j)+ψ_(i,j-1) ))/dy^2 )

	#solve psi from lec12.pdf
	psin = np.zeros((ny,nx),dtype=np.float64)
	psin[:,:] = psi.copy()
	#psin = np.zeros((ny,nx),dtype=np.float64)
	
	
	#error betweeen iterated solution and RHS at each gridpoint
	eps = np.zeros((ny,nx),dtype=np.float64)
	
	#Maybe not needed for some implementations
	#Nabla^2psi
	nabla2psi = np.zeros((ny,nx),dtype=np.float64)

	
	##JACOBI METHOD
	dtau = (1/4)*(0.5*dx**2+0.5*dy**2)
	nit = 0
	#while True:
	for r in range(4000): #pseudo-time
	
		#update iteration counter
		nit +=1
		
		
		#Make copy of current iterate
		#We should probably do a copy, because, otherwise, we are doing intermediate changes...
		#i.e we are changing values that will be used later, that shouldn't....
		#So, that is, to really follow the algorithm, we need to copy the values...
		#At least, most likely... if it was a double for loop, then we would have to....
		#BUT, MAYBE, we don't have to... maybe numpy is smart enough to figure it out
		psin[:,:] = psi.copy()
		
		
		
		#Calculate intermediate array
		#1/dx**2(psi[1:ny-1,2:nx]+psi[1:ny-1,0:nx-2]-2*psi[1:ny-1,1:nx-1]) +
		#1/dy**2(psi[2:ny,1:nx-1]+psi[0:ny-2,1:nx-1]-2*psi[1:ny-1,1:nx-1])
		#nabla2psi[1+Ngc:ny-1-Ngc,Ngc+1:nx-1-Ngc]= ((1/dx**2)*(psi[1+Ngc:ny-1-Ngc,2+Ngc:nx-Ngc]
		#						+psi[1+Ngc:ny-1-Ngc,0+Ngc:nx-2-Ngc]-2*psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc])
		#						+(1/dy**2)*(psi[2+Ngc:ny-Ngc,1+Ngc:nx-1-Ngc]
		#						+psi[0+Ngc:ny-2-Ngc,1+Ngc:nx-1-Ngc]-2*psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]))
		
		#Calculate error
		#eps[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc] = nabla2psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]-zeta[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]
		
		
		#Update iteration of streamfunction
		#psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc] = psin[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]+eps[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]/(2/dx**2+2/dy**2)
		
		
		#Interior points
		#psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(
		#		+(psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2
		#		+(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2
		#		-zeta[1:ny-1,1:nx-1])
				
		psi[1:ny-1,1:nx-1] = 0.25*(psin[1:ny-1,0:nx-2]+psin[1:ny-1,2:nx]+psin[2:ny,1:nx-1]+psin[0:ny-2,1:nx-1])-(dx**2/4)*zeta[1:ny-1,1:nx-1]
		
		#print((dx**2/4)*zeta[1:ny-1,1:nx-1])
		#https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
		#psi[1:ny-1,1:nx-1] = ((1/dx**2)*(psi[1:ny-1,2:nx]+psi[1:ny-1,0:nx-2]-2*psi[1:ny-1,1:nx-1])
		#						+(1/dy**2)*(psi[2:ny,1:nx-1]+psi[0:ny-2,1:nx-1]-2*psi[1:ny-1,1:nx-1])
		#						-zeta[1:ny-1,1:nx-1])/(2/dx**2+2/dy**2)
		
		#Bottom BCs
		#Dirichlet BCs, psi = 0
		#Not strictly needed, since we only update interior points anyway
		#psi[0,:] = 0
		#psi[ny-1,:] = 0
		
		#Jacobi interior
		#psi[1:ny-1,1:nx-1] = ((1.0/4.0)*(
		#						+psin[1:ny-1,2:nx]
		#						+psin[1:ny-1,0:nx-2]
		#						+psin[2:ny,1:nx-1]
		#						+psin[0:ny-2,1:nx-1])
		#					-(dx**2/4)*zeta[1:ny-1,1:nx-1])
		
		
		
		#x = 0 boundary
		#Tager et note ud af Holtons bog, og her, skipper vi psin[1:ny-1,-1] og bruger psin[1:ny-1,-2]
		#psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
		#		+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-2])/dx**2
		#		+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
		#		-zeta[1:ny-1,0])
				
		#psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
		#		+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-1])/dx**2
		#		+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
		#		-zeta[1:ny-1,0])		
		
		####
		psi[1:ny-1,0] = 0.25*(psin[1:ny-1,1]+psin[1:ny-1,nx-1]+psin[2:ny,0]+psin[0:ny-2,0])-(dx**2/4)*zeta[1:ny-1,0]
		#periodic BCs
		#psi[1:ny-1,nx-1] = psi[1:ny-1,0]
		psi[1:ny-1,nx-1] = 0.25*(psin[1:ny-1,0]+psin[1:ny-1,nx-2]+psin[2:ny,nx-1]+psin[0:ny-2,nx-1])-(dx**2/4)*zeta[1:ny-1,nx-1]
		
		
		#print(psi)
		
		
		#https://www3.nd.edu/~gtryggva/CFD-Course2010/2010-Lecture-11.pdf
		#for boundaries use (1/3) if eg f = 0 because then there's only 3 terms
		
		#I am misisng some boundary conditions here...
		#I think, on TOP of Dirichtlet conditions psi = 0 on y = 0 and y = ny
		#we should ALSO have Neumann conditions, maybe?
		#dpsi/dn = 0
		#differentiate psi wrt. the normal vector, should be 0...
		#psi = 0 on BC gives, no-cross flow,
		##or is it the other way around...
		#And dpsi/dn = 0 would give no-slip condition, right?
		#
		
		
		
		#Old
		#psi[1:ny-1,0] = psin[1:ny-1,0]+dtau*(
		#		+(psin[1:ny-1,1]-2*psin[1:ny-1,0]+psin[1:ny-1,-1])/dx**2
		#		+(psin[2:ny,0]-2*psin[1:ny-1,0]+psin[0:ny-2,0])/dy**2
		#		-zeta[1:ny-1,0])		
		
		#x = L boundary
		#Enten så gør den periodic sådan her
		#psi[1:ny-1,-1] = psi[1:ny-1,0]
		#eller sådan her, Holton
		#psi[:,-1] = psi[:,0]
		
		#Eller gør den periodic sådan her
		#psi[1:ny-1,-1] = psin[1:ny-1,-1]+dtau*(
		#+(psin[1:ny-1,0]-2*psin[1:ny-1,-1]+psin[1:ny-1,-2])/dx**2
		#+(psin[2:ny,-1]-2*psin[1:ny-1,-1]+psin[0:ny-2,-1])/dy**2
		#-zeta[1:ny-1,-1])
		
		
		
		
		#Psi no-flow across boundary GHOST CELLS
		#boundary at x = L, i'm gonna set equal to x = 0...
		#Boundary for psi, maybe i should remove these!
		#psi[:,-1] = 0 #right boundary
		#psi[:,0] = 0	#left boundary
		
		
		#What if i don't set it to 0? afaik it should be constant along the edges, but maybe not 0...
		#Det har i hvert fald bestemt en effekt om man sætter det til 0 eller ej..
		#Probably not needed to put it to 0, if I'm just careful the other places,
		#But, to follow the algorithm completely, I here set the y boundarys upper and lower to 0
		#psi[0,:] = 0  #south boundary
		#psi[-1,:] = 0 #north boundary
		


		#eps[1:ny-1,1:nx-1] = nabla2psi[1:ny-1,1:nx-1]-zeta[1:ny-1,1:nx-1]
		#print(eps)
		#if (np.abs(eps) < epstol).all():
		#	print("Break Poisson algorithm at", nit)
		#	break
		
		
	
	#I think here, we should return psi, instead of doing it implicitly....
	#so, change the way the method works...
	return psi


def gradHolton(u,v,psi,dx,dy,nx,ny):
		# % Take forward differences on top and bottom edges
	# dely(1,:) = (P(2,:) - P(1,:))/dy;
	# dely(n,:) = (P(n,:) - P(n-1,:))/dy;

	# % Take centered differences on interior points

	# dely(2:n-1,:) = (P(3:n,:)-P(1:n-2,:))/(2*dy);
	# % compute x component of gradient 

	# % Take cyclic differences on left and right boundaries

	# delx(:,1)=(P(:,2)-P(:,p-1))/(2*dx);
	# delx(:,p)= delx(:,1);

	# %DAN u=-delx
	# %DAN v=dely

	# % take centered differences on interior points

	# delx(:,2:p-1)= (P(:,3:p)-P(:,1:p-2))/(2*dx);
	
	
	#Forward differences top and bottom
	v[0,:] = (psi[1,:]-psi[0,:])/dy
	v[ny-1,:] = (psi[ny-1,:]-psi[ny-2,:])/dy
	
	#Centered differences interior
	u[:,1:nx-1] = -(psi[:,2:nx]-psi[:,0:nx-2])/(2*dx)
	v[1:ny-1] = (psi[2:ny,:]-psi[0:ny-2,:])/(2*dy)
	
	#Cyclic differences left and right
	#First take centered differences at x = 0, using the first and last point
	u[:,0] = -(psi[:,1]-psi[:,nx-1])/(2*dx)
	u[:,nx-1] = u[:,0]
	
	#If I want symmetry at the periodic boundary...
	#then, there's a difference, if I do calculate it from nx, and copy that to 0, vs calculate from 0, and copy to nx
	#i think that must be asymmetry
	return u,v
	
	
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d
import matplotlib

#np.set_printoptions(precision=2)

Backend = matplotlib.get_backend()
print(Backend)#We use TkAgg
#matplotlib.use("wx") #No wx
matplotlib.use("TkAgg")
from pylab import *
thismanager = get_current_fig_manager()
#thismanager.window.SetPosition((500, 0))
thismanager.window.wm_geometry("-1500+0")



Lx = 6000 #km
Ly = 3000 #km
nx = 65
ny = 65
pi = 3.141592
#x = np.linspace(0,Lx,nx)
#y = np.linspace(0,Ly,ny)
#X = np.linspace(-Lx/2,Lx/2,nx)
X = np.linspace(0,Lx,nx)
Y = np.linspace(-Ly,Ly,ny)
x,y = np.meshgrid(X*1000,Y*1000)
#k = 2*pi/(Lx*1000)
#m = pi/(Ly*1000)
dx = 1000*Lx/(nx-1)
dy = 2*1000*Ly/(ny-1)
U0 = 20				#zonal wind

R_earth = 6381*10**3 #meters
Omega_earth = 1/(24*60*60)
phi0 = np.pi/4 #45degrees
#beta = df/dy = 2*Omega_Earth*cos(phi0)/R_earth
beta = 2*Omega_earth*np.cos(phi0)/R_earth

print("coriolis beta ", beta)
zeta0 = np.zeros((ny,nx),dtype=np.float64)
zeta = zeta0.copy()
zetanew = zeta0.copy()
psi = np.zeros((ny,nx),dtype=np.float64)
dypsi = np.zeros((ny,nx),dtype=np.float64)
dxpsi = np.zeros((ny,nx),dtype=np.float64)
u = np.zeros((ny,nx),dtype=np.float64)
v = np.zeros((ny,nx),dtype=np.float64)
epstol = 1e-8
Ngc = 0

k1 = 2*pi/(Lx*1000)
k2 = 1*pi/(Ly*1000)
w = -k1/(k1**2+k2**2)
psi[:,:] = np.cos(k1*x)*np.sin(k2*y)




psiexact = np.zeros((ny,nx),dtype=np.float64)
psiexact[:,:] = np.cos(k1*x-w*0)*np.sin(k2*y)
#https://www.wolframalpha.com/input/?i=d%5E2%2Fdx%5E2cos%28ax-wt%29sin%28by%29%2Bd%5E2%2Fdy%5E2cos%28ax-wt%29sin%28by%29
zetaexact = -(k1**2+k2**2)*(np.cos(k1*x-w*0)*np.sin(k2*y))
#https://www.wolframalpha.com/input/?i=d%2Fdx+cos%28ax-wt%29sin%28by%29
uexact = -k1*np.sin(k2*y)*np.sin(k1*x-w*0)

#Oh actually
#I could ALSO calculate the velocities,
#And even the vorticity, from this expression!!
#So, if i can get zeta... would be very nice...
#Then, from the REAL zeta values, try to reconstruct psi....
#

#
#Btw, does the Rossby solution make physical sense? Because psi = gf/h or something... will this expression EVER be negative??


#WE HAVE A PROBLEM YES. The boundaries are NOT psi = 0, wrong setup!!
#OR, at least, it doesn't fulfil psi = 0 like this...
#NOW... try to give the poisson solver initial psi = 0, NOT psi = exact



print("psi[0,0:5] bottom",psiexact[0,0:5])
print("psi[ny-1,0:5] top",psiexact[ny-1,0:5])

u,v = gradHolton(u,v,psi,dx,dy,nx,ny)

# #Rescale velocity
# #done e.g here https://dspace.library.uvic.ca/bitstream/handle/1828/1218/Final%20THESIS.pdf?sequence=1
# #let's find max velocity
# vspeed = np.sqrt(u*u+v*v)
# vspeedmax = np.max(vspeed) #around 1e-6m/s
# vspeedScale = 5 #5m/s
# #k*vspeedmax = vspeedScale
# #k = vspeedScale/vspeedmax
# kspeedscale = vspeedScale/vspeedmax
# u = kspeedscale*u
# v = kspeedscale*v
#We can calculate zeta if given initial psi, just forward method
#Right now, doesn't calculate edges,
#so just by default, zeta = 0 on edges
for i in range(1,nx-1):
	for j in range(1,ny-1):
		#u[i,j] = -(psi[i+1,j]-psi[i-1,j])/(2*dy)
		#v[i,j] = (psi[i,j+1]-psi[i,j-1])/(2*dx)
		zeta0[j,i] = (v[j,i+1]-v[j,i-1])/(2*dx)-(u[j+1,i]-u[j-1,i])/(2*dy)

zeta0[1:ny-1,0] = 	(v[1:ny-1,1]-v[1:ny-1,nx-1])/(2*dx)-(u[2:ny,0]-u[0:ny-2,0])/(2*dy)
zeta0[1:ny-1,nx-1] = zeta0[1:ny-1,0]
#zeta = 0 on BCs, i think, for FULL SLIP...
#But, maybe the rossby wave packet defined, is NOT good with zeta = 0 on BCs,
#maybe those initial conditions do not satisfy the BCs...

#psi = PsiSolvePoissonJacobiChannel(psi,zeta,dx,dy,nx,ny,epstol,Ngc,)




#initial streamfunction
psizeros = np.zeros((ny,nx),dtype=np.float64)
#psi0 = PsiSolvePoissonJacobiChannel(psizeros,zeta0,dx,dy,nx,ny,epstol,Ngc)
psi0 = PsiSolvePoissonJacobiChannel(psizeros,zetaexact,dx,dy,nx,ny,epstol,Ngc)
print("zeta0 from fundamental forward definition")
print(zeta0)

print("zetaexact from Rossby wave")
print(zetaexact)

print("psi0 from poisson solver")
print(psi0)

print("psizeros input from poisson solver")
print(psizeros)
#u0,v0 = CalcVelocity(psi0,u,v,dx,dy)
u0,v0 = gradHolton(u,v,psi,dx,dy,nx,ny)

print("dx**2",dx**2)


print("dx**2*zetaexact",dx**2*zetaexact)
#plot initial poisson solution of nabla^2psi = zeta
fig00 = plt.figure() #If i do pcolor, then no need for 3d projection
thismanager = get_current_fig_manager()
thismanager.window.wm_geometry("-1500+0")
#ax00 = fig00.gca(projection='3d')
ax00 = fig00.gca()
#ax00.plot_surface(x, y, psi0)#, rstride=3, cstride=3, color='black')
C = ax00.contour(x/1000,y/1000,psi0,4,colors='black')
ax00.set_title('Initial Psi from Poisson solver nabla2psi = xi')
ax00.quiver(x/1000,y/1000,u,v)
ax00.set_xlabel('x')
ax00.set_ylabel('y')
plt.show()


fig003 = plt.figure(3) #If i do pcolor, then no need for 3d projection
thismanager = get_current_fig_manager()
thismanager.window.wm_geometry("-1500+0")
#ax003 = fig003.gca(projection='3d')
ax003 = fig003.gca()
#ax003.plot_surface(x, y, psi)#, rstride=3, cstride=3, color='black')
C = ax003.contour(x/1000,y/1000,psiexact,4,colors='black')
ax003.set_title('Initial Psi Rossby wave packet exact expression')
ax003.quiver(x/1000,y/1000,u,v)
ax003.set_xlabel('x')
ax003.set_ylabel('y')
plt.show()



# #plot initial conditions
# fig0 = plt.figure(1)
# thismanager = get_current_fig_manager()
# thismanager.window.wm_geometry("-1500+0")
# ax0 = plt.gca()
# ax0.set_title('Barotropic Vorticity Equation Initial configuraton from Poisson')
# ax0.set_xlabel('x')
# ax0.set_ylabel('y')
# #ax.quiver(x/1000,y/1000,u,v)
# #C = ax0.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
# C = ax0.contour(x/1000,y/1000,zeta0,8,colors='black')
# C = ax0.contour(x/1000,y/1000,psi0,4,colors='black')
# plt.clabel(C, fontsize=10, inline=1,fmt = '%1.0f')
# ax0.quiver(x/1000,y/1000,u0,v0)
# plt.savefig('BVEboundingboxinitial.png')
# plt.show()



# fig = plt.figure(2)
# thismanager = get_current_fig_manager()
# thismanager.window.wm_geometry("-1500+0")
# ax = plt.gca()
# ax.set_title('Barotropic Vorticity Equation final config')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# #ax.quiver(x/1000,y/1000,u,v)
# #C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
# #Czeta = ax.contour(x/1000,y/1000,zeta,6,colors='black')
# Cpsi = ax.contour(x/1000,y/1000,psi,4,colors='black')
# #plt.clabel(Czeta, fontsize=10, inline=1)#,fmt = '%1.0f')
# plt.clabel(Cpsi, fontsize=10, inline=1)#,fmt = '%1.0f')
# ax.quiver(x/1000,y/1000,u,v)
# ax.legend()
# plt.savefig('BVEboundingboxfinal.png')
# plt.show()




def FFTPoisson(nx,ny,Lx,dy,zeta):
	# function psi = stream1(Nxl,Nyl,Lx,dy,zeta)
	# % Invert vorticity for streamfunction.
	# % Nx is number of grid points in x, Ny number in y.
	# % Lx and Ly are domain lengths in x and y
	# % zeta is the vorticity specified on the Nx by Ny grid.
	# % Use cyclic conditions in x with FFT, Finite Difference in y.
	# psiT=zeros(size(zeta));
	# psi=psiT;
	# k = [0:Nxl-1 ] * (2*pi/(Lx*1000));  % Fourier wavenumber operators
	# K2 = k.^2;
	# dy2 = dy^2;
	# K2dy2 = K2*dy2;
	# % take fft of vorticity in x direction
	# fzeta = fft(zeta,Nxl,2);
	# % now invert for streamfunction
	# for s = 1:Nxl/2;
		# %   Coefficients for tridiagonal solver
		# A = -(K2dy2(s) +2);
		# B(2:Nyl) = dy2*fzeta(2:Nyl,s);      % forcing term
		# B(2) = B(2);                        % boundary condition for y = 0
		# e = ones(Nyl,1);
		# % Define the tridiagonal matrix for finite difference equation
		# M = spdiags([e A*e e], -1:1, Nyl-1,Nyl-1);
		# % solve M*psiT = B for psiT by matrix inversion
		# psiT(2:Nyl,s) = M\B(2:Nyl).';
		# psiT(1,s) = 0;
	# end
	# psi(:,1:Nxl) = 2*real(ifft(psiT,Nxl,2)); % grid point streamfunction
	# psi(:,Nxl+1) = psi(:,1); 
	# %DAN Nxl = Nx-1
	# %psi(:,Nx-1+1) = psi(:,1);
	# %psi(:,Nx) = psi(:,1);

	psiT = np.zeros_like(zeta)
	nxl = nx-1
	nyl = ny-1
	k = np.array([i for i in range(nx-1)])*(2*np.pi/(Lx*1000))
	K2 = k**2
	dy2 = dy**2
	K2dy2 = K2*dy2
	#take fft of zeta in x direction
	fzeta = np.fft.fft(zeta, n = nxl, axis=1, norm=None) #fft.fft(a, n=None, axis=-1, norm=None)
	#print("fzeta shape")
	#print(fzeta.shape)
	for s in range(0,int(nxl/2)):	
		A = -(K2dy2[s]+2)
		B = np.zeros(nyl)
		#print(B.shape)
		#print(B)
		B[1:nyl] = dy2*fzeta[1:nyl,s]
		#e = np.ones((nyl,1))
		#MAKE M Matrix, for
		#M*psiT = B
		#Solve for psiT
		M = np.eye(nyl-1,k=-1)+np.eye(nyl-1)*A+np.eye(nyl-1,k=1)
		#print(M)
		Minv = np.linalg.inv(M)
		#print("Minv*B")
		#print(Minv*B[1:nyl])
		psiT[1:nyl,s] = np.dot(Minv,B[1:nyl])
		psiT[0,s] = 0 #probably boundary condition, on bottom
	psi[:,0:nxl] = 2*np.real(np.fft.ifft(psiT,nxl,1))#2*np.real(np.fft.ifft(psiT,nxl,2))
	psi[:,nx-1] = psi[:,0]
	
	return psi
	
psiFFT = FFTPoisson(nx,ny,Lx,dy,zetaexact)

#plot initial poisson solution of nabla^2psi = zeta
fig00 = plt.figure() #If i do pcolor, then no need for 3d projection
thismanager = get_current_fig_manager()
thismanager.window.wm_geometry("-1500+0")
#ax00 = fig00.gca(projection='3d')
ax00 = fig00.gca()
#ax00.plot_surface(x, y, psi0)#, rstride=3, cstride=3, color='black')
C = ax00.contour(x/1000,y/1000,psiFFT,4,colors='black')
ax00.set_title('Initial Psi from Poisson solver nabla2psi = xi FFT')
ax00.quiver(x/1000,y/1000,u,v)
ax00.set_xlabel('x')
ax00.set_ylabel('y')
plt.show()