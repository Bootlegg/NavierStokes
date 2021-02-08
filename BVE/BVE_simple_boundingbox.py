"""
Barotropic Vorticity Equation
Holton pg 464-465
dzeta/dt = -F(x,y,t)
F(x,y,t) = d/dx...+d/dy....+beta*v_psi
Hence, we need a minus sign infront of cartesian beta plane term
In the end, that is, because F will get a negative sign itself, 

Overall algorithm
Initiate Zeta
1) Compute Psi from Zeta
2) Compute u,v from Psi
3) Update Zeta
4) Repeat

Holton p 119
eq 4.41 is called BVE


Holton pg 42
for synoptic scale motions, convenient to consider disturbance at 45 latitude, where f0 = 2*Omega*sin(45) = 10^(-4)s
Holton pg 160
Midlatitude beta.plane approximation
we expand the Coriolis parameter in a Taylor series about a latitude phi0
retaining only first two terms yield
f = f0 + beta*y
where beta = df/dy = 2*Omega*cos(phi0)/a
and y = 0 at phi0



https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
600-500hPa is pretty good for BVE, that part of the atmosphere has low divergence (BVE assumes 0 divergence)



https://apps.dtic.mil/dtic/tr/fulltext/u2/a065117.pdf
BVE can be used to model atmospheric planetary waves in the troposphere
Her har de nabla^2psi = eta-f
#Nok fordi eta-f = zeta


https://pdfs.semanticscholar.org/c937/d58642e376a098b2b319783cd121c7fbbfe9.pdf
Here,
dzeta/dt = -J(psi,nabla^2psi+f)+F+D(psi)

nabla^2 is Laplace operator
J is Jacobi operator
F is forcing function
D(psi) is diffusion term, usually -kappa*nabla^2psi (bottom-friction case)



https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
streamfunction psi = gh/f


Rossby-Haurwitz waves, This is when,
#The STREAMFUNCTION psi, is equal to a spherical harmonic mode,
So i think this is when doing it on a sphere


https://dspace.library.uvic.ca/bitstream/handle/1828/1218/Final%20THESIS.pdf?sequence=1
Solution to linearized equations, represents a wave travelling westward
but, only an eastward moving wave, called Kelvin wave, is retained



#Time increment,
https://maths.ucd.ie/~plynch/eniac/CFvN-1950.pdf
Here, they talk about dt = 15min or less, something to do with gravitational waves of 300m/s
Ahh and then of course, since it's a forecast, let's say we integrate 24h forward...
so with 15 minutes, it's around 100 cycles


Domain size,
https://maths.ucd.ie/~plynch/eniac/CFvN-1950.pdf
If you want to forecast some area, you should actually integrate a domain slightly larger than that area
because what happens outside the domain will influence at the boundaries...
But the speed at which outside influence travels it not that great, so you only need to integrate a slightly larger domain




http://www.m-hikari.com/ams/ams-2013/ams-49-52-2013/koomsubsiriAMS49-52-2013.pdf
Here, centered differences for interior points
At top and bottom edges, they have, forward differences
and left and right edges, they have cyclic conditions

domain is X = 6000km, Y = 3000km,
d=dx=dy=93,750m
dt = 900s, integrate up to 10days,

initial vorticity:
zeta0 = 10^-4*exp(-2(k^2x^2+m^2y^2))
where 
k=2pi/(6*10^6)
m=pi/(3*10^6)

they ALSO plot the psi streamfunction, to see the gradient etc



http://empslocal.ex.ac.uk/people/staff/dbs202/cag/courses/MTMW14/notes2006.pdf
psi = g/f * ø
where ø is geopotential
that is probably the geostrophic relation
soo... i should be able to also calculate geopotential ø = psi*f/g
also
zeta = dv/dx - du/dy
and since u = -dpsi/dx and v = dpsi/dy
then you get

zeta = d^2psi/dx^2 + d^2psi/dy^2
so then
nabla^2psi = zeta

you start with 500hPa geopoential ø

then, 
BCs
free-slip on north and south channel, hence v = 0 on north and south, no flow in or out of boundary
also, probably, du/dy = 0 no shear on the boundaries

this means, constant stream function on northern and southern boundary, so e.g psi = 0
and on northern and southern boundary zeta = 0

zeta = dv/dx - du/dy
zeta = 0 - 0
zeta = 0 on northern and southern boundary




solve poisson zeta = nabla^2psi
we need 2 boundary conditions, we have them, they are e.g psi = 0 on both
But I guess the problem is, it could be psi = A on northern border and psi = B on sourthen body
and apparantly we need the correct ones from the initial geopotential field 500hPa



BCs:
For bounding box, we need to set psi = 0 on ALL sides
Probably also some conditions on zeta

barotropic_model.m and grad.m from Holton scripts!
in grad.m, they ALSO use forward differences, for velocity, i.e differentiation psi...
so Holton, uses Full-Slip boundary conditions on north and south!
and it's implemented, with forward differences for u = -dpsi/dy


https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
#no flow across boundary, on south and northern, means v = 0
#then, you need psi constant in x direction....
#because v = dpsi/dx...
#so, along south and northern channel, i will have eg psi = 0 and v = 0
#at least... MAYBE also u = 0...
tangential, you can either have full slip, or no-slip
if no-slip, then ALSO u = 0
then, u0 = uNY+1 = 0, i'm pretty sure...
i.e the buffer zones, are equal to 0
BUT.. zeta is NOT equal to 0...
#Does this come automatically, from my jacobian atm??? or??
zeta(0,NY+1) = 2(psi(1,NY+1)-psi(0,NY))/dy^2





We need to make sure to calculate the BCs of psi first, because they are later used to calculate velocities...
#So make sure psi etc first satisfies periodic or hard wall boundaries etc



Robert-Asselin filters
https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf



Domains:
I should probably make buffer/ghost gridcells, +-1 on each side at least,
and as CFvN-1950 noted, you usually simulate slightly larger than what you are interested in
so both effects, grid is slightly larger than area of interest, AND we add additional cells to help with BCs





Iniital Conditions:
http://www.math.ualberta.ca/ijnam/Volume-10-2013/No-3-13/2013-03-04.pdf
Ahh, jeg tror man skal vel egentlig have initial velocity field....
faktisk
zeta = du/dy - dv/dx 
den giver jo ikke nogen unique requirement for speeds, altså, speeds på 2m/s og speeds på 400m/s,
de kan vidst faktisk give samme relative vorticity zeta
and btw, they set initial conditions in term of STREAM function.... they put Rossby Wave on STREAMFUNCTION
NOT as vorticity zeta
Yeah, også Holton, Rossby-Haurwits R-H waves de er vidst defined ud fra streamfunction psi, og IKKE zeta!


Ahh ja, det er også det som Holton siger right
Man har jo 500hPa geopotential field, og det er jo mere directly related til STREAMFUNCTION psi!



btw, det er vel som shallow water equations,
horizontal velocities er same HELE vejen op gennem atmosfære, Kinda unrealistic


barotropic_model.m from Holton, uses vorticity as initial condition, and THEN computes psi streamfunction
So I think even Holton file, probably skips some steps and doesn't QUITE follow what is more real


Dan Krog
"""

def FFTPoisson(nx,ny,Lx,dy,zeta,psi):
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
	k = np.array([i for i in range(nx-2)])*(2*np.pi/(Lx*1000)) #nx-2 vs nx-1(virkede vidst)
	K2 = k**2
	dy2 = dy**2
	K2dy2 = K2*dy2
	#take fft of zeta in x direction
	#fzeta = np.fft.fft(zeta, n = nxl, axis=1, norm=None) #fft.fft(a, n=None, axis=-1, norm=None)
	fzeta = scipy.fft(zeta, n = nxl, axis=1, norm=None) #fft.fft(a, n=None, axis=-1, norm=None)
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
	psi[:,0:nxl] = 2*np.real(np.fft.ifft(psiT,n=nxl,axis=1))#2*np.real(np.fft.ifft(psiT,nxl,2))
	psi[:,nx-1] = psi[:,0]
	
	return psi
def poisson_periodic(p,b, dx, dy):
	pn = np.empty_like(p)
	
	for q in range(500):
		pn = p.copy()
		p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
					(pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
					(2 * (dx**2 + dy**2)) -
					dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

		# Periodic BC Pressure @ x = 2
		p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 +
					(pn[2:, -1] + pn[0:-2, -1]) * dx**2) /
					(2 * (dx**2 + dy**2)) -
					dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])

		# Periodic BC Pressure @ x = 0
		p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 +
					(pn[2:, 0] + pn[0:-2, 0]) * dx**2) /
					(2 * (dx**2 + dy**2)) -
					dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])
        
		# Wall boundary conditions, pressure
		#p[-1, :] =p[-2, :]  # dp/dy = 0 at y = 2
		#p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
		
		
		#No in-flow, v = dpsi/dx = 0
		p[-1,:] = 0
		p[0,:] = 0
    
	return p

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
	
	
	delx = np.zeros_like(v)
	dely = np.zeros_like(u)
	
	#Forward differences top and bottom
	dely[0,:] = -(psi[1,:]-psi[0,:])/dy
	dely[ny-1,:] = -(psi[ny-1,:]-psi[ny-2,:])/dy
	
	#Centered differences interior
	delx[:,1:nx-1] = (psi[:,2:nx]-psi[:,0:nx-2])/(2*dx)
	dely[1:ny-1] = -(psi[2:ny,:]-psi[0:ny-2,:])/(2*dy)
	
	#Cyclic differences left and right
	#First take centered differences at x = 0, using the first and last point
	delx[:,0] = (psi[:,1]-psi[:,nx-1])/(2*dx)
	#dely[:,nx-1] = dely[:,0]# HOLTON version
	delx[:,nx-1] = (psi[:,0]-psi[:,nx-2])/(2*dx) #MY version
	
	#If I want symmetry at the periodic boundary...
	#then, there's a difference, if I do calculate it from nx, and copy that to 0, vs calculate from 0, and copy to nx
	#i think that must be asymmetry
	return delx,dely


def PoissonSolverPeriodicBoundaries(psi,nx,ny):
	
	return 0

def HardWallBoundaries(u,v,psi,zeta,dx,dy):

	#f[y,x], boundary conditions for velocity
	#no-flow across boundary
	u[:,0] = 0 #western border
	u[:,-1] = 0 #eastern border
	v[0,:] = 0 #southern border
	v[-1,:] = 0 #northern border
	
	#Full-Slip Forward differences on N,S,E,W BC, like Holton grad.m and some other papers...
	#Try for simplicity...
	u[0,:] = -(psi[1,:] - psi[0,:])/dy #south
	u[-1,:] = -(psi[-1,:] - psi[-2,:])/dy #north
	v[:,0] = (psi[:,1] - psi[:,0])/dx #west
	v[:,-1] = (psi[:,-1] - psi[:,-2])/dx #east
	
	
	#no slip,
	#u[0,:] = 0 #southern border
	#u[-1,:] = 0 #northern border
	#v[:,0] = 0 #western border
	#v[:,-1] = 0 #eastern border
	

	
	
	#Actually, maybe already HERE; i should put in BCs...
	#because Jacobian will use them!...
	
	#Maybe also BCs before first timestep... like maybe plot initial zeta....
	#as 3dplot...
	#check, how FAST does the instability grow
	#Is it already after 2nd step iteration
	#because maybe i should do it BEFORE poisson equation too, since poisson equation will use Zeta...
	#also, try to look at psi solution, from lower z-vertical value......
	#because maybe i can't see if there's a problem in psi, due to scaling of vertical plot axis...
	#so try to change axis....
	#maybe psi actually DOES have also a discontinuity instability at the edge, but it's just very SMALL
	#compared to the large inner mountain peak....
	#so it LOOOKs like boudnary is basically smooth at zero, but inreality, if we zoom in,
	#maybe psi is also kinda wrong near edges....
	

	#zeta Boundary Conditions 
	#Here I haven't set boundary conditions yet..

	#according to https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
	#maybe I can get their result from the jacobian too
	
	#Full-slip zeta BC
	#since u0=u1 and uNY=uNY+1, zeta = -du/dy = 0 on northern and southern BCs
	zeta[0,:] = 0 #south
	zeta[-1,:] = 0 #north
	zeta[:,0] = 0 #west
	zeta[:,-1] = 0 #east
	
	
	#No-slip zeta BC
	#GHOST CELLS OR BOUNDARY???
	#zeta = -du/dy
	#dv/dx = 0, because v = 0, due to no-flow across boundary
	#south, along the entire edge
	#zeta[0,:] = (2/dy**2)*(psi[1,:]-psi[0,:])
	#zeta[0,:] = -(u[1,:]-u[0,:])/dy #south
	#north, along the entire edge
	#zeta[-1,:] = (2/dy**2)*(psi[-1,:]-psi[-2,:])
	#zeta[-1,:] = -(u[-1,:]-u[-2,:])/dy #north
	
	
	#zeta = dv/dx
	#-du/dy = 0, because u = 0, due to no-flow across boundary
	#west, along entire egde
	#zeta[:,0] = (2/dx**2)*(psi[:,1]-psi[:,0]) #west
	#zeta[:,0] = (v[:,1]-v[:,0])/dx #west
	#east along entire edge
	#zeta[:,-1] = (2/dx**2)*(psi[:,-1]-psi[:,-2]) #east
	#zeta[:,-1] = (v[:,-1]-v[:,-2])/dx #east
	
	
	
	#2nd version no-slip BC condition on vorticity
	#https://web.math.princeton.edu/~weinan/papers/cfd5.pdf
	#zeta[0,:] = (2/dy**2)*psi[1,:]
	#zeta[-1,:] = (2/dy**2)*psi[-2,:] #maybe this one should have -minus infront??
	
	#zeta[:,0] = (2/dx**2)*psi[:,1]
	#zeta[:,-1] = (2/dx**2)*psi[:,-2] #maybe this one should have -minus infront??
	
	#Free-slip BCs (also called full-slip sometimes, i think)
	#https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
	#u0 = u1 and uNY+1 = uNY
	#aka, on north and south, velocity is not zero,
	
	
	#No-flow across boundary
	#https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
	#on north and south, v = 0
	#then psi0 = psiNY+1 = constant in x, 
	
	
	return u,v,zeta

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
	for r in range(2000): #pseudo-time
	
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
		
		
		#bottom and top, dpsi/dn = 0
		#psi[ny-2,:] = psi[ny-1,:]
		#psi[1,:] = psi[0,:]
		#or, reverse, maybe?
		#psi[ny-1,:] = psi[ny-2,:]
		#psi[0,:] = psi[1,:]
		
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


def PsiSolvePoissonJacobi(psi,zeta,dx,dy,nx,ny,epstol,Ngc):
	"""
	Solve Poisson Equation
	Holton page 465
	
	psi is streamfunction
	zeta is vorticity
	
	nabla^2psi = zeta
	
	where we are solving for psi the streamfunction, and zeta is the relative vorticity
	
	
	hmm i think zeta is relative vorticity right, but we actually use another Zeta thing too???
	
	hmm... zeta er fatisk absolute vorticity
	og relative vorticity er Xi!
	hmm ok...
	Holton bruger faktisk zeta til absolut vorticity....
	Men, nogen steder siges det,
	Spørgsmålet er om poisson skal løses for
	nabla^2psi = xi
	eller
	nabla^2psi = zeta
	
	fordi, zeta = xi+f
	tror jeg...
	
	
	ahh, absolute vorticity i Holton er
	eta = zeta+f
	så zeta ER relative vorticity
	så (11.14) er aboslute vorticity
	og jeg tror indeed at (13.26) på en måde er set ud fra relative vorticity zeta, fordi at f part er sat i F(x,y,t)
	så det er nok fra absolute vorticity zeta+f, men så har man lavet det til en PDE for zeta
	
	
	
	
	
	We set BCs at y = top and y = bottom, but not x, i think
	Following lec12.pdf here
	https://www.twitch.tv/greekgodx
	
	From lec12.pdf,
	Solving poisson equation means putting some boundary condition on psi, likely
	Because otherwise you could add any constant to psi, and it would still be able to give some curvature
	nabla^2psi = f...
	so it makes sense, to give some unique solution, we impose boundary conditions, 
	but i think, for the sake of calculating velocity, it also wouldn't matter if we added some constant to psi
	
	
	right now i'm doing psi[1:ny-1,1:nx-1]
	#So that means i'm not setting psi[bottom,x] but I AM setting psi[top,x]
	#So i think it should be psi[1:ny-2,1:nx-1]
	#Or maybe, x should also be changed...
	But, remember in Python, 1:ny-1, then 1 is included, but ny-1 is excluded, so in effect, we're doing 1:ny-2
	
	
	
	Solve Poisson Equation,
	nabla^2psi = zeta
	Elliptic PDE
	It is the steady state solution of the pseudo-time problem,
	dpsi/dt = nabla^2psi - zeta
	
	FTCS = Forward Time, Centered Space
	
	#psi1[m,n] = psi[m,n]+dt/dx^2 (psi[m-1,n]+psi[m+1,n]+psi[m,n-1]+psi[m,n+1]-4psi[m,n])-dt*zeta[m,n]
	
	#stability requires, dt <= 1/4 * dx^2
	
	Jacobi Method is choosing dt = 1/4 * dx^2
	
	
	
	
	
	
	https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
	nabla^2Ø - G0 = 0
	1) initial field, eg Ø = 0, which is what I have by psin = np.zeros()
	2) compute error e0 for each grid point... so it be practically an array...
	
	
	
	Should solve Poisson WITH beta plane
	http://www.math.ualberta.ca/ijnam/Volume-10-2013/No-3-13/2013-03-04.pdf
	nabla^2psi = zeta-beta*y
	
	Både Holton og Durran de solver for nabla^2psi = zeta
	de solver for RELATIVE vorticity
	
	Yeah, jeg kan se, forskellige steder bruger forskellig notation.....
	http://empslocal.ex.ac.uk/people/staff/dbs202/cag/courses/MTMW14/notes2006.pdf
	her bruger de Xi til relative vorticity,
	men det er SAMME PDE system overall, det er så nabla^2psi = xi
	så, overall, er det altid nabla^2psi = relative vorticity
	
	
	
	Boundary Conditions
	Lec12, Dirichlet BCs psi = 0
	http://twister.caps.ou.edu/CFD2003/Phillips_NLInstablity.pdf
	at boundaries, psi=0 and xi=0 for all time
	
	http://empslocal.ex.ac.uk/people/staff/dbs202/cag/courses/MTMW14/notes2006.pdf
	Free-slip BCs at north and south, means v = 0, and du/dy = 0, 
	this leads to psi = constant and xi = 0 on north and south.
	
	
	"""
	#Compute streamfunction from vorticity
	#ζ_(i,j)^1=(((ψ_(i+1,j)-2ψ_(i,j)+ψ_(i-1,j) ))/dx^2 +((ψ_(i,j+1)-2ψ_(i,j)+ψ_(i,j-1) ))/dy^2 )

	#solve psi from lec12.pdf
	psin = np.zeros((ny,nx),dtype=np.float64)
	psin[:,:] = psi.copy()
	
	#error betweeen iterated solution and RHS at each gridpoint
	eps = np.zeros((ny,nx),dtype=np.float64)
	
	#Maybe not needed for some implementations
	#Nabla^2psi
	nabla2psi = np.zeros((ny,nx),dtype=np.float64)

	
	##JACOBI METHOD
	#dtau = (1/4)*(0.5*dx**2+0.5*dy**2)
	nit = 0
	while True:
	#instead for r in range(500): #pseudo-time
	
		#update iteration counter
		nit +=1
		
		
		#Make copy of current iterate
		#We should probably do a copy, because, otherwise, we are doing intermediate changes...
		#i.e we are changing values that will be used later, that shouldn't....
		#So, that is, to really follow the algorithm, we need to copy the values...
		#At least, most likely... if it was a double for loop, then we would have to....
		#BUT, MAYBE, we don't have to... maybe numpy is smart enough to figure it out
		psin = psi.copy()
		
		
		
		#Calculate intermediate array
		#1/dx**2(psi[1:ny-1,2:nx]+psi[1:ny-1,0:nx-2]-2*psi[1:ny-1,1:nx-1]) +
		#1/dy**2(psi[2:ny,1:nx-1]+psi[0:ny-2,1:nx-1]-2*psi[1:ny-1,1:nx-1])
		nabla2psi[1+Ngc:ny-1-Ngc,Ngc+1:nx-1-Ngc]= ((1/dx**2)*(psi[1+Ngc:ny-1-Ngc,2+Ngc:nx-Ngc]
								+psi[1+Ngc:ny-1-Ngc,0+Ngc:nx-2-Ngc]-2*psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc])
								+(1/dy**2)*(psi[2+Ngc:ny-Ngc,1+Ngc:nx-1-Ngc]
								+psi[0+Ngc:ny-2-Ngc,1+Ngc:nx-1-Ngc]-2*psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]))
		
		#Calculate error
		eps[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc] = nabla2psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]-zeta[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]
		
		
		#Update iteration of streamfunction
		psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc] = psin[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]+eps[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]/(2/dx**2+2/dy**2)
		
		####
		psi[1:ny-1,0] = 0.25*(psin[1:ny-1,1]+psin[1:ny-1,nx-1]+psin[2:ny,0]+psin[0:ny-2,0])-(dx**2/4)*zeta[1:ny-1,0]
		#periodic BCs
		#psi[1:ny-1,nx-1] = psi[1:ny-1,0]
		psi[1:ny-1,nx-1] = 0.25*(psin[1:ny-1,0]+psin[1:ny-1,nx-2]+psin[2:ny,nx-1]+psin[0:ny-2,nx-1])-(dx**2/4)*zeta[1:ny-1,nx-1]
		#Interior points
		#psi[1:ny-1,1:nx-1] = psin[1:ny-1,1:nx-1]+dtau*(
		#		+(psin[1:ny-1,2:nx]-2*psin[1:ny-1,1:nx-1]+psin[1:ny-1,0:nx-2])/dx**2
		#		+(psin[2:ny,1:nx-1]-2*psin[1:ny-1,1:nx-1]+psin[0:ny-2,1:nx-1])/dy**2
		#		-zeta[1:ny-1,1:nx-1])
		
		
		#https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
		#psi[1:ny-1,1:nx-1] = ((1/dx**2)*(psi[1:ny-1,2:nx]+psi[1:ny-1,0:nx-2]-2*psi[1:ny-1,1:nx-1])
		#						+(1/dy**2)*(psi[2:ny,1:nx-1]+psi[0:ny-2,1:nx-1]-2*psi[1:ny-1,1:nx-1])
		#						-zeta[1:ny-1,1:nx-1])/(2/dx**2+2/dy**2)
		
		
		
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
		
		
		#Psi no-flow across boundary ACTUAL BOUNDARY
		#psi[:,-1+Ngc] = 0 #right boundary
		#psi[:,0+Ngc] = 0	#left boundary
		#psi[0+Ngc,:] = 0  #south boundary
		#psi[-1-Ngc,:] = 0 #north boundary
		
		#
		
		#No-slip psi streamfunction BC
		#Maybe I need additional BCs due to the no-slip condition
		#no-slip, hence, on sourthern boundary,
		#hence, u0 = 0 on sourthern boundary, and ALSO v0 = 0 on sourthern boundary
		#now, u0 = 0, but, u1 != 0 likely...
		#so, on southern boundary, there would be, du/dy...
		#hence, there WOULD be a psi and a zeta, right
		#ie, for no-slip, then zeta is NOT zero right, at boundary
		#yep... for no-slip, zeta is NOT zero at boundary, and i think, psi is ALSO a bit more complicated
		#maybe,
		#psi[:,1] = psi[:,0]+u[:,0]*dy
		#something like that..
		#actually hmm it might be good enough.... i actually think it's okay,
		#psi[:,1] will come from solving poisson equation, 
		
		

		
		#lec12, on sourthern boundary, n = (0,1), hence nablapsi.n = 0 = dpsi/dy = 0
		#yeah... that will give, u = -dpsi/dy = 0 on boundary....
		#but that almost means, that, it should be, psi0 = psi1 = 0...
		#like, psi0,i = psi1,i = 0...
		#So i should actually enfore the psi to be zero even further out, i THINK...
		#because ONE thing is ghost cells, those are not the ACTUAL boundary, i think...
		#they are a pseudo-boundary... the ACTUAL boundary is even further inside,
		#so i think i need to enfore psi = 0 even further inside that what i currently have...
		#also, hence, the Poisson solver, should actually solve LESS area than what it doesn
		#I need to reduce the area by +-1 in each direction i think...
		#like, there's BOUNDARY conditions, AND, ghost cell conditions...
		#so just think like, psi[:,0] is not an ACTUAL BOUNDARY condition, it's actually a GHOST CELL CONDITION
		#so i need to ALSO set psi[:,1] as the ACTUAL boundary condition, AND change the algorithms...
		


		#so... if i want u = 0 on southern boundary, I would need,
		#u = 0 = -(psi(1)-psi(-1))/dy
		#so maybe I DO need ghost cells or something...
		#and then I would have psi(-1) = psi(0) = 0
		#psi(0) = 0 from dirichlet psi = 0 condition on BCs
		#psi(-1) = 0 from no-slipp conditions...
		#
		#Like, det er lidt som, it's a bit like,
		#hvis man har en vorticity-streamfunction formulation
		#så SKAL man jo nærmest have en staggered grid
		#almost JUST HAPPENS to be that way...
		#because, u = dpsi/dx... so... it's like... i mean maybe you could FORCE it to be unstaggered,
		#but, by it's very nature, it seems staggered
		#so, how do we deal with that then
		
		
		#psi[:,-1-1] = 0 #right boundary
		#psi[:,0+1] = 0	#left boundary
		#psi[0+1,:] = 0  #south boundary
		#psi[-1-1,:] = 0 #north boundary
		
		
		#Now me tihnking, ghost cells or just complicated boundary conditions
		#Well actually, ghost cells are also used in eg Semi-Lagrangian CFD...
		#so, it's actually maybe rather common in some ways...
		#and it also speaks to the idea of simulating a larger domain, it's sort of same philosophy
		#so imo, might aswell just get used to it, and start adding ghost cells a bit here...
		#
		
		

		
		#print(eps)
		if (np.abs(eps) < epstol).all():
			print("Break Poisson algorithm at", nit)
			break
		
		
	
	#I think here, we should return psi, instead of doing it implicitly....
	#so, change the way the method works...
	return psi


def divfluxHolton(psi,u,v,dx,dy,ny,nx):

# % Take forward differences on top and bottom edges
# dfly(1,:) = 0*(P(2,:).*v(2,:) - P(1,:).*v(1,:))/dy;
# dfly(n,:) = 0*(P(n,:).*v(n,:) - P(n-1,:).*v(n-1,:))/dy;

# % Take centered differences on interior points

# dfly(2:n-1,:) = (P(3:n,:).*v(3:n,:)-P(1:n-2,:).*v(1:n-2,:))/(2*dy);
# % compute x component of flux 

# % Take cyclic differences on left and right boundaries

# dflx(:,1)=(P(:,2).*u(:,2)-P(:,p-1).*u(:,p-1))/(2*dx);
# dflx(:,p)= dflx(:,1);

# % take centered differences on interior points

# dflx(:,2:p-1)= (P(:,3:p).*u(:,3:p)-P(:,1:p-2).*u(:,1:p-2))/(2*dx);


	#Forward difference on top and bottom
	dflx = np.zeros((ny,nx),dtype=np.float64)
	dfly = np.zeros((ny,nx),dtype=np.float64)
	
	
	#TOp and bottom edges, forward differences
	dfly[0,:] = 0*(psi[1,:]*v[1,:]-psi[0,:]*v[0,:])/dy
	dfly[ny-1,:] = 0*(psi[ny-1,:]*v[ny-1,:]-psi[ny-2,:]*v[ny-2,:])/dy
	
	#Interior points, centered differences
	dfly[1:ny-1,:] = (psi[2:ny,:]*v[2:ny,:]-psi[0:ny-2,:]*v[0:ny-2,:])/(2*dy)
	
	
	#Cyclic difference on boundaries
	dflx[:,0] = (psi[:,1]*u[:,1]-psi[:,ny-1]*u[:,ny-1])/(2*dx)
	dflx[:,ny-1] = dflx[:,0] #HOLTON version
	#dflx[:,ny-1] = (psi[:,0]*u[:,0]-psi[:,ny-2]*u[:,ny-2])/(2*dx) #MY version
	#Interior points, centered differences
	
	dflx[:,1:ny-1] = (psi[:,2:ny]*u[:,2:ny]-psi[:,0:ny-2]*u[:,0:ny-2])/(2*dx)
	
	return dflx,dfly


def JacobianHolton(u,v,zeta,psi,beta,ny,nx):
	"""
	To update relative vorticity zeta

	J(p,q) = dp/dx * dq/dy - dp/dy * dq/dx

	One approximation is,
	Jhat(p,q) = (d2xp)(d2yq)-(d2yp)(d2xq)


	For BVE, J(psi, nabla^2psi)
	So, that corresponds to using u,v and zeta = nabla^2psi


	Also, it should probably actually be,
	Lec12,
	J(psi, f+zeta)
	#So, maybe J(psi, f + nabla^2psi)


	"""

	#Follows Flux version, Holton pg 465
	Fmn = np.zeros((ny,nx),dtype=np.float64)
	Fmn[1:ny-1,1:nx-1] = (
		(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)
		+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)
		+beta*v[1:ny-1,1:nx-1])
		
		
	
	
	#AT THE BORDERS
	##Fmn[0,:] = 
	Fmn[1:ny-1,0] = (
		(u[1:ny-1,1]*zeta[1:ny-1,1]-u[1:ny-1,nx-1]*zeta[1:ny-1,nx-1])/(2*dx)
		+(v[2:ny,0]*zeta[2:ny,0]-v[0:ny-2,0]*zeta[0:ny-2,0])/(2*dy)
		+beta*v[1:ny-1,0])
		
	Fmn[1:ny-1,nx-1] = (
		(u[1:ny-1,0]*zeta[1:ny-1,0]-u[1:ny-1,nx-2]*zeta[1:ny-1,nx-2])/(2*dx)
		+(v[2:ny,nx-1]*zeta[2:ny,nx-1]-v[0:ny-2,nx-1]*zeta[0:ny-2,nx-1])/(2*dy)
		+beta*v[1:ny-1,nx-1])
	
	#Mean vorticity is conserved, except for errors from time diffencing (Holton pg 466)

	return Fmn
	
	
	
def ArakawaJacobian(psi,zeta,dx,dy,nx,ny,beta,v):
	"""
	To update relative vorticity zeta
	
	J(p,q) = dp/dx * dq/dy - dp/dy * dq/dx
	
	One approximation is,
	Jhat(p,q) = (d2xp)(d2yq)-(d2yp)(d2xq)
	
	
	For BVE, J(psi, nabla^2psi)
	So, that corresponds to using u,v and zeta = nabla^2psi
	
	
	Also, it should probably actually be,
	Lec12,
	J(psi, f+zeta)
	#So, maybe J(psi, f + nabla^2psi)
	
	
	
	Final THESIS
	the Arakawa Jacobian must vanish at the walls? y =+-Y
	
	
	"""
	
	
	#https://dspace.library.uvic.ca/bitstream/handle/1828/1218/Final THESIS.pdf?sequence=1
	#J1,J2,J3 aka Jpp, Jpx, Jxp
	J1 = np.zeros((ny,nx))
	J2 = np.zeros((ny,nx))
	J3 = np.zeros((ny,nx))
	JA = np.zeros((ny,nx))
	for j in range(1,ny-1):
		for i in range(1,nx-1):
			J1[j,i] = (1/(4*dx**2))*(
								(psi[j,i+1]-psi[j,i-1])*(zeta[j+1,i]-zeta[j-1,i])
								-(zeta[j,i+1]-zeta[j,i-1])*(psi[j+1,i]-psi[j-1,i]))
								
								
			J2[j,i] = (1/(4*dx**2))*(psi[j,i+1]*(zeta[j+1,i+1]-zeta[j-1,i+1])
								-psi[j,i-1]*(zeta[j+1,i-1]-zeta[j-1,i-1])
								-psi[j+1,i]*(zeta[j+1,i+1]-zeta[j+1,i-1])
								+psi[j-1,i]*(zeta[j-1,i+1]-zeta[j-1,i-1])
								)
								
			J3[j,i] = (1/(4*dx**2))*(zeta[j+1,i]*(psi[j+1,i+1]-psi[j+1,i-1])
								-zeta[j-1,i]*(psi[j-1,i+1]-psi[j-1,i-1])
								-zeta[j,i+1]*(psi[j+1,i+1]-psi[j-1,i+1])
								+zeta[j,i-1]*(psi[j+1,i-1]-psi[j-1,i-1])
								)
	
	
	JA = (1/3)*(J1+J2+J3) + beta*v[:,:]
	Fmn = JA
	
	return Fmn




def nabladotv(u,v,nx,ny,dx,dy):
	"""
	Divergence should be 0 for BVE,
	nabla.v =0
	
	Test it out
	Should maybe be centered differences, since psi streamfunction is centered differences
	And maybe on edges, do forward differences
	"""
	
	
	divv = (u[:,1:nx-1]-u[:,0:nx-2])/dx + (v[1:ny-1,:]-v[0:ny-2,:])/dy

	return divv
	



def SimpleFilter(zetanew,zeta,zetaold):
	"""
	Filters out leapfrog computational mode
	Simple filter, from lec12.pdf
	
	It actually filters the INTERMEDIATE step, zeta, instead of zetanew!
	And then the filtered zeta, will be used for the next timestep
	"""
	zeta = 0.5*(zetanew+zetaold)
	
	return zeta
	
	


def RobertAsselinFilter(zetanew,zeta,zetaold):
	"""
	Robert-Asselin filter
	Often used with Leap-frog scheme
	http://weather.ou.edu/~ekalnay/NWPChapter3/Ch3_2_4.html
	
	
	We filter on zeta, NOT zetanew... it's the current intermediate step, zeta, that will be filtered

	
	http://www.scottsarra.org/math/papers/Shawn%20Cheeks%20Math%20Capstone%202015.pdf
	First do regular leapfrog integration
	Then apply the filter
	F(t)tbar = F(t) + 0.5*nu*(F(t-dt)bar - 2*F(t) + F(t+dt))
	
	
	#https://www.mathematics.pitt.edu/sites/default/files/SurveyTimeFilters1_technical-report.pdf
	#RA-leapfrog
	#u is filtered, v is unfiltered
	#vnew = uold + 2*dt*F(v)
	#u = v + 0.5*nu*(vnew-2*v+uold)
	#So, you can actually do it, without much intrusion at all, and I don't think another array is needed

	
	
	"""
	

	nu_filter = 0.1
	
	
	#apply filter to zeta
	#un = vn + (nu/2)*(vn+1-2vn+un-1)
	#zeta[1:ny-1,1:nx-1] = zetav[1:ny-1,1:nx-1] + (nu_filter/2)*(zetavnew[1:ny-1,1:nx-1]-2*zetav[1:ny-1,1:nx-1]+zetaold[1:ny-1,1:nx-1])
	zeta[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1] + (nu_filter/2.0)*(zetanew[1:ny-1,1:nx-1]-2.0*zeta[1:ny-1,1:nx-1]+zetaold[1:ny-1,1:nx-1])
		
	
	#now, we have
	#filtered zeta
	#unfiltered zetaold
	#unfiltered zetanew

	#Next time, 
	#zetaold = filtered zeta
	#zeta = unfiltered zetanew
	#
	
	
	#https://maths.ucd.ie/~plynch/LECTURE-NOTES/NWP-2004/NWP-CH03-2-3-P4.pdf
	#Unf = Un + gamma*(Un+1-2Un+Un-1f)
	#So, for this only, I probably really only need Un-1f 
	#Actually, on RHS, only 3 unique matrices are needed...
	#Maybe then the filter doesn't actually require extra memory data storage
	#Or maybe it DOES require some intermediate matrix data storage... but, maybe not 2, as I do now..
	
	
	return zeta



def NoSlipBCvorticity(zeta,psi,dx,dy):

	#No-slip zeta BC
	#GHOST CELLS OR BOUNDARY
	#zeta = -du/dy
	#dv/dx = 0, because v = 0, due to no-flow across boundary
	#south, along the entire edge
	zeta[0,:] = (2/dy**2)*(psi[1,:]-psi[0,:])
	#zeta[0,:] = -(u[1,:]-u[0,:])/dy #south
	#north, along the entire edge
	zeta[-1,:] = (2/dy**2)*(psi[-1,:]-psi[-2,:])
	#zeta[-1,:] = -(u[-1,:]-u[-2,:])/dy #north
	
	
	#zeta = dv/dx
	#-du/dy = 0, because u = 0, due to no-flow across boundary
	#west, along entire egde
	zeta[:,0] = (2/dx**2)*(psi[:,1]-psi[:,0]) #west
	#zeta[:,0] = (v[:,1]-v[:,0])/dx #west
	#east along entire edge
	zeta[:,-1] = (2/dx**2)*(psi[:,-1]-psi[:,-2]) #east
	#zeta[:,-1] = (v[:,-1]-v[:,-2])/dx #east
	
	return zeta
	
	
	


def AddOrography(zeta,u,v,h,f,H,dx,dy):
	"""
	Centered differences for dh/dx and dh/dy
	"""
	oroterm = -(f/H)*(u[1:ny-1,1:nx-1]*(h[:,2:nx]-h[:,0:nx-2])/(2*dx)+v[1:ny-1,1:nx-1]*(h[2:ny,:]-h[0:ny-2,:])/(2*dy))
	
	
	return oroterm




def CalcVelocity(psi,u,v,dx,dy):
	"""
	Should take u,v as input, or initiate u,v here...
	ug,vg geostrophic
	
	ug = -dpsi/dy
	vg = dpsi/dx
	
	these velocities depend on streamfunction psi being correct, especially near BCs...
	BCs are set outside a bit...
	#so v[0,:] eg is outside this function...
	
	
	"""
	
	#psi[y,x] array order
	#for i in range(1,nx-1):
	#	for j in range(1,ny-1):
	#		u[j,i] = -(psi[j+1,i]-psi[j-1,i])/(2*dy)
	#		v[j,i] = (psi[j,i+1]-psi[j,i-1])/(2*dx)
	

	u[1:ny-1,1:nx-1] = -(psi[2:ny,1:nx-1]-psi[:ny-2,1:nx-1])/(2*dy)
	v[1:ny-1,1:nx-1] = (psi[1:ny-1,2:nx]-psi[1:ny-1,0:nx-2])/(2*dx)
		
	return u,v



def CalcVelocityPeriodicBoundaries(psi,u,v,dx,dy):
	"""
	Should take u,v as input, or initiate u,v here...
	ug,vg geostrophic
	
	ug = -dpsi/dy
	vg = dpsi/dx
	
	these velocities depend on streamfunction psi being correct, especially near BCs...
	BCs are set outside a bit...
	#so v[0,:] eg is outside this function...
	"""
	
	#psi[y,x] array order
	#for i in range(1,nx-1):
	#	for j in range(1,ny-1):
	#		u[j,i] = -(psi[j+1,i]-psi[j-1,i])/(2*dy)
	#		v[j,i] = (psi[j,i+1]-psi[j,i-1])/(2*dx)
	
	u[1:ny-1,0] = -(psi[2:ny,1:nx-1]-psi[:ny-2,1:nx-1])/(2*dy)
	#u[1:ny-1,1:nx-1] = -(psi[2:ny,1:nx-1]-psi[:ny-2,1:nx-1])/(2*dy)
	#v[1:ny-1,1:nx-1] = (psi[1:ny-1,2:nx]-psi[1:ny-1,0:nx-2])/(2*dx)
		
	return u,v

if __name__ == "__main__":

	import numpy as np
	#from pylab import *
	import matplotlib.pyplot as plt
	import matplotlib.animation as animation
	from mpl_toolkits.mplot3d import axes3d
	import matplotlib
	import scipy
	
	#np.set_printoptions(precision=2)

	Backend = matplotlib.get_backend()
	print(Backend)#We use TkAgg
	#matplotlib.use("wx") #No wx
	matplotlib.use("TkAgg")
	from pylab import *
	thismanager = get_current_fig_manager()
	#thismanager.window.SetPosition((500, 0))
	thismanager.window.wm_geometry("-1500+0")

	
	#Er det mon nogen leapfrog scheme conditions jeg ikke overholder
	
	Lx = 6000#6000 #km
	Ly = 3000#3000 #km
	nx = int(3*65)
	ny = int(3*65)
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
	U0 = 0#20#20				#zonal wind
	
	R_earth = 6381*10**3 #meters
	Omega_earth = 2*np.pi/(24*60*60)
	phi0 = np.pi/4 #45degrees
	#beta = df/dy = 2*Omega_Earth*cos(phi0)/R_earth
	beta = 2*Omega_earth*np.cos(phi0)/R_earth
	
	print("coriolis beta ", beta)
	#mid-latitude beta plane approximation
	#holton page 160
	#f = f0+beta*y
	#with y = 0 at phi0
	#so, absolute vorticity = relative vorticity + f
	
	
	#pg 86
	#Coriolis term
	f0 = 2*Omega_earth*np.sin(phi0)
	f = f0 + beta*y
	
	
	#beta = 1e-11#1.62*10**(-11)	#he set 0 infront?
	Av4 = 10**(-6)
	A = 10**(-4)
	#initial vorticity and streamfunction
	
	
	#Poisson solver error
	epstol = 1e-8



	
	#time integration parameters					#hours
	time_end = 48*60*60			#seconds
	dt = 900 #seconds
	#dt = 0.8*dx/10
	print("TIMESTEP",dt)
	t = 0
	nt = 0
	
	#not NEEDED, but,could psi array should probably be nx+1,ny+1, or, velocity should be nx-1,ny-1
	#Due to centered differences
	psi = np.zeros((ny,nx),dtype=np.float64)
	dypsi = np.zeros((ny,nx),dtype=np.float64)
	dxpsi = np.zeros((ny,nx),dtype=np.float64)
	u = np.zeros((ny,nx),dtype=np.float64)
	v = np.zeros((ny,nx),dtype=np.float64)
	#dfly = np.zeros((ny,nx),dtype=np.float64)
	#dflx = np.zeros((ny,nx),dtype=np.float64)





	
	
	
	#===============================================
	#Initial conditions from psi Rossby-haurwitz
	
	#Rossby wave-packets
	#psi[:,:] = np.cos(k1*x-w*t)*sin(k2*y)
	#psi[:,:] = np.cos(k*x)*np.sin(m*y)
	#http://www.math.ualberta.ca/ijnam/Volume-10-2013/No-3-13/2013-03-04.pdf
	#Rossby wave packet satisfies Dirichlet Boundary condition, but those are specific...
	#So I think I need simpler IC...
	
	#https://dspace.library.uvic.ca/bitstream/handle/1828/1218/Final%20THESIS.pdf?sequence=1
	#Barotropic Rossby Waves psi(x,y) is prop to RE(exp(i(kx+ly-wt)))
	#1)Barotropic Rossby Waves psi(x,y) = A*RE(exp(i(kx+ly-wt)))
	#2)Planetary Wave Packets psi(x,y) = a*cos(k1x-wt)*sin(k2y)
	#2) solve equatorial BVE exactly
	
	#Okay, I want, psi = 0 on all edges, so cos() and sin() must be 0 on edges...
	#Now, sin(0) = 0, and cos(0)=1, so they are opposite in their usage, and thus k1,k2 should be different prolly
	#Should probably plot initial conditions... of psi...
	
	#2) has v = 0, on walls Y = +-
	#IF, k2 = m1*pi/Y, where m1 is some multiple
	
	#psi = g/f * ø
	#Or, zeta = 10^-4s
	#
	k1 = 2*pi/(Lx*1000)
	k2 = 1*pi/(Ly*1000)
	Amag = 10**(7)
	psi[:,:] = Amag*np.cos(k1*x)*np.sin(k2*y) #Det er solution til equatorial BVE? EBVE? Den har ingen beta? Eller er beta sat til = 1? natural units?
	#men, de bruger xi i 2.74, mens de bruger zeta i 4.1...
	#så måske, zeta = xi+f, men ikke sure
	#YUP! OKAY!
	#De bruger zeta = xi + f i 2.87
	#Så, jeg tror, 4.1, det er MED coriolis f
	#Og yep, de solver furthermore, poisson MINUS y...
	#3.1-3.3, bruger zeta = xi+f
	
	#og section 3.3 Free Equatiorial bruger de zeta, ikke xi...okay, så det er MED coriolis det her...
	#men IKKE mean flow u0 tror jeg dog...
	#De viser at Free EBVE kan skrives på 3 måder...
	
	w = U0*k1-beta*k1/(k1**2+k2**2)

	#for jacobi method
	psizeros = np.zeros((ny,nx),dtype=np.float64)
	print("ZETA MAG")
	print((k1**2+k2**2)*Amag)
	psi_exact_initial = np.zeros((ny,nx),dtype=np.float64)
	psi_exact_initial[:,:] = np.cos(k1*x-w*0)*np.sin(k2*y)*Amag
	#https://www.wolframalpha.com/input/?i=d%5E2%2Fdx%5E2cos%28ax-wt%29sin%28by%29%2Bd%5E2%2Fdy%5E2cos%28ax-wt%29sin%28by%29
	zeta_exact_initial = -(k1**2+k2**2)*(np.cos(k1*x-w*0)*np.sin(k2*y))*Amag
	zeta_exact_initial[0,:] = 0
	zeta_exact_initial[ny-1,:] = 0
	
	#https://www.wolframalpha.com/input/?i=d%2Fdx+cos%28ax-wt%29sin%28by%29
	uexact = -k1*np.sin(k2*y)*np.sin(k1*x-w*0)
	
	#Hmm, not quite 0 at the borders, so maybe i enforce it? Maybe numerical precision??
	#psi[0,:] = 0
	#psi[-1,:] = 0
	#psi[:,0] = 0
	#psi[:,-1] = 0
	
	v,u = gradHolton(u,v,psi,dx,dy,nx,ny)
	
	#==============================================
	#Add orography
	H = 8000
	h = np.zeros((ny,nx),dtype=np.float64)
	h0 = 6000
	x0m = 1000*Lx/2
	y0m = 0#1000*Ly/2
	#f = 1e-4
	r = 500*1000 #500 km radius
	k1_oro = 0.5*pi/(Lx*1000)
	k2_oro = 0.5*pi/(Ly*1000)
	h[:,:] = h0*np.exp(-((x-x0m)**2+(y-y0m)**2)/(2*r**2))+h0*np.exp(-((x-x0m/2)**2+(y-y0m)**2)/(2*r**2))
	
	#h[:,:] = h0*np.sin(k1_oro*x)*np.sin(k2_oro*y)
	#h[:,:] = h0*np.cos()**2
	h[0,:] = 0
	h[ny-1,:] = 0
	h[:,0] = 0
	h[:,nx-1] = 0
	#Best if orography == 0 on boundaries, according to notes2006.pdf
	
	
	#=============================================
	#Full bounded box
	# #start with psi0
	# #then
	# #u0 = -dpsi0/dy
	# #v0 = dpsi0/dx
	# u,v = CalcVelocity(psi,u,v,dx,dy)
	
	
	# #Rescale initial speed
	
	
	# #f[y,x], boundary conditions for velocity
	# #no-flow across boundary
	# u[:,0] = 0 #western border
	# u[:,-1] = 0 #eastern border
	# v[0,:] = 0 #southern border
	# v[-1,:] = 0 #northern border
	# #then
	# #zeta0 = du0/dy - dv0/dx
	# #actually pretty straight forward, also equilevant to nabla^2psi0 = zeta0 the forward poisson problem...
	# #
	# #and then you can also rescale the velocity, to give e.g 5m/s
	
	
	# #Full-Slip Forward differences on N,S,E,W BC, like Holton grad.m and some other papers...
	# #Try for simplicity...
	# u[0,:] = -(psi[1,:] - psi[0,:])/dy #south
	# u[-1,:] = -(psi[-1,:] - psi[-2,:])/dy #north
	# v[:,0] = (psi[:,1] - psi[:,0])/dx #east
	# v[:,-1] = (psi[:,-1] - psi[:,-2])/dx #west
	#==========================================================
	
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
	
	
	
	#We have 3 sets of zetas here, one is initial,zeta0
	#zeto0 is set to a gaussian it
	#zeta0 = np.array(A*np.exp(-2*(k**2*x**2+m**2*y**2)),dtype=np.float64)
	zeta0 = np.zeros((ny,nx),dtype=np.float64)
	zeta = zeta0.copy()
	zetanew = zeta0.copy()
	
	#We can calculate zeta if given initial psi, just forward method
	#Right now, doesn't calculate edges,
	#so just by default, zeta = 0 on edges
	for i in range(1,nx-1):
		for j in range(1,ny-1):
			#u[i,j] = -(psi[i+1,j]-psi[i-1,j])/(2*dy)
			#v[i,j] = (psi[i,j+1]-psi[i,j-1])/(2*dx)
			zeta0[j,i] = (v[j,i+1]-v[j,i-1])/(2*dx)-(u[j+1,i]-u[j-1,i])/(2*dy)
	
	zeta0[1:ny-1,0] = 	(v[1:ny-1,1]-v[1:ny-1,nx-1])/(2*dx)-(u[2:ny,0]-u[0:ny-2,0])/(2*dy)
	#zeta0[1:ny-1,nx-1] = zeta0[1:ny-1,0] #HOLTON style
	zeta0[1:ny-1,nx-1] = (v[1:ny-1,0]-v[1:ny-1,nx-2])/(2*dx)-(u[2:ny,nx-1]-u[0:ny-2,nx-1])/(2*dy) #MY STYLE
	#zeta = 0 on BCs, i think, for FULL SLIP...
	#But, maybe the rossby wave packet defined, is NOT good with zeta = 0 on BCs,
	#maybe those initial conditions do not satisfy the BCs...
	
	
	#We have 3 sets of zetas here, one is initial, zeta0
	#zeto0 is set to a gaussian it
	#zeta0 = np.array(A*np.exp(-2*(k**2*x**2+m**2*y**2)),dtype=np.float64)
	#zeta = zeta0.copy()
	#zetanew = zeta0.copy()

	zeta = zeta_exact_initial.copy()
	zetanew = zeta_exact_initial.copy()

	
	
	print("INITIAL BOUNDARY CONDITIONS")
	print("u[0,:5]", u[0,:5])
	print("u[-1,:5]", u[-1,:5])
	print("u[:5,0]", u[:5,0])
	print("u[:5,-1]", u[:5,-1])
	print("v[0,:5]", v[0,:5])
	print("v[-1,:5]", v[-1,:5])
	print("v[:5,0]", v[:5,0])
	print("v[:5,-1]", v[:5,-1])
	print("psi[0,:5]", psi[0,:5])
	print("psi[-1,:5]", psi[-1,:5])
	print("psi[:5,0]", psi[:5,0])
	print("psi[:5,-1]", psi[:5,-1])
	

	
	#Filter arrays... probably not needed
	zetav = np.zeros((ny,nx),dtype=np.float64)
	zetavnew = np.zeros((ny,nx),dtype=np.float64)
	
	
	#Ghost cells
	Ngc = 0
	
	#Boundary conditions from beginning
	#zeta = NoSlipBCvorticity(zeta,psi,dx,dy)
	
	
	
	
	#===================================================================================
	#Some arrays for enstrophy, energy, etc
	
	EnstrophyList = []
	MeanKEList = []
	MeanZetaList = []

	
	
	
	#====================================================================================
	#Main loop. Integrate forward in time
	
	
	#------------------------------------------------
	#First step
	#Euler forward is easy to implement...
	#If you can accept the large first error
	#Alternatively, you can use, Euler-Backwards (Matsuno) step on the first
	#https://maths.ucd.ie/~plynch/LECTURE-NOTES/NWP-2004/NWP-CH03-2-3-P4.pdf
	#You're basically doing a leapfrog scheme with half-step too i think...
	
	#Euler Forward
	
	t+= dt
	
	psi = psi+U0*(Ly/2*1000 -y)
	v,u = gradHolton(u,v,psi,dx,dy,nx,ny)
			#Try, full-slip Bcs,
	u[0,:] = u[1,:]
	u[ny-1,:] = u[ny-1,:]
	#============================================
	#FULL BOUNDED BOX
	# #Calculate streamfunction psi
	# psi = PsiSolvePoissonJacobi(psi,zeta,dx,dy,nx,ny,epstol,Ngc)
	
	# #calculating velocities
	# #only in interior
	# u,v = CalcVelocity(psi,u,v,dx,dy)
	
	# #f[y,x], boundary conditions for velocity
	# #no-flow across boundary
	# u[:,0] = 0 #western border
	# u[:,-1] = 0 #eastern border
	# v[0,:] = 0 #southern border
	# v[-1,:] = 0 #northern border
	
	# #Full-Slip Forward differences on N,S,E,W BC, like Holton grad.m and some other papers...
	# #Try for simplicity...
	# u[0,:] = -(psi[1,:] - psi[0,:])/dy #south
	# u[-1,:] = -(psi[-1,:] - psi[-2,:])/dy #north
	# v[:,0] = (psi[:,1] - psi[:,0])/dx #west
	# v[:,-1] = (psi[:,-1] - psi[:,-2])/dx #east
	

	# #Full-slip zeta BC
	# #since u0=u1 and uNY=uNY+1, zeta = -du/dy = 0 on northern and southern BCs
	# zeta[0,:] = 0 #south
	# zeta[-1,:] = 0 #north
	# zeta[:,0] = 0 #west
	# zeta[:,-1] = 0 #east
	#==================================================

	#Updating relative vorticity zeta
	#Fmn = ArakawaJacobian(psi,zeta,dx,dy,nx,ny,beta,v)
	#integralFmn = np.sum(Fmn)
	#print("Fmn integral ",integralFmn)
	
	
	#zeta[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1] - dt*Fmn #for normal jacobian
	#zetanew[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1] - dt*Fmn[1:ny-1,1:nx-1] #for Arakawa jacobian
	
	
	#Holton way
	#zetan(2:Nyl,1:Nxl) = zeta(2:Nyl,1:Nxl)-beta*dt*v(2:Nyl,1:Nxl) ...
    #        - dt*(dflx(2:Nyl,1:Nxl)+dfly(2:Nyl,1:Nxl))...
	dflx,dfly = divfluxHolton(psi,u,v,dx,dy,ny,nx)
	#zetanew[1:ny-1,:] = zeta[1:ny-1,:]-dt*(beta*v[1:ny-1,:]+dflx[1:ny-1,:]+dfly[1:ny-1,:])
	#zetanew[1:ny-,0:nx-2] = zeta[1:ny-2,0:nx-2]-dt*(beta*v[1:ny-2,0:nx-2]+dflx[1:ny-2,0:nx-2]+dfly[1:ny-2,0:nx-2])
	
	#zetanew[1:ny-1,0:nx-1] = zeta[1:ny-1,0:nx-1]-dt*(beta*v[1:ny-1,0:nx-1]+dflx[1:ny-1,0:nx-1]+dfly[1:ny-1,0:nx-1])
	Fmn = JacobianHolton(u,v,zeta,psi,beta,ny,nx)
	
	#RA-leapfrog
	#new_filter = old + 2*dt*F(current_filter)
	#vn+1 = un-1 + 2*dt*F(vn)
	#zetavnew[1:ny-1,1:nx-1] = zetaold[1:ny-1,1:nx-1] + 2*dt*(-Fmnv[1:ny-1,1:nx-1])
	
	#Indices should probably change, if doing periodic...
	#Only if bounded box, then zeta = 0 on boundaries pretty much...
	#but with periodic x... or, might actually be fine...
	#zetanew[1:ny-1,0:nx-1] = zetaold[1:ny-1,0:nx-1] + 2*dt*(-Fmn[1:ny-1,0:nx-1])#+
	#zetanew[1:ny-1,0:nx-1] = zeta[1:ny-1,0:nx-1] + dt*(-Fmn[1:ny-1,0:nx-1])
	zetanew[1:ny-1,0:nx] = zeta[1:ny-1,0:nx] + dt*(-Fmn[1:ny-1,0:nx])

	#Update bottom and top
	#zetanew[ny-1,:] = zetanew[ny-2,:]
	#zetanew[0,:] = zetanew[1,:]
	
	
	#Full slip bottom and top
	zetanew[0,:] = -(u[1,:]-u[0,:])/dy
	zetanew[ny-1,:] = -(u[ny-1,:]-u[ny-2,:])/dy
	
	#calc 0th index
	#zetanew[1:ny-1,0] = ...
	#Periodic BCs
	#zetanew[1:ny-1,nx-1] = zetanew[1:ny-1,0] 
	#zetanew[:,nx-1] = zetanew[:,0]
	#Make more awesome BCS
	
	
	print("zeta S edge", zeta[0,0:2],zeta[0,nx-2:nx])
	print("zeta N edge", zeta[-1,0:2],zeta[-1,nx-2:nx])
	print("zeta E edge", zeta[0:4,-1])
	print("zeta W edge", zeta[0:4,0])
	
	print("psi S+1 edge", psi[1,0:4])
	print("psi N-1 edge", psi[-2,0:4])
	
	
	vspeed = np.sqrt(u*u+v*v)
	print("max speed ", np.max(vspeed))
	
	
	#start main loop
	while t < time_end:
		
		#Do a copy, but just assigning may also be work. Depends on algorithms used.
		#old zeta, is equal to the intermediate zetacurrent we had in previous step
		zetaold = zeta.copy()
		#Current zeta, is equal to the zeta we just had...
		zeta = zetanew.copy()


		print(t,nt)
		nt += 1
		t+= dt
		
		#Calculate streamfunction psi
		#psi = PsiSolvePoissonJacobi(psi,zeta,dx,dy,nx,ny,epstol,Ngc,)
		
		#zeta har -(dx**2/4)*zeta... hvis dx**2 blvier stor, måske det skaber instability
		#print("dx**2*zeta")
		#print(dx**2*zeta)
		#Okay, what are the usual values of zeta?
		#It's usually, like, 10^-5 or something right? So yeah, that will explode...
		#so... that could then also mean, that dx should be smaller...
		#stuff like that...
		#so maybe smaller dx,
		#and find out what zeta scale usually is...
		#and maybe try fast fourier transform method.. for chanel problem...
		#fast fourier so weird man...
		#what about...
		#Making a "safer" convergence for Jacobi??? like, smaller pseudo timestep?
		
		#Hmmm... can I take an example, and explode it
		#Like, take Lorena Barbas example with Poisson, and explode it, either with large b values, or large dx and dy values...
		
		#psi = PsiSolvePoissonJacobiChannel(psizeros,zeta,dx,dy,nx,ny,epstol,Ngc,)
		psi = poisson_periodic(psi, zeta, dx, dy)
		#psi = FFTPoisson(nx,ny,Lx,dy,zeta,psi)
		psi = psi+U0*(Ly/2*1000 -y)
		#psi = PsiSolvePoissonJacobi(psi,zeta,dx,dy,nx,ny,epstol,Ngc)
		#print("psi")
		#print(psi)
		#psi = psi + U0*(Ly/2*1000 -y)
		#No-slip?
		#psi[:,-1-1] = 0 #right boundary
		#psi[:,0+1] = 0	#left boundary
		#psi[0+1,:] = 0  #south boundary
		#psi[-1-1,:] = 0 #north boundary
		


		#calculating velocities
		#only in interior
		#u,v = CalcVelocity(psi,u,v,dx,dy)
		
		
		#Bounday conditions
		#u,v,zeta = HardWallBoundaries(u,v,psi,zeta,dx,dy)

		
		v,u = gradHolton(u,v,psi,dx,dy,nx,ny)
		
		#Full-slip Bcs,
		u[0,:] = u[1,:]
		u[ny-1,:] = u[ny-1,:]
		

		
		
		#dt = 0.8*dx/(10+np.max(vspeed))
		#print("TIMESTEP",dt)
		#Updating relative vorticity zeta
		#forward euler, should be leapfrog
		#Fmn = (\
		#	(u[1:ny-1,2:nx]*zeta[1:ny-1,2:nx]-u[1:ny-1,0:nx-2]*zeta[1:ny-1,0:nx-2])/(2*dx)\
		#	+(v[2:ny,1:nx-1]*zeta[2:ny,1:nx-1]-v[0:ny-2,1:nx-1]*zeta[0:ny-2,1:nx-1])/(2*dy)\
		#	+beta*v[1:ny-1,1:nx-1])
		#Fmn = Jacobian(u,v,zeta,psi,beta)
		
		#============================
		#Holton way
		#zetan(2:Nyl,1:Nxl) = zeta0(2:Nyl,1:Nxl) -beta*2*dt*v(2:Nyl,1:Nxl)...
        #    -2*dt*(dflx(2:Nyl,1:Nxl)+dfly(2:Nyl,1:Nxl))...
		
		
		dflx,dfly = divfluxHolton(psi,u,v,dx,dy,ny,nx)
		#zetanew[1:ny-1,:] = zetaold[1:ny-1,:]-2*dt*(beta*v[1:ny-1,:]+dflx[1:ny-1,:]+dfly[1:ny-1,:])
		#zetanew[1:ny-1,0:nx-1] = zetaold[1:ny-1,0:nx-1]-2*dt*(beta*v[1:ny-1,0:nx-1]+dflx[1:ny-1,0:nx-1]+dfly[1:ny-1,0:nx-1])
		
		#zetanew[1:ny-1,1:nx-1] = zetaold[1:ny-1,1:nx-1]-2*dt*(beta*v[1:ny-1,1:nx-1]+dflx[1:ny-1,1:nx-1]+dfly[1:ny-1,1:nx-1])
		
		
		#periodic bcs
		#zetanew[1:ny-1,nx-1] = zetanew[1:ny-1,0]
		#zetanew[:,nx-1] = zetanew[:,0] #HOLTON
		
		#==============================================
		#Fmn = ArakawaJacobian(psi,zeta,dx,dy,nx,ny,beta,v) Euler forward
		#Fmn = ArakawaJacobian(psi,zeta,dx,dy,nx,ny,beta,v) #Leapfrog current intermediate, actually just same


		#Euler-forward
		#zeta[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1] - dt*Fmn #for normal jacobian
		#zeta[1:ny-1,1:nx-1] = zeta[1:ny-1,1:nx-1] - dt*Fmn[1:ny-1,1:nx-1] #for Arakawa jacobian
		
		
		#Leapfrog
		#New = Old + 2*dt*current
		#zetanew[1:ny-1,1:nx-1] = zetaold[1:ny-1,1:nx-1] + 2*dt*(-Fmn[1:ny-1,1:nx-1])
		
		
		
		#Periodic Bcs
		#PROBLEM... I didn't calculate the 0th index!!!
		#Does my jacobian calculate the 0th index atm???
		#Just try with holton flux version...
		#zetanew[:,nx-1] = zetanew[:,0]
		
		#===============================================
		#Filtering zeta
		
		#ArakawaJacobian will use zetav basically, filtered, to give F(v)
		#Fmnv = ArakawaJacobian(psi,zetav,dx,dy,nx,ny,beta,v)
		#Fmn = ArakawaJacobian(psi,zeta,dx,dy,nx,ny,beta,v)
		Fmn = JacobianHolton(u,v,zeta,psi,beta,ny,nx)
		#print(Fmn)
		
		#RA-leapfrog
		#new_filter = old + 2*dt*F(current_filter)
		#vn+1 = un-1 + 2*dt*F(vn)
		#zetavnew[1:ny-1,1:nx-1] = zetaold[1:ny-1,1:nx-1] + 2*dt*(-Fmnv[1:ny-1,1:nx-1])
		
		#Indices should probably change, if doing periodic...
		#Only if bounded box, then zeta = 0 on boundaries pretty much...
		#but with periodic x... or, might actually be fine...
		#zetanew[1:ny-1,0:nx-1] = zetaold[1:ny-1,0:nx-1] + 2*dt*(-Fmn[1:ny-1,0:nx-1])#+
		#zetanew[1:ny-1,0:nx-1] = zetaold[1:ny-1,0:nx-1] + 2*dt*(-Fmn[1:ny-1,0:nx-1])
		
		#http://empslocal.ex.ac.uk/people/staff/dbs202/cat/courses/MTMW14/notes2004-d.pdf
		#eddy viscosity diffusion coefficient should be as small as possible.
		#try with 10^5 to 10^7 m^2/s
		#it is included to reduce enstropy cascade into smaller scales and avoid non-linear computational instability
		#https://en.wikipedia.org/wiki/Eddy_diffusion
		#Eddy diffusion, is diffusion due to turbulence/eddies etc... all eddies from very big gyres, to microscales. The size
		#of eddies decrease as kinetic energy decrease, until eddies are small enough for viscosity to turn it into heat
		#The microscopic processes for atmospheric mixing are too complex to model in detail, so we treat mixing with a macroscopic "eddy diffusion" processes
		#The diffusion rate can be parameterized at each pressure level, with coefficient K.
		#Called eddy diffusitivity, units m^2s^-1.
		
		
		#Is there orography induced vorticity? If i start with u[:,:] = 20, and zeta = 0? And thus, psi = 0
		#Then due to orography, do we get zeta,psi, etc
		
		nu_eddy=1000
		zetanew[1:ny-1,0:nx] = zetaold[1:ny-1,0:nx] + 2*dt*(-Fmn[1:ny-1,0:nx])
					
		zetanew[1:ny-1,1:nx-1] += -2*dt*nu_eddy*((zeta[1:ny-1,0:nx-2]-2*zeta[1:ny-1,1:nx-1]+zeta[1:ny-1,2:nx])/dx**2
						+(zeta[2:ny,1:nx-1]-2*zeta[1:ny-1,1:nx-1]+zeta[0:ny-2,1:nx-1])/dy**2)
		
		zetanew[1:ny-1,0] += -2*dt*nu_eddy*((zeta[1:ny-1,nx-1]-2*zeta[1:ny-1,0]+zeta[1:ny-1,1])/dx**2
						+(zeta[2:ny,0]-2*zeta[1:ny-1,0]+zeta[0:ny-2,0])/dy**2)
						
		zetanew[1:ny-1,nx-1] += -2*dt*nu_eddy*((zeta[1:ny-1,nx-2]-2*zeta[1:ny-1,nx-1]+zeta[1:ny-1,0])/dx**2
				+(zeta[2:ny,nx-1]-2*zeta[1:ny-1,nx-1]+zeta[0:ny-2,nx-1])/dy**2)
		
		#I chose centered differences... because the flow might come from left or right, and centered is most symmetric
		orofactor = 1
		Oroadd = np.zeros((ny,nx),dtype=np.float64)
		Oroadd[1:ny-1,1:nx-1] = -(orofactor*f[1:ny-1,1:nx-1]/H)*(u[1:ny-1,1:nx-1]*(h[1:ny-1,2:nx]-h[1:ny-1,0:nx-2])/(2*dx)+v[1:ny-1,1:nx-1]*(h[2:ny,1:nx-1]-h[0:ny-2,1:nx-1])/(2*dy))
		Oroadd[1:ny-1,0] = -(orofactor*f[1:ny-1,0]/H)*(u[1:ny-1,0]*(h[1:ny-1,1]-h[1:ny-1,nx-1])/(2*dx)+v[1:ny-1,0]*(h[2:ny,0]-h[0:ny-2,0])/(2*dy))
		Oroadd[1:ny-1,nx-1] = -(orofactor*f[1:ny-1,nx-1]/H)*(u[1:ny-1,nx-1]*(h[1:ny-1,0]-h[1:ny-1,nx-2])/(2*dx)+v[1:ny-1,nx-1]*(h[2:ny,nx-1]-h[0:ny-2,nx-1])/(2*dy))
		
		zetanew += Oroadd
		#+(2*dt*nu_eddy*((zeta[1:ny-1,0:nx-2]-2*zeta[1:ny-1,1:nx-1]+zeta[1:ny-1,:2:nx])/dx**2+(zeta[2:ny,1:nx-1]-2*zeta[1:ny-1,1:nx-1]+zeta[0:ny-2,:1:nx-1])/dy**2)
		
		#						-(f/H)*(u[1:ny-1,1:nx-1]*(h[1:ny-1,2:nx]-h[1:ny-1:,0:nx-2])/(2*dx)+v[1:ny-1,1:nx-1]*(h[2:ny,1:nx-1]-h[0:ny-2,1:nx-1])/(2*dy)))
		
		#print(u)
		#print(v)
		#Update bottom and top
		#zetanew[ny-1,:] = zetanew[ny-2,:]
		#zetanew[0,:] = zetanew[1,:]
		
		#Full slip bottom and top
		zetanew[0,:] = -(u[1,:]-u[0,:])/dy
		zetanew[ny-1,:] = -(u[ny-1,:]-u[ny-2,:])/dy
		
		#Periodic boundary conditions
		#Holton does zetan(:,Nx)=zetan(:,1);
		#zetanew[:,nx-1] = zetanew[:,0]
		#Make, better detailed version of periodic BC...
		#zetanew[1:ny-1,nx-1] = zetanew[1:ny-1,0]
		#Apply filter
		#zeta = RobertAsselinFilter(zetanew,zeta,zetaold)
		
		#No filter
		#zeta = zetanew
		
		
		#Simple filter
		#if nt == 20:
		
		#	zeta = SimpleFilter(zetanew,zeta,zetaold)
		
		
		
		
		#=========================================================
		#Output, terminal, files, etc
		
		print("u S edge", u[0,0:4])
		print("u N edge", u[-1,0:4])

		print("zeta S edge", zeta[0,0:4])
		print("zeta N edge", zeta[-1,0:4])
		print("zeta E edge", zeta[0:4,-1])
		print("zeta W edge", zeta[0:4,0])
		
		print("psi S+1 edge", psi[1,0:4])
		print("psi N-1 edge", psi[-2,0:4])
		
		
		
		vspeed = np.sqrt(u*u+v*v)
		print("max speed ", np.max(vspeed))

		#Conserved quantitity, Holton
		
		#integralFmn = np.sum(Fmn)
		#print("Fmn integral ",integralFmn)
		
		
		#Conserved quantities?
		#According to Holton,
		#Enstrophy?
		#Mean vorticity?
		#Kinetic Energy?
		#maybe only conserved in the limit as dt -> 0,
		#so, theoretically conserved, but maybe leapfrog introduces errors,
		#then you could say, it's not the Jacobian scheme error per say,
		#it's a leapfrog error, which we can't do much about in the jacobian
		

		#Mean kinetic energy should be conserved
		KEmean = np.sum(0.5*vspeed)/(ny*nx)
		Enstrophy = np.sum((zeta+f)*(zeta+f)/2)
		zetamean = np.sum(zeta)/(ny*nx)
		print("mean KE ", KEmean)
		print("Enstrophy ", Enstrophy)
		print("mean zeta ", zetamean, ", max |zeta|", np.max(np.abs(zeta)))
		
		EnstrophyList.append(Enstrophy)
		MeanKEList.append(KEmean)
		MeanZetaList.append(zetamean)
		
		
		
	#initial streamfunction
	
	#psi0 = PsiSolvePoissonJacobiChannel(psizeros,zeta0,dx,dy,nx,ny,epstol,Ngc)
	#psi0 = PsiSolvePoissonJacobiChannel(psizeros,zetaexact,dx,dy,nx,ny,epstol,Ngc)
	psi0 = FFTPoisson(nx,ny,Lx,dy,zeta_exact_initial,np.zeros_like(psi))
	
	#u0,v0 = CalcVelocity(psi0,u,v,dx,dy)
	v0,u0 = gradHolton(u,v,psi0,dx,dy,nx,ny)
	
	
	#BVE Analytical solution for Rossby wave-packets
	#k1 = 8*pi/(Lx*1000)
	#k2 = 1*pi/(Ly*1000)
	#w_dispersion = U0*k1-beta*k1/(k1**2+k2**2)
	#psi_exact = np.sin(k1*x-w_dispersion*t)*np.sin(k2*y)
	psi_exact_final = np.cos(k1*x-w*t)*np.sin(k2*y)*Amag #from https://www.researchgate.net/profile/Tarik_ALI_ZIANE/publication/237088519_The_Arakawa_Jacobian_method_and_a_fourth-order_essentially_nonoscillatory_scheme_for_the_beta-plane_barotropic_equations/links/0c96051b72e20549d2000000/The-Arakawa-Jacobian-method-and-a-fourth-order-essentially-nonoscillatory-scheme-for-the-beta-plane-barotropic-equations.pdf
	#The Rossby-Wave packets solve the equations with periodic BCs in x, and hard walls at y?
	
	zeta_exact_final = -(k1**2+k2**2)*(np.cos(k1*x-w*t)*np.sin(k2*y))*Amag
	
	#Final THESIS.pdf also has psi = a*cos(k1x-wt)sin(k2y), w=-k1/(k1^2+k2^2)
	#they solve EBVE exactly. at walls y=+-Y, they have v=0, if k2 is a multiple of pi/Y
	
	
	fmt = matplotlib.ticker.LogFormatterSciNotation()
	fmt.create_dummy_axis()
	
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
	ax00.set_xlabel('x/km')
	ax00.set_ylabel('y/km')
	plt.show()


	fig003 = plt.figure(3) #If i do pcolor, then no need for 3d projection
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	#ax003 = fig003.gca(projection='3d')
	ax003 = fig003.gca()
	#ax003.plot_surface(x, y, psi)#, rstride=3, cstride=3, color='black')
	C = ax003.contour(x/1000,y/1000,psi_exact_final,4,colors='black')
	ax003.set_title('Initial $\psi$ Rossby wave packet')
	ax003.quiver(x/1000,y/1000,u,v)
	ax003.set_xlabel('x/km')
	ax003.set_ylabel('y/km')
	plt.savefig('BVEInitialPsiRossbyWavePacket.png')
	plt.show()
	
	
	#Proxy artists for contour plots, seems official from matplotlib
	#plot initial conditions
	fig0 = plt.figure(1)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax0 = fig0.add_subplot(111)
	ax0.set_title('BVE Initial $\zeta$ configuraton')
	ax0.set_xlabel('x/km')
	ax0.set_ylabel('y/km')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax0.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	#Czeta = ax0.contour(x/1000,y/1000,zeta_exact_initial,8,colors='black',label="zeta_init")
	#Cpsi = ax0.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	#plt.legend()
	#ax0.legend()
	
	
	CS1 = ax0.contour(x/1000,y/1000,zeta_exact_initial,8,colors='black',label="zeta_init")
	#ax0.clabel(CS1, inline=1, fontsize=10)

	#CS2 = ax0.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	#ax0.clabel(CS2, inline=1, fontsize=10)

	#lines = [ CS1.collections[0], CS1.collections[-1], CS2.collections[0], CS2.collections[-1]]
	#labels = ['CS1_neg','CS1_pos','CS2_neg','CS2_pos']
	#plt.legend(lines, labels)
	
	#plt.clabel(C, fontsize=10, inline=1,fmt = '%1.0f')
	#ax0.legend([Czeta,Cpsi],["zeta_init","psi_init"])
	#ax0.legend()
	#Czeta = plt.contour(x/1000,y/1000,zeta_exact_initial,8,colors='black',label="zeta_init")
	#Cpsi = plt.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	#plt.legend()
	ax0.quiver(x/1000,y/1000,u0,v0)
	plt.savefig('BVEzetapsiinitial.png')
	plt.show()


	#plot orography/mountains
	fig001 = plt.figure() #If i do pcolor, then no need for 3d projection
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax001 = fig001.gca(projection='3d')
	ax001.plot_surface(x, y, h)#, rstride=3, cstride=3, color='black')
	ax001.set_title('Orography')
	ax001.set_xlabel('x/km')
	ax001.set_ylabel('y/km')
	ax001.set_zlabel('h/m')
	plt.savefig('BVEboundingboxorography.png')
	plt.show()

	

	
	
	fig = plt.figure(7)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax = fig.add_subplot(111)
	ax.set_title('BVE Exact $\psi$ final vs initial')
	ax.set_xlabel('x/km')
	ax.set_ylabel('y/km')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
	#Czeta = ax.contour(x/1000,y/1000,zeta,6,colors='black')
	Cpsi1 = ax.contour(x/1000,y/1000,psi_exact_final,2,colors='black',label="final")
	Cpsi2 = ax.contour(x/1000,y/1000,psi_exact_initial,2,colors='red',label="initial")
	#plt.clabel(Czeta, fontsize=10, inline=1)#,fmt = '%1.0f')
	#plt.clabel(Cpsi1, fontsize=10, inline=1)#,fmt = '%1.0f')
	#plt.clabel(Cpsi2, fontsize=10, inline=1)
	#CS1 = ax0.contour(x/1000,y/1000,zeta_exact_initial,8,colors='black',label="zeta_init")
	ax.clabel(Cpsi1, inline=1, fontsize=10, fmt=fmt)
	#CS2 = ax0.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	ax.clabel(Cpsi2, inline=1, fontsize=10, fmt=fmt)
	#lines = [ Cpsi1.collections[0], Cpsi1.collections[-1], Cpsi2.collections[0], Cpsi2.collections[-1]]
	#labels = ['Cpsi1_neg','Cpsi1_pos','Cpsi2_neg','Cpsi2_pos']
	lines = [ Cpsi1.collections[-1], Cpsi2.collections[-1]]
	labels = ['Final','Initial']
	plt.legend(lines, labels)
	#ax.quiver(x/1000,y/1000,u,v)
	#ax.legend()
	plt.savefig('BVEpsiexactfinalvsinitial.png')
	plt.show()
	
	
	fig = plt.figure(7)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax = fig.add_subplot(111)
	ax.set_title('BVE Exact $\zeta$ final vs initial')
	ax.set_xlabel('x/km')
	ax.set_ylabel('y/km')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
	#Czeta = ax.contour(x/1000,y/1000,zeta,6,colors='black')
	Cpsi1 = ax.contour(x/1000,y/1000,zeta_exact_final,2,colors='black',label="final")
	Cpsi2 = ax.contour(x/1000,y/1000,zeta_exact_initial,2,colors='red',label="initial")
	#plt.clabel(Czeta, fontsize=10, inline=1)#,fmt = '%1.0f')
	#plt.clabel(Cpsi1, fontsize=10, inline=1)#,fmt = '%1.0f')
	#plt.clabel(Cpsi2, fontsize=10, inline=1)
	#CS1 = ax0.contour(x/1000,y/1000,zeta_exact_initial,8,colors='black',label="zeta_init")
	ax.clabel(Cpsi1, inline=1, fontsize=10, fmt=fmt)
	#CS2 = ax0.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	ax.clabel(Cpsi2, inline=1, fontsize=10, fmt=fmt)
	#lines = [ Cpsi1.collections[0], Cpsi1.collections[-1], Cpsi2.collections[0], Cpsi2.collections[-1]]
	#labels = ['Cpsi1_neg','Cpsi1_pos','Cpsi2_neg','Cpsi2_pos']
	lines = [ Cpsi1.collections[-1], Cpsi2.collections[-1]]
	labels = ['Final','Initial']
	plt.legend(lines, labels)
	#ax.quiver(x/1000,y/1000,u,v)
	#ax.legend()
	plt.savefig('BVEzetaexactfinalvsinitial.png')
	plt.show()

	
	#BVE final
	#plt.switch_backend('QT4Agg')
	fig = plt.figure(2)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax = plt.gca()
	ax.set_title('BVE final $\psi$ numerical vs initial exact')
	ax.set_xlabel('x/km')
	ax.set_ylabel('y/km')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta*10**7,8,colors='black')
	#Czeta = ax.contour(x/1000,y/1000,zeta,6,colors='black')
	Cpsi1 = ax.contour(x/1000,y/1000,psi,2,colors='black')
	Cpsi2 = ax.contour(x/1000,y/1000,psi_exact_initial,2,colors='red',label="initial")
	#plt.clabel(Czeta, fontsize=10, inline=1)#,fmt = '%1.0f')
	#plt.clabel(Cpsi, fontsize=10, inline=1)#,fmt = '%1.0f')
	ax.clabel(Cpsi1, inline=1, fontsize=10, fmt=fmt)
	#CS2 = ax0.contour(x/1000,y/1000,psi0,4,colors='black',label="psi_init")
	ax.clabel(Cpsi2, inline=1, fontsize=10, fmt=fmt)
	lines = [ Cpsi1.collections[-1], Cpsi2.collections[-1]]
	labels = ['Final','Initial']
	plt.legend(lines, labels)
	ax.quiver(x/1000,y/1000,u,v)
	#ax.legend()
	plt.savefig('BVEpsifinalvsinitial.png')
	plt.show()
	
	
	#plot zeta final
	fig01 = plt.figure() #If i do pcolor, then no need for 3d projection
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	#ax01 = fig01.gca(projection='3d')
	ax01 = fig01.gca()
	#ax01.plot_surface(x, y, zeta)#, rstride=3, cstride=3, color='black')
	Czeta1 = ax01.contour(x/1000,y/1000,zeta,2,colors='black')
	Czeta2 = ax01.contour(x/1000,y/1000,zeta_exact_final,2,colors='red')
	ax01.quiver(x/1000,y/1000,u,v,color="black")
	ax01.clabel(Czeta1, inline=1, fontsize=10, fmt=fmt)
	ax01.clabel(Czeta2, inline=1, fontsize=10, fmt=fmt)
	lines = [ Czeta1.collections[-1], Czeta2.collections[-1]]
	labels = ['Final num','Final exact']
	plt.legend(lines, labels)
	ax01.set_title('$\zeta$ final numerical')
	ax01.set_xlabel('x/km')
	ax01.set_ylabel('y/km')
	plt.savefig('BVEzetafinal.png')
	plt.show()
	
	
		#plot zeta final
	fig01 = plt.figure() #If i do pcolor, then no need for 3d projection
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	#ax01 = fig01.gca(projection='3d')
	ax01 = fig01.gca()
	#ax01.plot_surface(x, y, zeta)#, rstride=3, cstride=3, color='black')
	#Czeta1 = ax01.contour(x/1000,y/1000,zeta,2,colors='black')
	Czeta2 = ax01.contour(x/1000,y/1000,zeta_exact_initial,2,colors='red')
	ax01.clabel(Czeta1, inline=1, fontsize=10)
	#ax01.clabel(Czeta2, inline=1, fontsize=10)
	#lines = [ Czeta1.collections[-1], Czeta2.collections[-1]]
	#labels = ['Final','Initial']
	#plt.legend(lines, labels)
	ax01.quiver(x/1000,y/1000,u,v,color="black")
	ax01.set_title('$\zeta$ initial Rossby wave packet')
	ax01.set_xlabel('x/km')
	ax01.set_ylabel('y/km')
	plt.savefig('BVEzetainitial.png')
	plt.show()
	
	
	#plot zeta 
	fig1, ax2 = plt.subplots(constrained_layout=True)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax2.set_title("Zeta final contourf")
	CS = ax2.contourf(x, y, zeta)
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	plt.show()
		
	#plot psi final
	fig02 = plt.figure() #If i do pcolor, then no need for 3d projection
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	#ax02 = fig02.gca(projection='3d')
	ax02 = fig02.gca()
	#ax02.plot_surface(x, y,psi)#, rstride=3, cstride=3, color='black')
	Cpsi1 = ax02.contour(x/1000,y/1000,psi,2,colors='black')
	Cpsi2 = ax02.contour(x/1000,y/1000,psi_exact_final,2,colors='red')
	ax02.clabel(Cpsi1, inline=1, fontsize=10, fmt=fmt)
	ax02.clabel(Cpsi2, inline=1, fontsize=10, fmt=fmt)
	lines = [ Cpsi1.collections[-1], Cpsi2.collections[-1]]
	labels = ['Final num','Final exact']
	plt.legend(lines, labels)
	ax02.set_title('$\psi$ final')
	ax02.set_xlabel('x/km')
	ax02.set_ylabel('y/km')
	plt.savefig('BVEpsifinal.png')
	plt.show()
	
	
	
	
	fig2 = plt.figure(3)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	ax2 = plt.gca()
	ax2.set_title('Barotropic Vorticity Equation velocity final')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	ax2.quiver(x/1000,y/1000,u,v)
	plt.show()

	#===============================================
	#Statistics
	#Should try these WITH and WITHOUT filter, to see if there's a difference...
	#
	fig30 = plt.figure(30)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	plt.plot(EnstrophyList)
	plt.title("Enstrophy")
	plt.xlabel("t/s")
	plt.savefig('BVEenstrophy.png')
	plt.show()
	
	fig31 = plt.figure(31)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	plt.plot(MeanKEList)
	plt.title("Mean Kinetic Energy")
	plt.xlabel("t/s")
	plt.savefig('BVEmeanKE.png')
	plt.show()
	
	fig32 = plt.figure(32)
	thismanager = get_current_fig_manager()
	thismanager.window.wm_geometry("-1500+0")
	plt.plot(MeanZetaList)
	plt.title("Mean Vorticity")
	plt.xlabel("t/s")
	plt.savefig('BVEMeanVorticity.png')
	plt.show()

	#https://www.ocean.washington.edu/courses/oc512/rossby-waves-gfd109-lec5a-07.pdf
	#Rossby invented beta-plane approximation
	#to explain wavy flow of westerly winds
	#They are often stationary waves, whose westward phase propagation is cancelled by zonal eastward flow
	#OFTEN, but not always
		
	#https://en.wikipedia.org/wiki/Rossby_wave
	#Rossby waves are due to rotating fluids
	#Observed both in atmosphere and oceans, due to rotation of planet
	#Atmospheric Rossby-waves are mostly high-altitude winds. The jet stream, meanders etc
	#Meanders that become very pronounced, makes cold air or warm air detach, making low-strength cyclones/anticyclones
	#Rossby waves partially explain, why eastern continental boundaries are colders, than western continental boundaries.
	#I.e, Eastern Canada, Northeast US, are colder, than western Europe, like UK, Norway, Denmark etc
	
	#https://oceanservice.noaa.gov/facts/rossby-wave.html
	#Atmospheric Rossby-waves actually primarily form, due to Earths geography??
	#
	
	
	#https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.2153-3490.1950.tb00336.x
	#Neumann brugte stereographic projection, of en plane som ligger tangent på north pole.
	#They integrate slightly larger than forecast region, to account for boundary effects propagating into the domain
	#Kan det passe, at de brugte noget med constant velocity nærmest??
	#Quasi-Geostrophic, non-divergent barotropic vorticity equation
	
	
	#NNeumann etc, de havde jo kun y som boundary med "outside domain", pga stereographic projection
	#så, det er kun ved y = 0, at de faktisk har inflow vs outflow etc
	
	
	
	
	#https://www.mi.uni-hamburg.de/arbeitsgruppen/theoretische-meteorologie/personen/lunkeit-frank/numerik/dokumente/barotrop.pdf
	#Full -slip boundary conditions
	#ah ja,
	#u = -dpsi/dy
	#vi ønsker, for full slip, u0 = u1, basically...
	#og, i gradHolten, så ved edges, der laver vi jo forward og backward differences, men ikke centered...
	#så, hvis psi[0,:] = psi[1,:] og psi[ny-1,:] = psi{ny-2,:]... then, that should work, to give u velocity at edges...
	#I can print out velocity, see if it is full-slip or not#
	#no wait... still weird... it shouldn'nt be dpsi/dy = 0, then...
	#Yes...
	#Full slip is NOT dpsi/dy = 0, i think!!!
	#it's rather, du/dy = 0.. which means, zeta = 0.. so, it's more like, zeta[0,:] = zeta[1,:] i think!
	#No, or actually just, zeta = 0, right??
	#so, I kind of need completely different BCs coded elsewhere...
	#
	
	
	
	#=====================
	#Holton script
	#Poisson solver, has psi(:,Nxl+1) = psi(:,1);
	#Which is psi(:,Nx) = psi(:,1)
	#So, it does calculate for psi at x = 0, then extends that to x = L
	
	#Velocity is supposed to be 0 at y boundaries
	#
	
	#Holton does zetan(:,Nx)=zetan(:,1);
	#In main loop, zetan periodic in x
	
	
	#I think, I need to remove, coriolis f, from zeta, before doing PoissonSolver??? Maybe???
	#nabla^2psi = zeta-f , right?
	
	#already the INITIAL config has errors, before any timesteps
	
	
	
	
	#Make symmetric calc velocity and divflux... 