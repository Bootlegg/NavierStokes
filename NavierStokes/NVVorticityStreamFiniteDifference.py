"""
Based on 
http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture05.pdf

Navier-Stokes in Vorticity-Stream formulation, Finite-Difference

w = dv/dx - du/dy

u = dpsi/dy
v = -dpsi/dx

automatically satisefies, incompressibility, non-divergence,
du/dx + dv/dy = 0


Poisson equation, Elliptic equation
nabla^2psi = -w



BCs:
no-flow across boundary, and, 
left and right,
u = 0 = dpsi/dy
psi = constant

top and bottom
v = 0 = dpsi/dx
psi = constant

because boundaries meet, same constant on all boundaries, psi = constant 


We enforce no-slip boundary conditions.
right and left,
v = 0 = -dpsi/dx
bottom boundary 
u = 0 = dpsi/dy
top boundary
u = Uwall = dpsi/dy

vorticity boundary conditions come from psi BCs
basically the Poisson equation
left and right:
w = -d^2psi/dx^2

top and bottom
w = -d^2psi/dy^2



Dan Krog

"""



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
	
	
	
	
	Should solve Poisson WITH beta plane?
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
	psin = psi.copy()
	
	#error betweeen iterated solution and RHS at each gridpoint
	eps = np.zeros((ny,nx),dtype=np.float64)
	
	#DO I NEED THIS???
	nabla2psi = np.zeros((ny,nx),dtype=np.float64)
	
	#I should, imo, be explicit and set BCs psi = 0....
	
	##JACOBI METHOD
	#dtau = 0.5*0.5*(0.5*dx**2+0.5*dy**2)
	nit = 0
	while True:
	#for r in range(500): #pseudo-time
	
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
		
		
		#Update iteration
		psi[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc] = psin[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]+eps[1+Ngc:ny-1-Ngc,1+Ngc:nx-1-Ngc]/(2/dx**2+2/dy**2)
		
		
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
		#Jeg tror, vi skal tage et note ud af Holtons bog, og her, skipper vi psin[1:ny-1,-1] og bruger psin[1:ny-1,-2]
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
		psi[:,-1] = 0 #right boundary
		psi[:,0] = 0	#left boundary
		
		
		#What if i don't set it to 0? afaik it should be constant along the edges, but maybe not 0...
		#Det har i hvert fald bestemt en effekt om man sætter det til 0 eller ej..
		#Probably not needed to put it to 0, if I'm just careful the other places,
		#But, to follow the algorithm completely, I here set the y boundarys upper and lower to 0
		psi[0,:] = 0  #south boundary
		psi[-1,:] = 0 #north boundary
		
		
		#Psi no-flow across boundary ACTUAL BOUNDARY
		#psi[:,-1+Ngc] = 0 #right boundary
		#psi[:,0+Ngc] = 0	#left boundary
		#psi[0+Ngc,:] = 0  #south boundary
		#psi[-1-Ngc,:] = 0 #north boundary
		
		#
		
		#No-slip psi streamfunction BC
		#Maybe I need additional BCs due to the no-slip condition??
		#no-slip, hence, on sourthern boundary,
		#hence, u0 = 0 on sourthern boundary, and ALSO v0 = 0 on sourthern boundary
		#now, u0 = 0, but, u1 != 0 likely...
		#so, on southern boundary, there would be, du/dy...
		#hence, there WOULD be a psi and a zeta, right?
		#ie, for no-slip, then zeta is NOT zero right, at boundary?
		#yep... for no-slip, zeta is NOT zero at boundary, and i think, psi is ALSO a bit more complicated
		#maybe,
		#psi[:,1] = psi[:,0]+u[:,0]*dy
		#something like that..
		#actually hmm it might be good enough.... i actually think it's okay,
		#psi[:,1] will come from solving poisson equation, 
		
		
		#But, should maybe plot the solution... to see if the boundaries are somewhat satisfied...
		#so make a 3d plot of the psi field...
		
		
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
		#så SKAL man jo nærmest have en staggered grid?
		#almost JUST HAPPENS to be that way...
		#because, u = dpsi/dx... so... it's like... i mean maybe you could FORCE it to be unstaggered,
		#but, by it's very nature, it seems staggered
		#so, how do we deal with that then??
		
		
		#psi[:,-1-1] = 0 #right boundary
		#psi[:,0+1] = 0	#left boundary
		#psi[0+1,:] = 0  #south boundary
		#psi[-1-1,:] = 0 #north boundary
		
		
		#Now me tihnking, ghost cells or just complicated boundary conditions?
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







def nabladotv(u,v,nx,ny,dx,dy):
	"""
	Divergence should be 0 for BVE,
	nabla.v =0
	
	Test it out
	Should maybe be centered differences??
	And maybe on edges, do forward differences
	"""
	
	
	divv = (u[:,1:nx-1]-u[:,0:nx-2])/dx + (v[1:ny-1,:]-v[0:ny-2,:])/dy

	return divv
	




def NoSlipBCvorticity(zeta,psi,dx,dy):

	#No-slip zeta BC
	#GHOST CELLS OR BOUNDARY???
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
	
	
	





def CalcVelocity(psi,dx,dy):
	"""
	Should probably also take u,v as input, or initiate u,v here...
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
	#		u[j,i] = (psi[j+1,i]-psi[j-1,i])/(2*dy)
	#		v[j,i] = -(psi[j,i+1]-psi[j,i-1])/(2*dx)
	
	u[1:ny-1,1:nx-1] = (psi[2:ny,1:nx-1]-psi[0:ny-2,1:nx-1])/(2*dy)
	v[1:ny-1,1:nx-1] = -(psi[1:ny-1,2:nx]-psi[1:ny-1,0:nx-2])/(2*dx)		
			
	return u,v

if __name__ == "__main__":

	import numpy as np
	#from pylab import *
	import matplotlib.pyplot as plt
	import matplotlib.animation as animation
	from mpl_toolkits.mplot3d import axes3d

	
	
	print("ARRAY TEST")
	atest = np.zeros((4,4))-3
	print(np.abs(atest))
	if (atest<4).all():
		
		print("ths works .all()")
		
	if (atest<4).any():
		print("ths works .any()")
	
	

	Lx = 1
	Ly = 1
	nx = 22
	ny = 22
	#x = np.linspace(0,Lx,nx)
	#y = np.linspace(0,Ly,ny)
	X = np.linspace(0,Lx,nx)
	Y = np.linspace(0,Ly,ny)
	x,y = np.meshgrid(X,Y)

	#initial vorticity and streamfunction
	
	
	#Poisson solver error
	epstol = 1e-4

	Uwall = 0.5
	Re = 1000 #Re = uL/nu  where u is flow speed, L is length, nu is kinematic viscosity
	dt = 0.001
	h = 1.0/(nx-1)
	dx = h
	dy = h

	#time integration parameters					#hours
	time_end = 5000*dt				#seconds
	
	t = 0
	nt = 0
	
	#not NEEDED, but, I think psi array should probably be nx+1,ny+1, or, velocity should be nx-1,ny-1
	psi = np.zeros((ny,nx),dtype=np.float64)
	dypsi = np.zeros((ny,nx),dtype=np.float64)
	dxpsi = np.zeros((ny,nx),dtype=np.float64)
	u = np.zeros((ny,nx),dtype=np.float64)
	v = np.zeros((ny,nx),dtype=np.float64)
	w = np.zeros((ny,nx),dtype=np.float64)
	#dfly = np.zeros((ny,nx),dtype=np.float64)
	#dflx = np.zeros((ny,nx),dtype=np.float64)
	
	
	# num1 = 125
	# num2 = 28
	# num3 = -25
	# print(num1)
	# print(num2)
	# print(num3)
	# print("AVERAGE", (num1+num2+num3)/3)
	
	# print("START")
	# print("num1 is", num1)
	# print("num2 is", num2)
	# print("num3 is", num3)
	# print("The average of num1,num2,and num3 is", (num1+num2+num3)/3)
	
	# print("DONE")
	
	
	
	#Tracer transport initial conditions
	S = np.zeros((ny,nx),dtype=np.float64)
	S0 = 0.01
	r = 0.01
	S[:,:] = S0*np.exp(-((x-0.2)**2+(y-0.95)**2)/(2*r**2))

	
	#Ghost cells
	Ngc = 0
	
	#Boundary conditions from beginning
	#zeta = NoSlipBCvorticity(zeta,psi,dx,dy)
	
	
	
	#====================================================================================
	#Main loop. Integrate forward in time
	
	
	#First step
	#Euler forward is easy to implement...
	#If you can accept the large first error?
	#Alternatively, you can use, Euler-Backwards (Matsuno) step on the first?
	#https://maths.ucd.ie/~plynch/LECTURE-NOTES/NWP-2004/NWP-CH03-2-3-P4.pdf
	#You're basically doing a leapfrog scheme with half-step too i think...
	


	while t < time_end:
		
		#Do a copy, but just assigning may also be work
		#old zeta, is equal to the intermediate zetacurrent we had in previous step

		#First solve for psi
		psi = PsiSolvePoissonJacobi(psi,-w,h,h,nx,ny,epstol,Ngc)
		
		
		
		#Put BCs before updating vorticity, probably best
		#vorticity no slip BCs
		#when it goes with the axis, it gets a minus-, and +plus when going against the axis
		#hmm, -w = d^2psi/dx^2 + d^2psi/dy^2 right? SO maybe we have double negative signs...
		#w[1:ny-1,0] = -(psi[1:ny-1,0]-psi[1:ny-1,1])*2/h**2#left i think these left and right ALSO need Uwall term!!!
		#w[1:ny-1,-1] = -(psi[1:ny-1,-1]-psi[1:ny-1,-2])*2/h**2#right
		#w[0,1:nx-1] = (psi[0,1:nx-1]-psi[1,1:nx-1])*2/h**2 #bottom
		#w[-1,1:nx-1] = -(psi[-1,1:nx-1]-psi[-2,1:nx-1])*2/h**2-Uwall*2/h#top
		
		psi[:,0] = 0
		psi[:,-1] = 0
		psi[0,:] = 0
		psi[-1,:] = 0
		
		
		##WHAT I THINK
		#sf = streamfunction in his code... hmm, recheck these...
		w[0,1:nx-1] = -(psi[1,1:nx-1]-psi[0,1:nx-1])*2/h**2 #bottom, forward difference
		w[-1,1:nx-1] = (psi[-1,1:nx-1]-psi[-2,1:nx-1])*2/h**2-Uwall*2/h#top, backward difference
		w[1:ny-1,0] = -(psi[1:ny-1,1]-psi[1:ny-1,0])*2/h**2#left i think these left and right ALSO need Uwall term!!!
		w[1:ny-1,-1] = -(psi[1:ny-1,-1]-psi[1:ny-1,-2])*2/h**2#right #backward difference
		#maybe +Uwall???
		
		
		
		#Update vorticity w... #+- flip
		w[1:ny-1,1:nx-1] = w[1:ny-1,1:nx-1]+dt*(
						+(psi[2:ny,1:nx-1]-psi[0:ny-2,1:nx-1])*(w[1:ny-1,2:nx]-w[1:ny-1,0:nx-2])/(4*h**2)
						-(psi[1:ny-1,2:nx]-psi[1:ny-1,0:nx-2])*(w[2:ny,1:nx-1]-w[0:ny-2,1:nx-1])/(4*h**2)
					+1/(Re)*(w[1:ny-1,2:nx]+w[1:ny-1,0:nx-2]+w[2:ny,1:nx-1]+w[0:ny-2,1:nx-1]-4*w[1:ny-1,1:nx-1])/h**2)

		
		
		
		#Tracer transport
		#add some chemical
		
		u,v = CalcVelocity(psi,h,h)
		#they changed u=-u and v = -v i think, different definition
		u[0,:] = (psi[1,:] - psi[0,:])/dy #south
		u[-1,:] = (psi[-1,:] - psi[-2,:])/dy #north
		v[:,0] = -(psi[:,1] - psi[:,0])/dx #west
		v[:,-1] = -(psi[:,-1] - psi[:,-2])/dx #east
		
		Dcoff = 1e-14
		S[1:ny-1,1:nx-1] = S[1:ny-1,1:nx-1] + dt*(
									-u[1:ny-1,1:nx-1]*(S[1:ny-1,2:nx]-S[1:ny-1,0:nx-2])/(2*h)
									-v[1:ny-1,1:nx-1]*(S[2:ny,1:nx-1]-S[0:nx-2,1:nx-1])/(2*h)
									+Dcoff*(S[1:ny-1,2:nx]+S[1:ny-1,0:nx-2]+S[2:ny,1:nx-1]+S[0:ny-2,1:nx-1]-4*S[1:ny-1,1:nx-1])/h**2)

		#adaptive timestep dt
		
		vspeed = np.sqrt(u*u+v*v)
		vspeedmax = np.max(vspeed)
		dt = 0.9*0.25*h*h/(vspeedmax+0.01)
		
		print(t, dt)
		nt += 1
		t+= dt #should actually be 2*dt for leapfrog
		
		
	
	#plot initial poisson solution of nabla^2psi = zeta
	fig00 = plt.figure() #If i do pcolor, then no need for 3d projection
	ax00 = fig00.gca()
	#ax00.plot_surface(x, y, psi)#, rstride=3, cstride=3, color='black')
	C = ax00.contour(x,y,psi,8,colors='black')
	ax00.set_title('Final Psi from Poisson nabla2psi = xi')
	ax00.set_xlabel('x')
	ax00.set_ylabel('y')
	plt.show()
	
	

	#plot final conditions
	fig0 = plt.figure(1)
	ax0 = plt.gca()
	ax0.set_title('Navier Stokes Final configuraton')
	ax0.set_xlabel('x')
	ax0.set_ylabel('y')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax0.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	C = ax0.contour(x,y,w,8,colors='black')
	#plt.clabel(C, fontsize=10, inline=1,fmt = '%1.0f')
	plt.savefig('NV2.png')
	plt.show()

	#plot final tracer conditions
	fig0S = plt.figure(100)
	ax0S = plt.gca()
	ax0S.set_title('Navier Stokes Final tracer configuraton')
	ax0S.set_xlabel('x')
	ax0S.set_ylabel('y')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax0.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	C = ax0S.contourf(x,y,S)
	#plt.clabel(C, fontsize=10, inline=1,fmt = '%1.0f')
	plt.savefig('NV2tracerS.png')
	plt.show()

	

	
	u,v = CalcVelocity(psi,h,h)
	#they changed u=-u and v = -v i think, different definition
	u[0,:] = (psi[1,:] - psi[0,:])/dy #south
	u[-1,:] = (psi[-1,:] - psi[-2,:])/dy #north
	v[:,0] = -(psi[:,1] - psi[:,0])/dx #west
	v[:,-1] = -(psi[:,-1] - psi[:,-2])/dx #east
	
	fig2 = plt.figure(3)
	ax2 = plt.gca()
	ax2.set_title('Navier Stokes')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	ax2.quiver(x,y,u,v)
	plt.show()


