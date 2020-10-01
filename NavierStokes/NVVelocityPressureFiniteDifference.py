"""
Dan Krog
based mostly on
http://www.fem.unicamp.br/~phoenics/SITE_PHOENICS/Apostilas/CFD-1_U%20Michigan_Hong/Lecture06.pdf

NV in Primitive variable Velocity/Pressure

conservations of mass, incompressible
nabla.u = 0

conservation of momentum,

du/dt = -u.nablau + alpha*nabla^2u - nablap / rho
du/dt = advection + diffusion + pressure




MAC grid, Marker-And-Cell
P pressure in the middle, they get the i,j indexs
u is to left and right of P
v is below and above P
So let's say we have a cell, then,

Pij
ui-1/2,j
ui+1/2,j
vi, j+1/2
vi, j-1/2


Oh cool, we actually change the Cell definition, based on which equation we look at
So, for X-momentum equation, the cell is centered on ui+1/2,j

FIRST TRY SIMPLE NO STAGGERED GRID


staggered grid, define, because we can't use fractional indices in the code
#define like this, and then be careful
u(i,j) = u(i+1/2,j)
v(i,j) = v(i,j+1)


#GHost Cells,
the Pressure ghost cells OUTSIDE BOUDNARY are ALWAYS = 0 I think, so no interpolation etc for them
the pressure ON BOUNDARY is not 0 though...



https://www.diva-portal.org/smash/get/diva2:1328776/FULLTEXT01.pdf
Divergence free, so constant density, so no shocks and no sound waves
But, we are still non-linear, so turbulence persists


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
	
	
	
	
	
	
	
	
	
	#TO DO
	So actually maybe also, we can't JUST do a simple double for loop both both u,v
	because one of them has nx+1 in one direction and the other has ny+1 in one direction etc...
	SO probably must be more careful...
	

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
	dt = 0.0001
	h = 1.0/(nx-1)
	dx = h
	dy = h
	
	alpha = 0.1
	rho = 1

	#time integration parameters					#hours
	time_end = 20*dt				#seconds
	
	t = 0
	nt = 0
	
	#not NEEDED, but, I think psi array should probably be nx+1,ny+1, or, velocity should be nx-1,ny-1
	# P = np.zeros((ny,nx),dtype=np.float64)
	# u = np.zeros((ny,nx),dtype=np.float64)
	# v = np.zeros((ny,nx),dtype=np.float64)
	# ut = np.zeros((ny,nx),dtype=np.float64)
	# vt = np.zeros((ny,nx),dtype=np.float64)
	# Au = np.zeros((ny,nx),dtype=np.float64) #advection
	# Av = np.zeros((ny,nx),dtype=np.float64) #advection
	# Du = np.zeros((ny,nx),dtype=np.float64) #diffusion
	# Dv = np.zeros((ny,nx),dtype=np.float64) #diffusion
	
	#Staggered Grid arrays
	P = np.zeros((ny+2,nx+2),dtype=np.float64)
	u = np.zeros((ny+2,nx+1),dtype=np.float64)
	v = np.zeros((ny+1,nx+2),dtype=np.float64)
	ut = np.zeros((ny+2,nx+1),dtype=np.float64)
	vt = np.zeros((ny+1,nx+2),dtype=np.float64)
	Au = np.zeros((ny+2,nx+1),dtype=np.float64) #advection
	Av = np.zeros((ny+1,nx+2),dtype=np.float64) #advection
	Du = np.zeros((ny+2,nx+1),dtype=np.float64) #diffusion
	Dv = np.zeros((ny+1,nx+2),dtype=np.float64) #diffusion
	
	
	
	
	nablau = np.zeros((ny,nx),dtype=np.float64)
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
	
	
	#Ghost cells
	#Ngc = 0
	
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
		#psi = PsiSolvePoissonJacobi(psi,-w,h,h,nx,ny,epstol,Ngc)
		
		
		
		#Put BCs before updating fields, probably best
		##Alternative BCs
		#ALL velocities at edges are 0, except the lid
		#no slip, and no-flow through
		#u[-1,:] = 1 #top
		#u[0,:] = 0
		#u[:,-1] = 0
		#u[:,0] = 0
		#v[-1,:] = 0
		#v[0,:] = 0
		#v[:,-1] = 0
		#v[:,0] = 0
		
		
		#STAGGERED GRID BCs, NOT ghost cells
		#no-slip, no-inflow
		# u[-2,:] = 1 #top Uwall
		# u[1,:] = 0 #bottom no slip
		# u[:,-2] = 0 #right no-inflow
		# u[:,1] = 0 #left no-inflow
		# v[-2,:] = 0 #top no-inflow
		# v[1,:] = 0 #bottom no-inflow
		# v[:,-2] = 0 #right no-slip
		# v[:,1] = 0 #left no-slip
		
		

		#STAGGERED GRID, CHOST CELLS!
		Ubottomwall = 0
		Utopwall = 1
		Urightwall = 0
		Uleftwall = 0
		Vtopwall = 0
		Vbottomwall = 0
		Vleftwall = 0
		Vrightwall = 0
		#Should not be :... should be more specific...
		#also, do reflection technique for Ghost Cells
		#u0 = 2*uwall - u1
		#ahh interesting... if we set Utopwall = 1, actually, NONE of the arrays will be set u[-1,:] = Utopwall
		#The Utopwall is actually an "inbetween" point... i.e, 
		#u[-1,:] is Ghost cell
		#u[-2,:] is NOT boundary, actually! It's ACTUALLY an INTERIOR point
		#So technically, on SOME boundaries, for SOME variables, they are actually NOT DEFINED ON THE BOUNDARY...
		#Except e.g P, I think this is defined everywhere, Ghost, Cell, AND interior...
		
		#Need to specific, WHICH u or v is defined ON which boundary??
		#And which side/variable is ONLY interior/ghost
		
		#TOP BOUNDARY:
		#u is ghost/interior,
		#v is boundary/interior
		#BOTTOM BOUNDARY:
		#u is ghost/interior
		#v is boundary/interior
		#LEFT BOUNDARY:
		#u is boundary/interior
		#v is ghost/interior
		#RIGHT BOUNDARY:
		#u is boundary/interior
		#v is ghost/interior
		
		#THe reflection technique is ONLY USED when the variable is ghost/interior, I think
		u[-1,:] = 2*Utopwall-u[-2,:] # top ghost cells
		u[0,:] = 2*Ubottomwall-u[-1,:] # #bottom ghost cell
		u[:,-1] = 0 # no inflow BOUNDARY 2*Urightwall-u[:,-2] # #right ghost
		u[:,0] = 0 #no inflow BOUNDARY #2*Uleftwall-u[:,1] #left no-inflow
		v[-1,:] = 0#no inflow BOUNDARY #2*Vtopwall-v[ #top no-inflow
		v[0,:] = 0 #no inflow BOUNDARY bottom #2*Vbottomwall-v[] #bottom no-inflow
		v[:,-1] = 2*Vrightwall-v[:,-2] #right ghost cell
		v[:,0] = 2*Vleftwall-v[:,1] #left ghost cell
		
		
		#Should be maybe ALSO.... put BCs on ut,vt???
		
		
		
		#A and D is integer j in y direction, but i+1/2 in x direction,
		#Advection
		#Right now doing it simple, take care that we really get staggered i+1/2 etc grid
		#Actually yes it's wrong..... most likely wrong, u[j,i+1] has specific meaning
		#for j in range(1,ny-1):
		#	for i in range(1, nx-1):
		#Staggered Grid, I think like this
		#We only want INTERIOR loop here.
		#Ie, NOT ghost cells, and NOT boundary cells
		for j in range(1,ny-1):
			for i in range(1, nx-1):
				
				#Put in also intermediates here actually
				#I think we do some linear interpolation actually, halfway midpoint interpolation
				ui1j = (u[j,i+1]+u[j,i])*0.5
				uij = (u[j,i]+u[j,i-1])*0.5
				uvi05j05 = 0.5*(u[j,i]+u[j+1,i])*0.5*(v[j,i]+v[j,i+1])
				uvi05jm05 = 0.5*(u[j,i]+u[j-1,i])*0.5*(v[j-1,i]+v[j-1,i+1])
				
				
				#v advection
				vij1 = (v[j+1,i]+v[j,i])*0.5
				vij = (v[j,i]+u[j-1,i])*0.5
				uvi05j05 = 0.5*(u[j,i]+u[j+1,i])*0.5*(v[j,i]+v[j,i+1]) #Actually just same as above
				uvim05j05 = 0.5*(u[j+1,i-1]+u[j,i-1])*0.5*(v[j,i]+v[j,i-1])
				
				
				Au[j,i] = (ui1j**2-uij**2)/h + (uvi05j05-uvi05jm05)/h
				Av[j,i] = (vij1**2-vij**2)/h + (uvi05j05-uvim05j05)/h
				
				#Non-staggered versions
				#Au[j,i] = (u[j,i+1]**2-u[j,i]**2)/h + (u[j+1,i+1]*v[j+1,i+1]-u[j-1,i+1]*v[j-1,i+1])/h
				#Av[j,i] = (u[j+1,i+1]*v[j+1,i+1]-u[j-1,i+1]*v[j-1,i+1])/h+(v[j+1,i]**2-v[j,i]**2)/h
				
				
		#Diffusion
		#We would actually need 2, I think?
		#D for u, and D for v?
		#for j in range(1,ny-1):
		#	for i in range(1, nx-1):
		#Staggered grid
		#We could do ONLY interior points, i.e NOT ghost cells and NOT Boundary cells,
		#but, keep it this way?
		#for j in range(1,ny):
		#	for i in range(1,nx):
		#		Du[j,i] = alpha*((u[j,i+1]-2*u[j,i]+u[j,i-1])/h**2+(u[j+1,i]-2*u[j,i]+u[j-1,i])/h**2)
		#		Dv[j,i] = alpha*((v[j,i+1]-2*v[j,i]+v[j,i-1])/h**2+(v[j+1,i]-2*v[j,i]+v[j-1,i])/h**2)
		#		for j in range(1,ny):
		#Diffusion u
		for j in range(1,ny):
			for i in range(1,nx-1):
				Du[j,i] = alpha*((u[j,i+1]-2*u[j,i]+u[j,i-1])/h**2+(u[j+1,i]-2*u[j,i]+u[j-1,i])/h**2)

		#Diffusion v
		for j in range(1,ny-1):
			for i in range(1,nx):
				Dv[j,i] = alpha*((v[j,i+1]-2*v[j,i]+v[j,i-1])/h**2+(v[j+1,i]-2*v[j,i]+v[j-1,i])/h**2)
		
		#momentum
		#u = u + dt(-A+D)+dt*nablaP / rho
		#v = v + dt(-A+D)+dt*nablaP / rho
		

		
		#splitting + project procedure... projection(fractional step) method
		#we do intermediate step, before reaching un+1
		#we do un -> ut
		#ut -> un+1
		#and some stuff changes etc
		#basically to enforce nabla.u = 0
		
		
		#FIRST INTEGRATE WITHOUT PRESSURE, THEN get PRESSURES FROM Poisson...	
		#for j in range(1,ny-1):
		#	for i in range(1,nx-1):
		#Staggered Grid, I think like thius
		#But i think, ONLY interior, right?
		#for j in range(1,ny-1):
		#	for i in range(1,nx-1):
		#		ut[j,i] = u[j,i] + dt*(-Au[j,i]+Du[j,i])
		#		vt[j,i] = v[j,i] + dt*(-Av[j,i]+Dv[j,i])
				
		#ut update
		#j ghost 0, ny+1 ghost -> interior 1, ny interior
		#i border 0, nx border -> interior 1, nx-1 interior
		for j in range(1,ny):
			for i in range(1,nx-1):
				ut[j,i] = u[j,i] + dt*(-Au[j,i]+Du[j,i])
		
		#vt update
		#j border 0, ny border -> interior 1, ny-1 interior
		#i ghost 0, nx+1 ghost -> interior 1, nx interior 
		for j in range(1,ny-1):
			for i in range(1,nx):
				vt[j,i] = v[j,i] + dt*(-Av[j,i]+Dv[j,i])
		
		
		#FIRST, I need to get the pressure from Poisson
		#THEN do projection
		#SOLVE POISSON
		npseudotime = 200
		#P[:,:] = 0
		pn = P.copy()
		for nit in range(npseudotime):
			pn = P.copy()
			#for j in range(1,ny-1):
			#	for i in range(1,nx-1):
			#Staggered Grid, I thin like this
			#ONLY INTERIOR NODES, NOT GHOST CELLS, and NOT boundary cells!
			for j in range(2,ny-1):
				for i in range(2,nx-1):
					#P[j,i] = (1/4)*(P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i])-(1/4)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
					P[j,i] = (1/4)*(pn[j,i+1]+pn[j,i-1]+pn[j+1,i]+pn[j-1,i])-(1/4)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
			
			#ALTERNATIVE PRESSURE BOUNDARIES
			#P[:,0] = P[:,1] #dp/dx = 0 at left boundary
			#P[0,:] = P[1,:] #dp/dy = 0 at bottom boundary
			#P[:,-1] = P[:,-2] #dp/dx = 0 at right boundary
			#P[-1,:] = 0 #Pressure is 0 at upper lid
			
			
			
			# #BCs from WPIdoc
			#I think the 1/3, 1/2 factor etc is because in interior points we do double centered differences
			#but at sides, we can only do centered in ONE direction, and FORWARD in another direction, so
			#We only get 1/3....
			#And at corners, we can ONLY do forward differences in BOTH directions, so it will be 1/2.... probably..
			#when evaluating the nabla^2
			
			NGC = 1
			#SIDES
			#Left side, i.e left INTERIOR
			i = 0+NGC
			#j = (0:ny)
			P[1:ny-1,i] = (1/3)*(P[1:ny-1,i+1]+P[1:ny-1,i-1]+P[2:ny,i]+P[0:ny-2,i])-(1/3)*h/dt*(ut[1:ny-1,i]-ut[j,i-1]+vt[1:ny-1,i]-vt[0:ny-2,i])
			
			#Right side, i.e right INTERIOR
			i = -1-NGC
			# #j = (0:ny)
			P[1:ny-1,i] = (1/3)*(P[1:ny-1,i+1]+P[1:ny-1,i-1]+P[2:ny,i]+P[0:ny-2,i])-(1/3)*h/dt*(ut[1:ny-1,i]-ut[j,i-1]+vt[1:ny-1,i]-vt[0:ny-2,i])
			
			#Bottom side, i.e bottom INTERIOR
			# #i = (0:nx)
			j = 0+NGC
			P[j,1:nx-1] = (1/3)*(P[j,2:nx]+P[j,0:nx-2]+P[j+1,1:nx-1]+P[j-1,1:nx-1])-(1/3)*h/dt*(ut[j,1:nx-1]-ut[j,0:nx-2]+vt[j,1:nx-1]-vt[j-1,1:nx-1])
			
			# #i = (0:nx), top INTERIOR
			j = -1-NGC
			P[j,1:nx-1] = (1/3)*(P[j,2:nx]+P[j,0:nx-2]+P[j+1,1:nx-1]+P[j-1,1:nx-1])-(1/3)*h/dt*(ut[j,1:nx-1]-ut[j,0:nx-2]+vt[j,1:nx-1]-vt[j-1,1:nx-1])
			
			
			#CORNER nodes IN the domain, not ghost chells...
			#so, i,j should actually not be 0,-1 etc. 
			#Ghost cell corners are (0,0) (0,nx+1), (ny+1,nx+1), (ny+1,0)
			#BOUNDARY corners are (1,1), (1,nx), (ny,nx), (ny,1)
			#P is actually NEVER defined ON the boundary
			#P is ONLY ghost/interior
			i = 1
			j = 1
			P[j,i] = (1/2)*(P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i])-(1/2)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
			
			i = -2
			j = -2
			P[j,i] = (1/2)*(P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i])-(1/2)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
			
			i = 1
			j = -2
			P[j,i] = (1/2)*(P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i])-(1/2)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
			
			i = -2
			j = 1
			P[j,i] = (1/2)*(P[j,i+1]+P[j,i-1]+P[j+1,i]+P[j-1,i])-(1/2)*h/dt*(ut[j,i]-ut[j,i-1]+vt[j,i]-vt[j-1,i])
			
			
			
		#PROJECTION, forward difference
		#for j in range(1,ny-1):
		#	for i in range(1,nx-1):
		#Staggered Grid method
		#ONLY INTERIOR POINTS HERE, Ghost Cells and Boundary Cells will be treated separately
		#for j in range(1,ny):
		#	for i in range(1,nx):
		#		u[j,i] = ut[j,i]-dt/(rho*h)*(P[j,i+1]-P[j,i])
		#		v[j,i] = vt[j,i]-dt/(rho*h)*(P[j+1,i]-P[j,i])
		
		#update u velocity
		for j in range(1,ny):
			for i in range(1,nx-1):
				u[j,i] = ut[j,i]-dt/(rho*h)*(P[j,i+1]-P[j,i])
		
		#update v velocity
		for j in range(1,ny-1):
			for i in range(1,nx):
				v[j,i] = vt[j,i]-dt/(rho*h)*(P[j+1,i]-P[j,i])
		
		
		#Now, more specifically, we want nabla.un+1 = 0
		#so, the updated velocity will satisfy incompressibility
		#SO... now we take the divergence of each side! Simple!
		
		#nabla.un+1 = nabla.ut -nabla*dt*nablap/rho
		#Now, enforce first time = 0, then, we get the equations,
		#0 = nabla.ut - dt/rho*nabla^2p
		
		#So basically, we are solving a poisson equation
		
		
		#Pressure boundary condition is not needed?
		#But, velocity is
		
		
		#Btw, maybe to be safe, make intermediate velocity array ut
		
		#===========================================================
		#BOUNDARY CONDITIONS
		#Reflection technique seems to be mostly done on GHOST cells actually...
		#they are done for cells e.g ny+1, so i think they are for ghost-cell purposes
		#also, you can tell what the consequence of refleciton technique is,
		#it's basically creating a velocity slightly larger than Uwall...
		#As if the velocity above the wall is even slightly larger thana t the wall and etc etc etc....
		#makes sense, in terms of viscosity etc!
		
		
		#no slip, error in slides, i think he swapped u and v
		#also, correct for ghost cells
		# #bottom boundary
		# v[:,0]= 0 
		# v[:,-1] = 0
		# U0wall = 0 
		# Unywall = 1
		# u[0,:] = 2*U0wall-u[1,:] #Bottom reflection technique
		# u[-1,:]= 2*Unywall-u[-2,:] #top reflection technique
		
		# #left boundary
		# u[:,0] = 0 #no slip
		
		
		# #right boundary
		# u[:,-1] = 0 #no slip
		
		
		# #Pressure boundary conditions, Neumann conditions from Projection
		# #These needs to be put into the Poisson solver i think actually...
		# #Corner
		# #(i,j() = (1,1) (nx,ny) (1,ny) (nx,1)
		# P[1,1] = (1/2)*beta*()
		
		
		# #vspeed = np.sqrt(u*u+v*v)
		# #vspeedmax = np.max(vspeed)
		# #dt = 0.9*0.25*h*h/(vspeedmax+0.01)
		

		
		
		
		
		print(t, dt)
		nt += 1
		t+= dt #should actually be 2*dt for leapfrog
		
		
	print("Some of P")
	print(P[:10,0:10])
	Xp = np.linspace(0,Lx,nx+2)
	Yp = np.linspace(0,Ly,ny+2)
	xp,yp = np.meshgrid(Xp,Yp)
	
	
	#plot initial poisson solution of nabla^2psi = zeta
	fig00 = plt.figure() #If i do pcolor, then no need for 3d projection
	ax00 = fig00.gca()
	#ax00.plot_surface(x, y, psi)#, rstride=3, cstride=3, color='black')
	C = ax00.contour(xp,yp,P,8,colors='black')
	ax00.set_title('Final Psi from Poisson nabla2psi = xi')
	ax00.set_xlabel('x')
	ax00.set_ylabel('y')
	plt.show()
		
	
	
	print("SIzes")
	print(xp.shape())
	print(yp.shape())
	print(u.shape())
	print(v.shape())
	fig2 = plt.figure(3)
	ax2 = plt.gca()
	ax2.set_title('Navier Stokes')
	ax2.set_xlabel('x')
	ax2.set_ylabel('y')
	#ax.quiver(x/1000,y/1000,u,v)
	#C = ax.contour(x/1000,y/1000,zeta0*10**7,8,colors='black')
	ax2.quiver(xp[1:ny,1:nx],yp[1:ny,1:nx],u[0:ny,:],v[:,0:nx])
	plt.show()


