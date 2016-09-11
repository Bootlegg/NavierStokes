from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy
#%matplotlib inline

f1 = open('u.dat','a')
f2 = open('v.dat','a')
f3 = open('p.dat','a')
f4 = open('par.dat','a')

#parameters
nx = 41
ny = 41
nt = 100#500
nit = 50
xmin = 0
xmax = 2
ymin = 0
ymax = 2
Lx = xmax-xmin
Ly = ymax-ymin

dx = Lx/(nx-1)
dy = Ly/(ny-1)

rho = 1
nu = 0.1
dt = 0.001

f4.write('{0:4.1f} {1:4.1f} {2:4.1f} {3:4.1f} {4:4.1f} {5:4.1f}'\
		.format(nt+1,dx,nx,Lx,ny,Ly))


#https://en.wikipedia.org/wiki/Reynolds_number
#reynolds number over 4000 giver turbulent?
#må prøve det... se om mit shit det virker til turbulence :P

#initialize arrays
p  = numpy.zeros((ny,nx))
b  = numpy.zeros((ny,nx))
u = numpy.zeros((ny,nx))
v = numpy.zeros((ny,nx))


x = numpy.linspace(xmin,xmax,nx)
y = numpy.linspace(ymin,ymax,ny) #y only goes to 1

X,Y = numpy.meshgrid(x,y)

#conditions
# u[-1,:] = 1
# u[0,:] = 0
# u[1:-2,0] = 0 #we don't y index to go up all the way otherwise we overwrite u[-1,:] = 0, we are at left boundary
# u[1:-2,-1] = 0 #right boundary

# v[-1,:] = 0 #v is zero at this boundary
# v[0,:] = 0
# v[1:-2,0] = 0 #we don't y index to go up all the way otherwise we overwrite u[-1,:] = 0, we are at left boundary
# v[1:-2,-1] = 0 #right boundary


# p[:,-1] = p[:,-2] #zero pressure differential at right boundary	
# p[0,:] = p[1,:] # dp/dy = 0, at y = 0, using forward differences			
# p[:,0] = p[:,1]  #zero pressure differential at left boundary
# p[-1,:] = 0 #zero pressure on that boundary, water does not push into it??	

#Write initial conditions to a file
for j in range(ny):
	for i in range(nx):
		f1.write('{0:20.5f}'.format(u[j,i]))
		f2.write('{0:20.5f}'.format(v[j,i]))
		f3.write('{0:20.5f}'.format(p[j,i]))
	f1.write("\n")
	f2.write("\n")
	f3.write("\n")

def bpoisson(b,rho,dt,u,v,dx,dy):
	#minus sign beware infront of rho
	b[1:-1,1:-1] = rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))\
					-((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2\
					-2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))\
					-((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)
					
					#this term goes up there... dunno where it's from
					#1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))
	# b[1:-1,1:-1] = rho*(1/dt*((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx)+(v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))-\
					# ((u[1:-1,2:]-u[1:-1,0:-2])/(2*dx))**2-\
					# 2*((u[2:,1:-1]-u[0:-2,1:-1])/(2*dy)*(v[1:-1,2:]-v[1:-1,0:-2])/(2*dx))-\
					# ((v[2:,1:-1]-v[0:-2,1:-1])/(2*dy))**2)
	return b

def ppoisson(p,dx,dy,b):
	pn = numpy.empty_like(p)
	pn = p.copy()
	for q in range(nit): #pseudo-time
		pn = p.copy()
		p[1:-1,1:-1] = (1.0/(2.0*(dy**2+dx**2)))*(dx**2*(pn[2:,1:-1]+pn[0:-2,1:-1])\
												+dy**2*(pn[1:-1,2:]+pn[1:-1,0:-2])\
												-b[1:-1,1:-1]*dx**2*dy**2)
												
		p[:,-1] = p[:,-2] #zero pressure differential at right boundary	
		p[0,:] = p[1,:] # dp/dy = 0, at y = 0, using forward differences			
		p[:,0] = p[:,1]  #zero pressure differential at left boundary
		p[-1,:] = 0 #zero pressure on that boundary, water does not push into it??

	return p
def CavityFlow(nt, u, v, dt, dx, dy, p, rho, nu):
	un = numpy.empty_like(u)
	vn = numpy.empty_like(v)
	b = numpy.zeros((ny, nx))
	
	for n in range(nt):
		un = u.copy()
		vn = v.copy()
		b = bpoisson(b,rho,dt,u,v,dx,dy)
		
		#Calculate pressure
		p = ppoisson(p,dx,dy,b)
		
		#east-west velocity
		u[1:-1,1:-1] = un[1:-1,1:-1]\
						-dt*un[1:-1,1:-1]*(un[1:-1,1:-1]-un[1:-1,0:-2])/dx\
						-dt*vn[1:-1,1:-1]*(un[1:-1,1:-1]-un[0:-2,1:-1])/dy\
						-(dt/rho)*(p[1:-1,2:]-p[1:-1,0:-2])/(2*dx)\
						+dt*nu*((un[1:-1,2:]-2*un[1:-1,1:-1]+un[1:-1,0:-2])/dx**2\
						+(un[2:,1:-1]-2*un[1:-1,1:-1]+un[0:-2,1:-1])/dy**2)
		#north-south velocity		
		v[1:-1,1:-1] = vn[1:-1,1:-1]\
						-dt*un[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[1:-1,0:-2])/dx\
						-dt*vn[1:-1,1:-1]*(vn[1:-1,1:-1]-vn[0:-2,1:-1])/dy\
						-(dt/rho)*(p[2:,1:-1]-p[0:-2,1:-1])/(2*dy)\
						+dt*nu*((vn[1:-1,2:]-2*vn[1:-1,1:-1]+vn[1:-1,0:-2])/dx**2\
						+(vn[2:,1:-1]-2*vn[1:-1,1:-1]+vn[0:-2,1:-1])/dy**2)
						
		#conditions
		u[0,:] = 0
		u[:,0] = 0
		u[:,-1] = 0
		u[-1,:] = 13.5

		
		#Now, in the tutorial, she doesn't care about "going all the way up"
		#u[1:-2,0] = 0 #we don't y index to go up all the way otherwise we overwrite u[-1,:] = 0, we are at left boundary
		#u[1:-2,-1] = 0 #right boundary

		v[0,:] = 0
		v[-1,:] = 0 #v is zero at this boundary
		v[:,0] = 0
		v[:,-1] = 0
		
		#Now, in the tutorial, she doesn't care about "going all the way up"
		#v[1:-2,0] = 0 #we don't y index to go up all the way otherwise we overwrite u[-1,:] = 0, we are at left boundary
		#v[1:-2,-1] = 0 #right boundary
		
		for j in range(ny):
			for i in range(nx):
				f1.write('{0:20.5f}'.format(u[j,i]))
				f2.write('{0:20.5f}'.format(v[j,i]))
				f3.write('{0:20.5f}'.format(p[j,i]))
			f1.write("\n")
			f2.write("\n")
			f3.write("\n")
		
	return u,v,p

u,v,p = CavityFlow(nt, u, v, dt, dx, dy, p, rho, nu)
fig = pyplot.figure(figsize=(11,7), dpi=100)
pyplot.contourf(X,Y,p,alpha=0.5)    ###plotting the pressure field as a contour
pyplot.colorbar()
pyplot.contour(X,Y,p)               ###plotting the pressure field outlines
pyplot.quiver(X[::2,::2],Y[::2,::2],u[::2,::2],v[::2,::2]) ##plotting velocity
pyplot.xlabel('X')
pyplot.ylabel('Y')

pyplot.show()

# fig = pyplot.figure() #If i do pcolor, then no need for 3d projection
# ax = fig.gca(projection='3d')
# #ax = fig.add_subplot(111,projection='3d')

# ax.set_xlim(0,Lx)
# ax.set_ylim(0,Ly)
# ax.set_title('Shallow Water n = 0')
# ax.set_ylabel('y [m]')
# ax.set_xlabel('x [m]')
# ax.set_zlim(dat.min(),dat.max())
# #ax.xaxis.set_ticks(0,Lx,4)
# #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
# ax.plot_surface(X, Y, dat[0], rstride=10, cstride=10)

# #line, = ax.plot_surface(X, Y, dat[0], rstride=2, cstride=2)



#contour ,cmap='jet'.... 
#pcolor looks like matlab equilevant contour
#http://matplotlib.org/examples/images_contours_and_fields/contourf_log.html
#This also looks close to Matlab
#http://stackoverflow.com/questions/15601096/contour-graph-in-python
#http://matplotlib.org/examples/pylab_examples/contourf_demo.html

# plt.pcolor(X,Y,dat[0],cmap='jet')
# plt.hold(False)
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# plt.colorbar()
# def animate(i): #i increment with 1 each step
	
	# ax.clear()
	# ax.set_zlim(dat.min(),dat.max())
	# ax.set_title('Shallow Water n = {0:3.0f}'.format(i))
	# ax.set_ylabel('y [m]')
	# ax.set_xlabel('x [m]')
	# ax.set_zlim(dat.min(),dat.max())
	# ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%5.1f'))
	# # Z = dat[i]
	# ax.plot_surface(X, Y, bottomtopo, rstride=2, cstride=2, alpha=0.3)
	# # line = ax.plot_surface(X, Y, dat[i], rstride=10, cstride=10)
	# plot=plt.contour(X,Y,dat[i],3,cmap='jet')
	# plt.pcolor(X,Y,dat[i],cmap='jet')
	# return line,

# anim = animation.FuncAnimation(fig,animate, frames = nstop, interval=500)#,blit=False)
#pyplot.show()

f1.close()
f2.close()
f3.close()
f4.close()