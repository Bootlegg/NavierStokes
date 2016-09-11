from matplotlib import cm
from matplotlib import pyplot
from pylab import *
from mpl_toolkits.mplot3d import axes3d
#import numpy as np
#import matplotlib.pyplot as plt
import matplotlib.animation as animation

#%reset
#%matplotlib qt



#dat files, like fortran
u = loadtxt("u.dat")
v = loadtxt("v.dat")
p = loadtxt("p.dat")
par = loadtxt("par.dat")


nstop = int(par[0])
dx = par[1]
nx = int(par[2])
Lx = int(par[3])
ny = int(par[4])
Ly = int(par[5])

print("--------------PARAMETERS FOR PLOTTING--------------")

X = linspace(0, Lx, nx)
Y = linspace(0, Ly, ny)
X,Y = meshgrid(X,Y)

u=u.reshape(nstop,nx,ny)
v=v.reshape(nstop,nx,ny)
p=p.reshape(nstop,nx,ny)






#http://stackoverflow.com/questions/30605870/animating-quiver-in-matplotlib
fig = pyplot.figure() #If i do pcolor, then no need for 3d projection

pyplot.contourf(X,Y,p[0,:,:],alpha=0.5) 
pyplot.colorbar() #hvis jeg vil bruge pyplot.colobar(), så skal den have en image
#først, siger den som error, så jeg skal også bruge pyplot.contourf(X,Y,p[0,:,:],alpha=0.5) 
#Så jeg kan ikke bruge ax.contourf med colorbar, tror jeg...
pyplot.xlabel('x')
pyplot.ylabel('y')
pyplot.title('Driven Lid')


try:
	pyplot.contour(X,Y,p[0,:,:])###plotting the pressure field outlines    
except ValueError:  #raised if `y` is empty.
	pass

pyplot.quiver(X[::2,::2],Y[::2,::2],u[2,::2,::2],v[2,::2,::2])
#pyplot.quiver(X,Y,u[1,:,:],v[1,:,:])



def animate(i): #i increment with 1 each step, the first are 0,0,1,2,3...
	fig.clear()
	
	pyplot.contourf(X,Y,p[i,:,:],alpha=0.5)
	pyplot.colorbar()
	pyplot.quiver(X[::2,::2],Y[::2,::2],u[i,::2,::2],v[i,::2,::2])
	pyplot.xlabel('x')
	pyplot.ylabel('y')
	pyplot.title('Driven Lid')
	#Er det noget med at hun skip hver 2. arrow måske, her?
	#For ikke at draw for mange arrows, måske?

	
	#Denne gør at exception...http://stackoverflow.com/questions/22903114/overcome-valueerror-for-empty-array
	#Kan måske gøre det samme til quiver?
	try:
		pyplot.contour(X,Y,p[i,:,:])###plotting the pressure field outlines    
	except ValueError:  #raised if `y` is empty.
		pass

	return None

anim = animation.FuncAnimation(fig, animate, frames = nstop, interval=500, blit=False)	

pyplot.show()


