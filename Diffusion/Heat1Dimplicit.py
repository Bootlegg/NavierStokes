import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#%reset
#%matplotlib qt

#Inspired by FDheat.pdf, problem_set_fd_implicit.pdf

#grid, conditions
nx = 200
xmin = 0
xmax = 10
L = xmax-xmin
dx = L/(nx-1) #FDheat.pdf, prolly cause we count a node at x=0, count it as 1
nt = 200
nite = 400
a = 0.15
#dt = 0.5*dx**2/a #this is approx 0.02
dt = 1

r = a*dt/dx**2
#FDheat.pdf
#alpha*dt/dx**2<1/2, gives condition on t...

x = np.linspace(xmin,xmax,nx)

p = np.zeros((nt+1,nx))
xi = np.zeros((nite+1,nx))
rx = np.zeros((nite+1,nx))
p[0,:] = np.exp(-(x-5)**2) #should learn to do the cosine hill, for better behaviour at boundary x=0,L

#This is finite difference, based on self paper
A = (1+2*r)*np.eye(nx)+(-r)*np.eye(nx,k=1)+(-r)*np.eye(nx,k=-1)
A[0,0] = 1
A[0,1] = 0
A[-1,-1] = 1
A[-1,-2] = 0
Ainv = np.linalg.inv(A)

#diagonal matrix, for inversion iteration
Adiag = (1+2*r)*np.eye(nx)
Adiag[0,0] = 1
Adiag[-1,-1] = 1
#Matrix L^-1
Adiaginv = np.linalg.inv(Adiag)

fig = plt.figure()
ax = fig.gca()
ax.set_xlim(0,L)
ax.set_ylim(0,1.1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('implicit Heat Equation,dt={}'.format(dt))
#ax.plot(x,p)
line, = ax.plot([],[])
	
def animate(i): #i increment with 1 each step

	#Here we do the actual calculation, need to change np.dot(Ainv,p)
	#first, x_1 = L^-1 * b
	xi[0,:] = np.dot(Adiaginv,p[i,:])
	#first, r_0 = b-L*x_1
	rx[0,:] = p[i,:]-np.dot(A,xi[0,:])
	#iterating x_2 = x_1 + L^-1*r_0
	for k in range(nite):
		#update x
		xi[k+1,:] = xi[k,:] + np.dot(Adiaginv,rx[k,:])
		#update r
		rx[k+1,:] = p[i,:]-np.dot(A,xi[k+1,:])

	p[i+1,:] = xi[-1,:] 
	#p[i+1,:] = np.dot(Ainv,p[i,:])
	#ax.plot(x,p)
	line.set_data(x,p[i+1,:])
	return line

anim = animation.FuncAnimation(fig,animate, frames = nt, interval=500)#,blit=False)

plt.show()
