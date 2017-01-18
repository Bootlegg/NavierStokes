from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy as np



class funcs:
	"""
	A class containing the functions/equations, derivatives and parameters.
	Perhaps superfluous, could just have seperate functions.
	"""
	
	#Number of iterations. Should change to a while loop and break from some convergence measure.
	nt = 40


	#Stepping for static gradient, the same step used for all variables
	dx = 0.001

	#increments for each functon/equation added, or updates once if all equations are added at once
	nequations = 0


	flist = np.array([])

	#def __init__(self,fu,listofdf):
	def __init__(self, listoffunctions):

		#Number of equatiomns
		funcs.nequations += len(listoffunctions)

		#Adding functions
		funcs.flist = np.append(funcs.flist, listoffunctions)
		

	def gradientf(f,x, xold = 0):
		"""
		Gradient of an n dim function, df/dx. Static gradient, constant stepping of variables, dx.
		"""
		xes1 = np.array([])
		for i in range(funcs.nequations):
			xes1 = np.append(xes1,x)
		
		xes1 = np.reshape(xes1, (funcs.nequations,funcs.nequations))
		for i in range(funcs.nequations):
			xes1[i,i] += funcs.dx


		
		dfvec = np.array([])
		for i in range(funcs.nequations):
			dfvec = np.append(dfvec, (f(xes1[i])-f(x))/funcs.dx)
		return dfvec


	def df(f,x,xold):
		"""
		Gradient of an n dim function, df/dx. Dynamic stepping of variables based on dx = x-xold

		"""

		#xes = np.array([[] for i in range(3)])
		# xes = np.array([])
		
		# x1 = [x[0]+dx,x[1]   ,x[2]   ]
		# x2 = [x[0]   ,x[1]+dx,x[2]   ]
		# x3 = [x[0]   ,x[1]   ,x[2]+dx]

		# d1fv = (f(x1)-f(x))/dx
		# d2fv = (f(x2)-f(x))/dx
		# d3fv = (f(x3)-f(x))/dx

		#dfvec = np.array([d1fv,d2fv,d3fv])
		#funcs.dfveclist = np.append(funcs.dfveclist, dfvec)
		
		dxvec = (x-xold)#np.absolute()

		xes = np.array([])
		for i in range(funcs.nequations):
			xes = np.append(xes,x)
		
		xes = np.reshape(xes, (funcs.nequations,funcs.nequations))
		for i in range(funcs.nequations):
			xes[i,i] += dxvec[i]

		dfvec = np.array([])
		for i in range(funcs.nequations):

			dfvec = np.append(dfvec, (f(xes[i])-f(xold))/dxvec[i])

		return dfvec

		#Doesnt look good:
		#dfvec, (f(xes1[i])-f(x))/dxvec[i]
		#dxvec = (x-xold)

		#Looks MAYBE good
		#dfvec, (f(xes1[i])-f(xold))/dxvec[i]
		#dxvec = (x-xold)

		#Det kan også være man kan lave en slags backwards step, så
		#f(x)-f(x-dx) / dx

	def ReshapeToJacobian(listofdfunc):
		return np.reshape(listofdfunc, (funcs.nequations,funcs.nequations))

#Really should also pass the guess to this class.
#Should make this as a call from a function, if possible, since then this script can be converted to a module.
funcs([lambda x: x[0]**2+x[1]**2-4+2*x[0]*x[1],
		lambda x: 3*x[1]+x[0],
		lambda x: 3*x[2]**2+2*x[1]*x[0]])

#Guess arrays
xold = np.zeros(funcs.nequations)+2
x = np.zeros(funcs.nequations)+10 #np.zeros(funcs.nequations)+3
xnew = np.zeros(funcs.nequations)


print("List of functions in class")
print(funcs.flist)



#We iterate the guesses using Newtons method, with Jacobi ,matrix
for n in range(funcs.nt-1):

	#This could also be called in the class
	Jacobian = np.array([])
	for i in range(funcs.nequations):
		Jacobian = np.append(Jacobian, funcs.gradientf(funcs.flist[i],x,xold))

	Jacobian = np.reshape(Jacobian, (funcs.nequations,funcs.nequations))	

	

	#Need to iterate this system of LINEAR equations myself. Non sparse matrix. Non diagonal.
	Jinv = np.linalg.inv(Jacobian)

	
	fvec = np.zeros(funcs.nequations)
	for i in range(funcs.nequations):
		fvec[i] = funcs.flist[i](x)
	
 	#Newtons Method here
	xnew = x-np.dot(np.array(Jinv),fvec)

	xold = x
	x = xnew



#iterated point
print('Final iterated point, THE SOLUTION')
print(x)

#function evaluated
print('Functions at the iterated point')

n = 0
for f in funcs.flist:
	print("f{} = {}".format(n,f(x)))
	n+=1
