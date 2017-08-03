"""

Use Newton's Method
Atm it doesn't work for n=1 i think.
Because of Jacobian matrix, maybe.

Dan Krog
"""

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import numpy as np
import UserInput



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
		"""
		A whole list of functions is initialized at once.
		You can also add one function at a time though, so it's flexible that way.
		"""
		#Number of equations
		funcs.nequations += len(listoffunctions)

		#Adding functions, we add to funcs.flist
		#Maybe, funcs.flist shouldn't be a numpy array, actually?
		#Numpy arrays are probably only faster if we do math on the items themselves.
		#In this case, we just call each element, because it's a function, not a number...
		#So, standard python list is probably faster?
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
		
		Returns a list?

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
		
		#Huh? Her laver jeg jo bare en copy af x, basically?
		#Så altså jeg kan vel bare sige xes = x.copy() right?
		xes = np.array([])
		for i in range(funcs.nequations):
			xes = np.append(xes,x)
		
		xes = np.reshape(xes, (funcs.nequations,funcs.nequations))
		for i in range(funcs.nequations):
			xes[i,i] += dxvec[i]
		
		
		
		#Here, we can probably remove the append way, like i did for Jacobian
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
# funcs([
		# lambda x: x[0]**2+x[1]**2-4+2*x[0]*x[1],
		# lambda x: 3*x[1]+x[0],
		# lambda x: 3*x[2]**2+2*x[1]*x[0]
		# ])

# #Guess arrays
# xold = np.zeros(funcs.nequations)+2
# x = np.zeros(funcs.nequations)+10 #np.zeros(funcs.nequations)+3
# xnew = np.zeros(funcs.nequations)


print("List of functions in class")
print(funcs.flist)


#I want this to be a function, actually
#We iterate the guesses using Newtons method, with Jacobi ,matrix
def NewtonMethod(xold,x,xnew):
	for n in range(funcs.nt-1):

		#This could also be called in the class
		#Jacobian = np.array([])
		#for i in range(funcs.nequations):
		
			#To uzywa static gradient
			#Jacobian = np.append(Jacobian, funcs.gradientf(funcs.flist[i],x,xold))
			
			
			#To uzywa dynamic gradient
			#dflist = funcs.df(funcs.flist[i],x,xold)
			#print(len(dflist))
			#Jacobian = np.append(Jacobian, dflist)

		#Jacobian = np.reshape(Jacobian, (funcs.nequations,funcs.nequations))	

		
		#Kan man gøre det uden Numpy append?
		#Jacobian = np.zeros(funcs.nequations*funcs.nequations)
		
		#for i in range(funcs.nequations):
		#	Jacobian[i*funcs.nequations:(i+1)*funcs.nequations] = funcs.df(funcs.flist[i],x,xold)
			
		#Jacobian = np.reshape(Jacobian, (funcs.nequations,funcs.nequations))


		#Kan man gøre det uden reshape? YEs!
		Jacobian = np.zeros((funcs.nequations,funcs.nequations))
		
		for i in range(funcs.nequations):
			Jacobian[i,:] = funcs.df(funcs.flist[i],x,xold)
		
		

		

		#Need to iterate this system of LINEAR equations myself. Non sparse matrix. Non diagonal.
		Jinv = np.linalg.inv(Jacobian)

		
		fvec = np.zeros(funcs.nequations)
		for i in range(funcs.nequations):
			
			fvec[i] = funcs.flist[i](x)
		#fvec[:] = funcs.flist[:](x)
		
		
		
		#Newtons Method here
		xnew = x-np.dot(Jinv,fvec) #Før var det x-np.dot(np.array(Jinv)*fvec), men Jinv er vel allerede np array?
		#xnew = x-Jinv*fvec

		xold = x
		x = xnew
	
	PrintSolution(x)
	
def PrintSolution(x):
	print('Final iterated point, THE SOLUTION')
	print(x)
	
	print('Functions at the iterated point')

	n = 0
	for f in funcs.flist:
		print("f{} = {}".format(n,f(x)))
		n+=1


if __name__=="__main__":
	
	
	
	#Here i import functions 
	#funcs([imported.list])
	#funcs([
	#	lambda x: x[0]**2+x[1]**2-4+2*x[0]*x[1],
	#	lambda x: 3*x[1]+x[0],
	#	lambda x: 3*x[2]**2+2*x[1]*x[0]
	#	])
	
	funcs(UserInput.ListOfFunctions)
	#Here i make guesses for newtons method
	#Guess arrays
	#Actually, xold need to be guessed, but, x doesn't need to be guessed?
	#Or maybe it does, because we use numerical derivatives... dunno
	xold = np.zeros(funcs.nequations)+2
	x = np.zeros(funcs.nequations)+10 #np.zeros(funcs.nequations)+3
	xnew = np.zeros(funcs.nequations)
	
	print("Initial guess")
	print(xold)
	
	# funcs([
	# lambda x: x[0]*x[1] + 3*x[0] - 20,
	# lambda x: x[0]*x[0]-10+4*x[1]
	# ])
	# xold = np.zeros(funcs.nequations)+2
	# x = np.zeros(funcs.nequations)+9 #np.zeros(funcs.nequations)+3
	# xnew = np.zeros(funcs.nequations)
	
	NewtonMethod(xold,x,xnew)
	
	
	#iterated point
	#PrintSolution(x)

	#function evaluated
	#print('Functions at the iterated point')

	#n = 0
	#for f in funcs.flist:
	#	print("f{} = {}".format(n,f(x)))
	#n+=1
	