from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np



x0 = 0
x1 = 10
y0 = 0
y1 = 10
x = np.linspace(x0,x1,11)
y = np.linspace(y0,y1,11)




def MakeCell(Cell,a,b):
	#y = ax + b
	dx = 1

	L = np.sqrt(dx**2+a**2*dx**2)

	angle = np.arctan(a+b)
	r = np.sqrt(2)*L/2
	
	for i in range(Nx):
		j = i
		
		xc = i
		yc = a*i+b
		
		#Jeg laver sikkert ikke disse rotations ordentligt, perhaps...
		#måske skal jeg full on bruge roationsmatrix
		
		#skal faktisk rotate vectors, kan man sige, men de beholder deres end point... maybe
		#Men, de bliver også lidt længere, pga y=ax+b
		#
		
		#Vi skal også plotte normal vectors
		
		
		#Man kan loop meget bedre over en dict end jeg gør nu
		
		#Lad os se hvad vi har lige nu
		#ved 45*, for (0,0)
		#så skal x0 = 0, right?
		#og x3 = 0
		#Og x1,x2 = r, right?
		x0 = xc + r*np.cos(5*np.pi/4+angle)
		x1 = xc + r*np.cos(7*np.pi/4+angle)
		x2 = xc + r*np.cos(1*np.pi/4+angle)
		x3 = xc + r*np.cos(3*np.pi/4+angle)
		
		y0 = yc + r*np.sin(5*np.pi/4+angle)
		y1 = yc + r*np.sin(7*np.pi/4+angle)
		y2 = yc + r*np.sin(1*np.pi/4+angle)
		y3 = yc + r*np.sin(3*np.pi/4+angle)
		
		
		Cell["{}".format(i)] = {
			"Node":[i,a*i+b],
			"Corners":[[x0,y0],[x1,y1],[x2,y2],[x3,y3]],
			"V":L*L}
			
		plt.plot([x0,x1,x2,x3,x0],[y0,y1,y2,y3,y0],color="black")
		plt.scatter(xc,yc,color="black")



def MakeSidesAndNormalVectors(Cell):
	for index in Cell:
		Side1=[
		[Cell[index]["Corners"][0][0],Cell[index]["Corners"][1][0]], #x
		[Cell[index]["Corners"][0][1],Cell[index]["Corners"][1][1]]  #y
		]

		Cell[index]["Side1"] = Side1
		
		Cell[index]["Normalvector1"] = [Cell[index]["Corners"][1][1]-Cell[index]["Corners"][0][1], #dy
								-(Cell[index]["Corners"][1][0]-Cell[index]["Corners"][0][0])]
		
		
		
		Side2=[
		[Cell[index]["Corners"][1][0],Cell[index]["Corners"][2][0]], #x
		[Cell[index]["Corners"][1][1],Cell[index]["Corners"][2][1]]  #y
		]
		
		Cell[index]["Side2"] = Side2
		
		
		Cell[index]["Normalvector2"] = [Cell[index]["Corners"][2][1]-Cell[index]["Corners"][1][1],
			-(Cell[index]["Corners"][2][0]-Cell[index]["Corners"][1][0])]
		
		Side3=[
		[Cell[index]["Corners"][2][0],Cell[index]["Corners"][3][0]], #x
		[Cell[index]["Corners"][2][1],Cell[index]["Corners"][3][1]]  #y
		]

		Cell[index]["Side3"] = Side3
		
		
		Cell[index]["Normalvector3"] = [
								Cell[index]["Corners"][3][1]-Cell[index]["Corners"][2][1], #dy
								-(Cell[index]["Corners"][3][0]-Cell[index]["Corners"][2][0]) #-dx
								]
		
		Side4=[
		[Cell[index]["Corners"][3][0],Cell[index]["Corners"][0][0]], #x
		[Cell[index]["Corners"][3][1],Cell[index]["Corners"][0][1]]  #y
		]

		Cell[index]["Side4"] = Side4
		
		Cell[index]["Normalvector4"] = [
								Cell[index]["Corners"][0][1]-Cell[index]["Corners"][3][1], #dy
								-(Cell[index]["Corners"][0][0]-Cell[index]["Corners"][3][0]) #-dx
								]
	
	
def PlotNormalVectors(Cell):
	for index in Cell:
		#JEg skal btw lige displace dem, jo...
		#Lige nu bliver de plottet fra (0,0) de skal displaces mere generally til side edges
		#De skal displaces til average af corners, tror jeg
		x1 = 0.5*(Cell[index]["Corners"][0][0]+Cell[index]["Corners"][1][0])
		y1 = 0.5*(Cell[index]["Corners"][0][1]+Cell[index]["Corners"][1][1])
		
		
		plt.plot([x1,x1+Cell[index]["Normalvector1"][0]],
		[y1,y1+Cell[index]["Normalvector1"][1]],color="black")
		
		
		
		x1 = 0.5*(Cell[index]["Corners"][1][0]+Cell[index]["Corners"][2][0])
		y1 = 0.5*(Cell[index]["Corners"][1][1]+Cell[index]["Corners"][2][1])
		
		
		plt.plot([x1,x1+Cell[index]["Normalvector2"][0]],
		[y1,y1+Cell[index]["Normalvector2"][1]],color="black")
		
		
			
		x1 = 0.5*(Cell[index]["Corners"][2][0]+Cell[index]["Corners"][3][0])
		y1 = 0.5*(Cell[index]["Corners"][2][1]+Cell[index]["Corners"][3][1])
		
		
		plt.plot([x1,x1+Cell[index]["Normalvector3"][0]],
		[y1,y1+Cell[index]["Normalvector3"][1]],color="black")
		
			
		x1 = 0.5*(Cell[index]["Corners"][3][0]+Cell[index]["Corners"][0][0])
		y1 = 0.5*(Cell[index]["Corners"][3][1]+Cell[index]["Corners"][0][1])
		
		
		plt.plot([x1,x1+Cell[index]["Normalvector4"][0]],
		[y1,y1+Cell[index]["Normalvector4"][1]],color="black")
		
		#plt.plot(Cell["{}".format(i)]["Normalvector2"][0],Cell["{}".format(i)]["Normalvector2"][1],color="black")
		#plt.plot(Cell["{}".format(i)]["Normalvector3"][0],Cell["{}".format(i)]["Normalvector3"][1],color="black")
		#plt.plot(Cell["{}".format(i)]["Normalvector4"][0],Cell["{}".format(i)]["Normalvector4"][1],color="black")
	
	
	


if __name__ == "__main__":


	Cell = {}
	Nx = 5

	a = 1
	b = 0
	#Bør faktisk include x-range.... fordi det er line segments...
	MakeCell(Cell,a,b)
	MakeSidesAndNormalVectors(Cell)
	PlotNormalVectors(Cell)
	
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("pseudo-1d grid")
	plt.show()
	
	for index in Cell:
		print(Cell[index])