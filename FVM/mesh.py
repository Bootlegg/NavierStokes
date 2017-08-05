from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np


#Lad os lave et 2d cartesian grid, med henblik på FVM terminology
#Lad os sige fra 10x10 grid, [0,10] i x,y altså
#Lad os sige nodes er (i,j), right, så fx (4,5)
#Altså, så må selve control volumes jo være
# [0.5,1.5,0.5,1.5] til (1,1)
# [1.5,3,0.5,1.5] til (2,1)

#Alternativt, så kunne (i,j) være corner af en control volume...
#så (1,1) er første corner,
#så må første node være (0.5,0.5)
#Er måske nemmere?
#Yea, tror jeg faktisk, lad os prøve det

x0 = 0
x1 = 10
y0 = 0
y1 = 10
x = np.linspace(x0,x1,11)
y = np.linspace(y0,y1,11)



#Måske skal vi have en naming convention?
#vi har 100 nodes her
#vi label dem 1,2,3,4.....100
#en dictionary

#ControlVolumes={"1":[Node=[0.5,0.5],Corners,Edges,Faces,Volume],"2":"[Node,Corners,Edges;Faces,Volume]}
ControlVolumes={
"1":{"Node":[0.5,0.5],"Corners":[[0,0],[1,0],[1,1],[0,1]],"V":1},"2":{"Node":[1.5,1.5]}
}


#Ahh,hmm.. noget går galt med navngivning her...
#Lige nu, så er Cell["ij"] den cell som har lower left corner (i,j)
#Så makes sense, men det er nok ikke særlig general for fx triangular meshes

nx = 4
ny = 4
Cell={}
for i in range(nx):
	for j in range(ny):
	
		#Evt, "Corners":{"1":[],"2":[],"3":[],"4":[]}
		Cell["{}{}".format(i,j)] = {"Node":[i+0.5,j+0.5], "Corners":[[i,j],[i+1,j],[i+1,j+1],[i,j+1]],"V":1}
		
		#Skal lave normalvectors
		#Nooow....
		#Hvis jeg har 2 corners...
		#Så laver man en straight line til disse corners
		#dx = x2-x1
		#dy = y2-y1
		#n1 = (-dy,dx)
		#n2 = (dy,-dx)
		

		
		Cell["{}{}".format(i,j)]["NormalVectors"] = [[],[],[],[]]
		
		dx = Cell["{}{}".format(i,j)]["Corners"][1][0]-Cell["{}{}".format(i,j)]["Corners"][0][0]
		dy = Cell["{}{}".format(i,j)]["Corners"][1][1]-Cell["{}{}".format(i,j)]["Corners"][0][1]
		Cell["{}{}".format(i,j)]["NormalVectors"][0] = [dy,-dx]
		
		
		dx = Cell["{}{}".format(i,j)]["Corners"][2][0]-Cell["{}{}".format(i,j)]["Corners"][1][0]
		dy = Cell["{}{}".format(i,j)]["Corners"][2][1]-Cell["{}{}".format(i,j)]["Corners"][1][1]
		Cell["{}{}".format(i,j)]["NormalVectors"][1] = [dy,-dx]
		
		
		dx = Cell["{}{}".format(i,j)]["Corners"][3][0]-Cell["{}{}".format(i,j)]["Corners"][2][0]
		dy = Cell["{}{}".format(i,j)]["Corners"][3][1]-Cell["{}{}".format(i,j)]["Corners"][2][1]
		Cell["{}{}".format(i,j)]["NormalVectors"][2] = [dy,-dx]
		
		
		dx = Cell["{}{}".format(i,j)]["Corners"][0][0]-Cell["{}{}".format(i,j)]["Corners"][3][0]
		dy = Cell["{}{}".format(i,j)]["Corners"][0][1]-Cell["{}{}".format(i,j)]["Corners"][3][1]
		Cell["{}{}".format(i,j)]["NormalVectors"][3] = [dy,-dx]
		
		
		#Så skal jeg vel faktisk også , add en height til Cell, den skal jo være i Cell Node, vi siger bare U
		
		Cell["{}{}".format(i,j)]["U"] = np.exp(-(i-x1/2)**2-(j-y1/2)**2)
		
		#Cool. Kan jeg plot disse heights?
		






#print(ControlVolumes2)
print(x,y)

print(Cell)