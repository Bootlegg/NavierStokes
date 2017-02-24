import numpy as np
import matplotlib.pyplot as plt


xdat = np.random.uniform(0,10,15)
ydat = np.random.uniform(0,10,15)

#x = np.linspace(0,10,5)
#y = np.linspace(0,10,5)

#X,Y = np.meshgrid(x,y)


points = [[xdat[i],ydat[i]] for i in range(len(xdat))]
#points = []
#for i in range(5):
#	for j in range(5):
#		points.append([X[i,j],Y[i,j]])
#points.extend([[0,i] for i in range(10)])
#points.extend([[10,i] for i in range(10)])
#points.extend([[i,0] for i in range(10)])
#points.extend([[i,10] for i in range(10)])


circles = []
fig, ax = plt.subplots()

#Ikke helt godt med disse 3 loops her, det blvier ikke rigtigt...
#Jeg skal nok lave combinations i stedet for, faktisk...
#Så jeg skal tage alle 3 points uniquely, lave sæt af 3 points som er unique...

#Add a,b,r
Circles = []

for point1 in points:
	for point2 in points:
		for point3 in points:
			if point1 != point2:
				if point1 != point3:
					if point2 != point3:

						Sx = 0.5*np.linalg.det([[point1[0]**2+point1[1]**2,point1[1],1],
												[point2[0]**2+point2[1]**2,point2[1],1],
												[point3[0]**2+point3[1]**2,point3[1],1]])

						Sy = 0.5*np.linalg.det([[point1[0],point1[0]**2+point1[1]**2,1],
												[point2[0],point2[0]**2+point2[1]**2,1],
												[point3[0],point3[0]**2+point3[1]**2,1]])

						a = np.linalg.det([[point1[0], point1[1],1],
											[point2[0],point2[1],1],
											[point3[0],point3[1],1]
											])

						b = np.linalg.det([[point1[0],point1[1],point1[0]**2+point1[1]**2],
											[point2[0],point2[1],point2[0]**2+point2[1]**2],
											[point3[0],point3[1],point3[0]**2+point3[1]**2]
											])

						r = np.sqrt((point1[0]-Sx/a)**2+(point1[1]-Sy/a)**2)
						
						#print(Sx/a,Sy/b,r)
						# circle1 = plt.Circle((Sx/a, Sy/a), r,color='r',fill=False)
						# ax.add_artist(circle1)
						if [Sx/a,Sy/a,r] not in Circles:
							Circles.append([Sx/a,Sy/a,r,point1,point2,point3])

				#PlotCircle = []

				#xcircle1 = np.linspace(a-r,a+r,20)
				#ycircle1 = b+np.sqrt(r**2-(xcircle1-a)**2)
				#ycircle2 = b-np.sqrt(r**2-(xcircle1-a)**2)


				#for i in range(len(xcircle1)): 
				#	PlotCircle.append([xcircle1[i],ycircle1[i]])
				#	PlotCircle.append([xcircle1[i],ycircle2[i]])
				
				#for i in range(len(PlotCircle)):
				#	plt.scatter(*PlotCircle[i], color = "black",s=4)


for circle in Circles:
	PointsInside = 0
	for point in points:
		if np.sqrt((circle[0]-point[0])**2+(circle[1]-point[1])**2) < circle[2]-0.00001:
			#print(circle)
			#Break stopper kun første loop, right?
			#
			PointsInside += 1
			break
			#Hmm, problemet er her, at jeg skal tjekke for ALLE points
			#Fordi jeg kan nå at draw circle, FØR den break..

	if PointsInside == 0:
		#circle1 = plt.Circle((circle[0], circle[1]), circle[2],color='black',fill=False,alpha=0.03333)
		#ax.add_artist(circle1)
		plt.plot([circle[3][0],circle[4][0]],[circle[3][1],circle[4][1]],c="black")
		plt.plot([circle[3][0],circle[5][0]],[circle[3][1],circle[5][1]],c="black")
		plt.plot([circle[4][0],circle[5][0]],[circle[4][1],circle[5][1]],c="black")
#print(Circles)
#circle1 = plt.Circle((5,5),1)

#ax.add_artist(circle1)
ax.scatter(xdat,ydat,c="black")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Delaunay Triangulation, random points")
plt.show()