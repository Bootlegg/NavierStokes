import numpy as np

def SetSideHeightAndVelocity(Cell,U,v,Nx):
	"""
	We take from node center and distribute it to the sides her.
	Based on upstream etc.
	Ahh, btw...
	Jeg skal jo være sikker på 2 neighbor cells har samme enighed om hvad height og velocity er!
	Så jeg skal måske add nogle neighbor list fx til hver cell...
	
	Det virker som om jeg kommer til at gøre dobbelt så mange calculation ift. hvad der er nødvendigt...
	men, som første proof-of-concept er det måske godt nok det her, så kan jeg tænke
	
	over hvordan det skal forbedres næste gang
	SÅ lærer jeg også måderne IKKE at gøre det på
	
	
	MakeFluxDensity() skal jo kaldes EFTER denne her, så problem skal solves her i denne function
	
	Jeg skal se om jeg er nødt til at access neighbor cells...
	Jeg kan jo faktisk bare sætte både Cell[i] og neighbor cell side height på samme tid
	Det er 2 heights jeg kan sætte på samme tid, når jeg først har set hvilken side height skal have...
	
	
	En Cell node vil give sin height til højre side af sin Cell.
	Så, det er "Side2" som skal have U value fra Cell[i]
	
	Lidt weirdness med at Nx = 20, men len(Cell) = 21, så der er lidt forvirring i disse if statements
	
	
	Jeg har nu 3 functions... men jeg kan måske refactor disse loops endnu mere...
	så jeg kun har en enkelt stor function, og derefter kan kalde constnat,linear,average,quadratic...
	
	
	"""
	print("Setting SideHeight")
	for i in Cell:
		Cell[i]["Side1"][2] = 0 #we don't use side1 
		Cell[i]["Side2"][2] = v
		Cell[i]["Side3"][2] = 0 #we don't use side3
		Cell[i]["Side4"][2] = v 
		
		if int(i) == 0:
			
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			Cell["0"]["Side4"][3] = Cell["{}".format(Nx-1)]["U"]
			Cell["0"]["Side2"][3] = Cell["0"]["U"]
		
		#elif int(i) == Nx-1:
			
		
		else:
			if v > 0:
				
				#Flow goes to the right, so either we up U into side2 or side4, not sure
				Cell[i]["Side1"][3] = 0 #we don't use side1 
				Cell[i]["Side2"][3] = Cell[i]["U"] #flow goes to the right....
				Cell[i]["Side3"][3] = 0 #we don't use side3
				Cell[i]["Side4"][3] = 0 #Cell[i]["U"]
				
				if int(i) < Nx-1:
					#print("True")
					#Cell["{}".format(int(i)+1)]["Side1"][3] = 0#Cell[i]["Side2"][3]
					#Cell["{}".format(int(i)+1)]["Side2"][3] = 0#Cell[i]["Side2"][3]
					#Cell["{}".format(int(i)+1)]["Side3"][3] = 0#Cell[i]["Side2"][3]
					Cell["{}".format(int(i)+1)]["Side4"][3] = Cell[i]["Side2"][3]
					
					
					#Cell["{}".format(int(i)+1)]["Side4"][3] = Cell[i]["U"]
					
					#print("{}".format(int(i)+1))
					#print(Cell[i]["Side2"][3])
					#print(Cell["{}".format(int(i)+1)]["Side4"][3])
					#print(Cell["{}".format(int(i)+1)]["Side4"][3])
					
					print(Cell[i]["Side2"][3])
					print(Cell["{}".format(int(i)+1)]["Side4"][3])
				#Cell[i] and Cell[i+1] share side/edge/surface, so we duplicate results.
				#if int(i) < Nx-1:
					
	for i in Cell:
		if int(i) < Nx-1:
			Cell["{}".format(int(i)+1)]["Side4"][3] = Cell[i]["Side2"][3]
	# print("Printing side4 height")
	#De printer alle fucking 0.... makes no sense dude...
	#JEg sætter dem jo!!! side 4!!!!
	
	
	
	# for i in Cell:
		# print(Cell[i]["Side4"][3])
		# print(Cell["{}".format(int(i))]["Side4"][3])
		# if int(i) < Nx-1:
			# print(Cell["{}".format(int(i)+1)]["Side4"][3])
			
	# print("Cell 0 sides and heights")
	# print(Cell["0"]["Side2"][2])
	# print(Cell["0"]["Side2"][3])
	# print(Cell["0"]["Side4"][2])
	# print(Cell["0"]["Side4"][3])
def SetSideHeightAndVelocityLinear(Cell,U,v,Nx):
	"""
	Linear interpolation
	
	Kan godt være jeg mangler at gøre den "symmetric":.. lige nu laver jeg kun linear interpolation på den ene side
	Quadratic var også symmetric, right?
	Skal være sikker på, at jeg i UpdateU, gør det ordentligt.
	Tager fra begge sides.
	
	"""
	print("Setting SideHeight")
	print("NX is {}".format(Nx))
	print("Len of cell is {}".format(len(Cell)))
	for i in Cell:
		Cell[i]["Side1"][2] = 0 #we don't use side1 
		Cell[i]["Side2"][2] = v
		Cell[i]["Side3"][2] = 0 #we don't use side3
		Cell[i]["Side4"][2] = v 
		
		if int(i) == int(0):
			dx = Cell["1"]["Node"][0] - Cell["0"]["Node"][0]
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			Cell["0"]["Side4"][3] = Cell["{}".format(Nx-1)]["U"] + dx/2*(Cell["0"]["U"]-Cell["{}".format(Nx-1)]["U"])/dx
			Cell["0"]["Side2"][3] = Cell["0"]["U"] + dx/2*(Cell["1"]["U"]-Cell["0"]["U"])/dx
		
		elif int(i) == int(Nx-1):
			
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			dx = Cell["1"]["Node"][0] - Cell["0"]["Node"][0]
			#print("dx = {}".format(dx))
			#Hmmm... faktisk, så hvis man udregne dx mellem de her to, så ville dne jo give hele grid som distance!
			#Så vi må lave en fake dx på en måde... vi bruger bare node "1" og node "0"
			
			Cell[i]["Side2"][3] = Cell[i]["U"] + dx/2*(Cell["0"]["U"]-Cell[i]["U"])/dx
			
			Cell[i]["Side4"][3] = Cell["{}".format(Nx-2)]["U"] + dx/2*(Cell[i]["U"]-Cell["{}".format(Nx-2)]["U"])/dx
		
		#elif int(i) == Nx-1:
			
		
		else:
			if v > 0:
				#print(i)
				#Flow goes to the right, so either we up U into side2 or side4, not sure
				
				dx = Cell["{}".format(int(i)+1)]["Node"][0] - Cell[i]["Node"][0]
				
				Cell[i]["Side1"][3] = 0 #we don't use side1 
				Cell[i]["Side2"][3] = Cell[i]["U"] + dx/2*(Cell["{}".format(int(i)+1)]["U"]-Cell[i]["U"])/dx
				Cell[i]["Side3"][3] = 0 #we don't use side3
				Cell[i]["Side4"][3] = Cell["{}".format(int(i)-1)]["U"] + dx/2*(Cell[i]["U"]-Cell["{}".format(int(i)-1)]["U"])/dx
				

	

			
	
def SetSideHeightAndVelocityQuadratic(Cell,U,v,Nx):
	"""
	Quadratic interpolation
	"""
	print("Setting SideHeight")
	print("NX is {}".format(Nx))
	print("Len of cell is {}".format(len(Cell)))
	for i in Cell:
		Cell[i]["Side1"][2] = 0 #we don't use side1 
		Cell[i]["Side2"][2] = v
		Cell[i]["Side3"][2] = 0 #we don't use side3
		Cell[i]["Side4"][2] = v 
		
		if int(i) == int(0):
			
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			Cell["0"]["Side4"][3] = Cell["{}".format(Nx-1)]["U"]
			Cell["0"]["Side2"][3] = (1.0/8)*(
												-Cell["{}".format(Nx-1)]["U"]
												+6*Cell["0"]["U"]
												+3*Cell["1"]["U"]
												)
		
		elif int(i) == int(Nx-1):
			
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			Cell[i]["Side4"][3] = Cell["0"]["Side2"][3]
			Cell[i]["Side2"][3] = (1.0/8)*(
												-Cell["{}".format(Nx-1)]["U"]
												+6*Cell["{}".format(Nx-1)]["U"]
												+3*Cell["0"]["U"]
												)
		
		#elif int(i) == Nx-1:
			
		
		else:
			if v > 0:
				print(i)
				#Flow goes to the right, so either we up U into side2 or side4, not sure
				Cell[i]["Side1"][3] = 0 #we don't use side1 
				Cell[i]["Side2"][3] = (1.0/8)*(
												-Cell["{}".format(int(i)-1)]["U"]
												+6*Cell[i]["U"]
												+3*Cell["{}".format(int(i)+1)]["U"]
												)
				Cell[i]["Side3"][3] = 0 #we don't use side3
				Cell[i]["Side4"][3] = 0 #Cell[i]["U"]
				

	
	#==============
	#Det her bør jeg lave til en separate function som jeg kan call!
	
	#Måske først calculate "Side2", og her copy dem til opposite side...
	#Men ikke sikkert at det er så simpelt at copy dem hen, når vi bruger quadratic interpolation...
	#Eller det må det vel være? Ikke sikker :S
	for i in Cell:
		#if int(i) == int(0):
		#	Cell["0"]["Side4"][3] = Cell["{}".format(Nx)]["U"]
			
		if int(i) < Nx-1:
			Cell["{}".format(int(i)+1)]["Side4"][3] = Cell[i]["Side2"][3]


def SetSideHeightAndVelocityBothSides(Cell,U,v,Nx):
	"""
	Piecewise constant
	works for BOTH v > 0 and v < 0!
	"""
	print("Setting SideHeight")
	print("NX is {}".format(Nx))
	print("Len of cell is {}".format(len(Cell)))
	for i in Cell:
		Cell[i]["Side1"][2] = 0 #we don't use side1 
		Cell[i]["Side2"][2] = v
		Cell[i]["Side3"][2] = 0 #we don't use side3
		Cell[i]["Side4"][2] = v 
		
		if int(i) == int(0):
			
			#Side4 ved Cell 0 skal have value fra last cell, når vi har v > 0
			#Cell["0"]["Side4"] = Cell["{}".format(Nx-1)]
			#denne her skal have "U" fra Cell["Nx-1"]
			Cell["0"]["Side4"][3] = 0.5*v*(Cell["0"]["U"]+Cell["{}".format(Nx-1)]["U"]) -0.5*np.abs(v)*(Cell["0"]["U"]-Cell["{}".format(Nx-1)]["U"])
			
			Cell["0"]["Side2"][3] = 0.5*v*(Cell["1"]["U"]+Cell["0"]["U"])
			-0.5*np.abs(v)*(Cell["1"]["U"]-Cell["0"]["U"])
		
		elif int(i) == int(Nx-1):
			Cell[i]["Side4"][3] = 0.5*v*(Cell["{}".format(Nx-1)]["U"]+Cell["{}".format(Nx-2)]["U"])-0.5*np.abs(v)*(Cell["{}".format(Nx-1)]["U"]-Cell["{}".format(Nx-2)]["U"])
			Cell[i]["Side2"][3] = 0.5*v*(Cell["0"]["U"]+Cell[i]["U"])-0.5*np.abs(v)*(Cell["0"]["U"]-Cell[i]["U"])
	
		else:
			if v > 0:
				print(i)
				#Flow goes to the right, so either we up U into side2 or side4, not sure
				Cell[i]["Side1"][3] = 0 #we don't use side1 
				Cell[i]["Side2"][3] = (0.5*v*(Cell["{}".format(int(i)+1)]["U"]+Cell[i]["U"])
				-0.5*np.abs(v)*(Cell["{}".format(int(i)+1)]["U"]-Cell[i]["U"]))
				Cell[i]["Side3"][3] = 0 #we don't use side3
				Cell[i]["Side4"][3] = (0.5*v*(Cell[i]["U"]+Cell["{}".format(int(i)-1)]["U"])
				-0.5*np.abs(v)*(Cell[i]["U"]-Cell["{}".format(int(i)-1)]["U"]))
				

	
	#==============
	#Måske først calculate "Side2", og her copy dem til opposite side...
	#Men ikke sikkert at det er så simpelt at copy dem hen, når vi bruger quadratic interpolation...
	#Eller det må det vel være? Ikke sikker :S
	#for i in Cell:
		#if int(i) == int(0):
		#	Cell["0"]["Side4"][3] = Cell["{}".format(Nx)]["U"]
			
		#if int(i) < Nx-1:
		#	Cell["{}".format(int(i)+1)]["Side4"][3] = Cell[i]["Side2"][3]			
