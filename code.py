"""
	EE2703
	Assignment 2
		Implementing spice in Python Part - 2
			(A1: Spice Part 2)
	09/02/2022
	Saravnan S V
	EE20B117
"""
import sys	# importing necessary modules			
import numpy 
from numpy import *

if (len (sys.argv)!=2):
	exit("Invalid input: Two arguments needed")	# displays an error message if the no. of arguments provided by the user is not two and exits the program

circuitfile = sys.argv[1]
if (not circuitfile.endswith ('.netlist')):
	exit ("Invalid file: Enter a .netlist file")	# displays an error message if the user entered file name is not a .netlist file and exits the program

class passive:	# class to represent passive elements in the circuit
	def __init__ (self, parameters, f=1):	# constructor to initalise class variables
		if (parameters[0][0]=='R'):
			self.impedance = float (parameters[3])	# impedance for a resistor
		if (parameters[0][0]=='L'):
			self.impedance = 2j*pi*f*float(parameters[3])	# impedance for an inductor
		if (parameters[0][0]=='C'):
			self.impedance = complex (0, -1/(2*pi*f*float(parameters[3])))	# impedance for a capacitor
		if (parameters[1][0]=='G'):
			self.nodes = [0, int(parameters[2])]
		elif (parameters[2][0]=='G'):	# assigning an alias '0' for GND
			self.nodes = [int(parameters[1]), 0]
		else:
			self.nodes = [int(parameters[1]), int(parameters[2])]	# extracting the nodes between which the element is connected

class active:	# class to represent active elements in the circuit												
	type = 0	# type 0 for current source (default)
	voltage = 0
	current = 0
	phase = 0	# represents phase of voltage/current
	nodes = [0, 0]
	def __init__ (self, parameters):
		if (parameters[3]=='ac'):	# voltage and current for AC circuit	
			self.phase = parameters[5]
			if (parameters[0][0]=='V'):
				self.type = 1	# type 1 for voltage source
				self.voltage = (float (parameters[4]))/2	# peak voltage = Vpp/2
			if (parameters[0][0]=='I'):
				self.current = (float (parameters[4]))/2	# peak current = Ipp/2
		if (parameters[3]=='dc'):	# voltage and current for DC circuit
			if (parameters[0][0]=='V'):
				self.type = 1		# type 1 for voltage source
				self.voltage = (float (parameters[4]))		# peak voltage
			if (parameters[0][0]=='I'):					 
				self.current = (float (parameters[4]))		# peak current
		if (parameters[1][0] == 'G'):
			self.nodes = [0, int (parameters[2])]			# assigning an alias '0' for GND
		elif (parameters[2][0] == 'G'):
			self.nodes = [int (parameters[1]), 0]				
		else:
			self.nodes = [int(parameters[1]), int(parameters[2])]	# extracting the nodes between which the element is connected

with open (circuitfile, "r") as f:	# opening the file in read only mode and assigning it to a file pointer f
	eachline = f.readlines()
	newline = []			# seperating the lines from the .netlist file
	for line in eachline:
		newline.append (line.split('#')[0].split('\n')[0])	# splitting and extracting ony the part before comments or \n
	try:
		start = newline.index ('.circuit')	# finding the position of .circuit
		last = newline.index ('.end')		# finding the position of .end
	except ValueError:
		exit ("Invalid circuit definition")	# checking if the file has .circuit and .end
	try:
		if (newline[last+1].split()[0]=='.ac'):
			parameters = newline[last+1].split()
			frequency = float (parameters[2])	# extracting frequency for AC circuit
	except IndexError:
		frequency = float (exp(-300))	# assigning a very low value as frequency for DC circuit
	
	passiveno = 0	# variable to store no. of passive elements; 0 by default 
	activeno = 0	# variable to store no. of active elements; 0 by default								
	passivecomponents = []	# list containing circuit parameters of passive elements 	
	activecomponents = []	# list contining circuit parameters of active elements 
	maxnode = 1	# variable to store the maximum node number; 1 by default

	for line in newline[start+1:last]:	# seperating the tokens starting from the line just before .end uptill the line below .circuit
		parameters = line.split()	# splitting the line into circuit parmeters
		count = len (parameters)
		if (count==4):													
			passivecomponents.append (passive (parameters, frequency))	# creating a class object for the passive element 
			maxnode = max (passivecomponents[passiveno].nodes[1], passivecomponents[passiveno].nodes[0], maxnode)
			passiveno += 1
		if ((count==5)or(count==6)):
			activecomponents.append (active (parameters))	# creating a class object for the active element 
			maxnode = max (activecomponents[activeno].nodes[1], activecomponents[activeno].nodes[0], maxnode)
			activeno += 1

	dimension = maxnode+activeno							
	A = array ([array([0.+0.j for k in range (dimension)]) for i in range (dimension)])
	b = array ([0.+0.j for i in range (dimension)])		# initialising the conductance and source matrices with zeroes as default

	for component in passivecomponents:	# updating the stamp of conductance matrix by analysing the passive elements one by one
			if (component.nodes[0]!=0):
				A[component.nodes[0]-1][component.nodes[0]-1] += 1/(component.impedance)
			if (component.nodes[1]!=0):
				A[component.nodes[1]-1][component.nodes[1]-1] += 1/(component.impedance)
				if (component.nodes[1]!=0):	
					A[component.nodes[0]-1][component.nodes[1]-1] -= 1/(component.impedance)
					A[component.nodes[1]-1][component.nodes[0]-1] -= 1/(component.impedance)
	
	CTSnode1 = []	# list containing positive nodes of voltage sources; (CTS: Current Through Source)
	CTSnode2 = []	# list containing negative nodes of voltage sources
	flag = maxnode
	for component in activecomponents:	# updating the source vector and stamp of conductance matrix  by analysing the active elements one by one
		if (component.type==0):		# current from positive to negative node of the voltage source is assumed to be positive 
			b[component.nodes[0]] -= component.current
			b[component.nodes[1]] += component.current
		if (component.type==1):		# voltage connected to first node is assumed to be positive 
			if (component.nodes[0]!=0):
				A[flag][component.nodes[0]-1] -= 1
				A[component.nodes[0]-1][flag] -= 1
			CTSnode1.append (component.nodes[0])
			if (component.nodes[1]!=0):
				A[flag][component.nodes[1]-1] += 1
				A[component.nodes[1]-1][flag] += 1
			CTSnode2.append (component.nodes[1])
			b[flag] = component.voltage
		flag += 1

	b = b.reshape(dimension,1)	# changing the dimension of source vector
	x = linalg.solve(A, b)		# obtaining the voltage vector						
	print ("\n")
	
	for i in range (maxnode):
		print ("Voltage at node %d is " %(i+1), end = ''); print (x[i])		# printing the node voltages 
	for i in range (0, len(CTSnode1)):
		index = maxnode+i
		print ("Current between nodes %d and %d is " %(CTSnode1[i], CTSnode2[i]), end = ''); print (x[index][0])	# printing the currents through voltage sources
  

		
			