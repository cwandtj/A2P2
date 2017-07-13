#!/usr/bin/python

from sys import argv
import json

file = argv[1]
object = argv[2]
molecule = argv[3]
geometry = argv[4]
distance = argv[5]
bond_length = argv[6]
auto = argv[7]
	  
# next following argv are index of coordinates in POSCAR 

def surface_analyzer(file,object,molecule,geometry,distance,bond_length,auto,*argv):

###########################################
# Reading information of atom from POSCAR #
###########################################

	line_list = [line.strip() for line in open(file)]  
	compact_line = [x for x in line_list if x != []]
	 
	for line in compact_line:								  # pinpoint where direct is & check Selective dynamics option
		if line.lower() == "direct" or line.lower() == "cartesian":
			coordi_type = line.lower()
			coordi_type_index = compact_line.index(line)
			coordi_start = compact_line.index(line) + 1		# line after line describing coordinate system type
	#print coordi_start 
   
	if coordi_start == 8:
		atom_species_index = coordi_start - 3
		atom_species = compact_line[coordi_start - 3]
		num_of_atom_index = coordi_start -2
		num_of_atom = compact_line[coordi_start - 2]
		#print "Selective dynamics not included"
	
	else: 
		atom_species_index = coordi_start - 4
		atom_species = compact_line[coordi_start - 4]
		num_of_atom_index = coordi_start - 3
		num_of_atom = compact_line[coordi_start - 3]
		#print "OPTION :: 'Selective dynamics' included"

	total_atom = sum(int(i) for i in num_of_atom.split())		  

	coordi_line = compact_line[coordi_start:coordi_start + total_atom]

############################################
# Direct Coodinate to Cartesian Coordinate #
############################################
 
	lat_vec = compact_line[2:5]	 
	real_line = []	


	if coordi_type == "direct":				  # direct to cartesian
		for line in coordi_line:
			splited_line = line.split()
			x_fact = float(splited_line[0]) * float(lat_vec[0].split()[0]) + float(splited_line[1]) * float(lat_vec[1].split()[0]) + float(splited_line[2]) * float(lat_vec[2].split()[0])
			y_fact = float(splited_line[0]) * float(lat_vec[0].split()[1]) + float(splited_line[1]) * float(lat_vec[1].split()[1]) + float(splited_line[2]) * float(lat_vec[2].split()[1])
			z_fact = float(splited_line[0]) * float(lat_vec[0].split()[2]) + float(splited_line[1]) * float(lat_vec[1].split()[2]) + float(splited_line[2]) * float(lat_vec[2].split()[2])
		 
			x_coordi = x_fact  *  float(compact_line[1])
			y_coordi = y_fact  *  float(compact_line[1])
			z_coordi = z_fact  *  float(compact_line[1])
			
			splited_line[0] = str(x_coordi)
			splited_line[1] = str(y_coordi)
			splited_line[2] = str(z_coordi)
			
			real_line.append(' '.join(splited_line))

	else:
		real_line = coordi_line
		  


	#print "Total number of atoms in this structure : ",total_atom	  


###################################################################
# 1. Sorting all coordinates by z-axis							#
# 2. Find the atom on the top and its z coordinate				#
# 3. Define if an atom is closer to the top atom than 1 angstrom  #
###################################################################


	sorted_coordi_line = sorted(real_line, key=lambda x: float(x.split()[2]))

	#print sorted_coordi_line

	top_coordi = sorted_coordi_line[-1].split()[2]	#  Z coordinates of an atom on the top 

	#print top_coordi

	surface_atoms = [atoms for atoms in sorted_coordi_line if float(top_coordi) - float(atoms.split()[2]) < 1.0]

	#print "\n--------------------------\n"
	#print "Surface atoms\n"
	#for i in surface_atoms:
	#	print i 

	#print "\n--------------------------"
	#print "Total number of atoms on the surface : ", len(surface_atoms) 


###############################################################################################
# Choosing absobate type and Writing its coodinates, considering input distance & bond length #
###############################################################################################

	absorbate_line = [] 
	oxygen_line = []
	hydrogen_line = []

								 ## for O, O2, H, H2 ##

	if object == "H" or object == "O":

		if object == "H":
			object_type = " H "

		elif object == "O":
			object_type = " O "


		if auto == "auto":
		
			#print "\n--------------------\n Auto mode :: Creating absorbate .... \n------------------\n"	 
	   
			if molecule.lower() == "atom": 
				
				for i in surface_atoms:
					a = i.split()
					b = float(a[2]) + float(distance) 
				   	
					a[2] = str(b)
					c = ' '.join(a)
					absorbate_line.append(c)
		
			elif molecule.lower() == "molecule":
		
				if geometry == "perpendicular":
		
					for i in surface_atoms:
						a = i.split()
						b = float(a[2]) + float(distance)
						c = float(a[2]) + float(distance) + float(bond_length)
					   
						a[2] = str(b)
						d = ' '.join(a)
						absorbate_line.append(d)
		
						a[2] = str(c)
						e = ' '.join(a)
						absorbate_line.append(e)
				 
				elif geometry == "parallel":
		
					for i in surface_atoms:
						a = i.split()
						centroid = float(a[2]) + float(distance)
						y_minus = float(a[1]) - float(bond_length)/2
						y_plus = float(a[1]) + float(bond_length)/2
		
						a[1] = str(y_minus)
						a[2] = str(centroid)
						d = ' '.join(a)
						absorbate_line.append(d)
		
						a[1] = str(y_plus)
						a[2] = str(centroid)
						e = ' '.join(a)
						absorbate_line.append(e)
		
		elif auto != "auto":
			#print "\n-------------------\n Customizing mode :: Checking coorinates ....... \n----------------------\n "
		
#			selected = int(auto)			 #  Get the coordinates in the form of list
		
#			if auto in surface_atoms:		  
#				print " You have selected one of the atoms on the surface"
#			else:
#				print " Warning !! check the atom you have selected, it may not be on the surface !! "
		
			if molecule.lower() == "atom":
					
				for i in argv[8:]:
					#print " POSCAR INDEX : ", i
					
					a = real_line[int(i)].split()
					b = float(a[2]) + float(distance)
		
					a[2] = str(b)
					c = ' '.join(a)
					absorbate_line.append(c)
		
			elif molecule.lower() == "molecule":
			   
				if geometry == "perpendicular":
		
					for i in argv[8:]:
						#print " POSCAR INDEX : ", i
						a = real_line[int(i)].split()
						b = float(a[2]) + float(distance)
						c = float(a[2]) + float(distance) + float(bond_length)
		
						a[2] = str(b)
						d = ' '.join(a)
						absorbate_line.append(d)
		
						a[2] = str(c)
						e = ' '.join(a)
						absorbate_line.append(e)
		
				elif geometry == "parallel":
		
					for i in argv[8:]:
						#print " POSCAR INDEX : ", i
						a = real_line[int(i)].split()
						centroid = float(a[2]) + float(distance)
						y_minus = float(a[1]) - float(bond_length)/2
						y_plus = float(a[1]) + float(bond_length)/2
		
						a[1] = str(y_minus)
						a[2] = str(centroid)
						d = ' '.join(a)
						absorbate_line.append(d)
		
						a[1] = str(y_plus)
						a[2] = str(centroid)
						e = ' '.join(a)
						absorbate_line.append(e)


								## for OH ##

	elif object == "OH":
	   
		object_type = " O H "

		if auto == "auto":

			molecule = "molecule"

			if geometry == "perpendicular":

				for i in surface_atoms:
					a = i.split()
					b = float(a[2]) + float(distance)
					c = float(a[2]) + float(distance) + float(bond_length)

					a[2] = str(b)
					d = ' '.join(a)
					oxygen_line.append(d)

					a[2] = str(c)
					e = ' '.join(a)
					hydrogen_line.append(e)

			elif geometry == "parallel":

				for i in surface_atoms:
					a = i.split()
					centroid = float(a[2]) + float(distance)
					y_minus = float(a[1]) - float(bond_length)/2
					y_plus = float(a[1]) + float(bond_length)/2

					a[1] = str(y_minus)
					a[2] = str(centroid)
					d = ' '.join(a)
					oxygen_line.append(d)

					a[1] = str(y_plus)
					a[2] = str(centroid)
					e = ' '.join(a)
					hydrogen_line.append(e)

		elif auto != "auto":
			#print "\n-------------------\n Customizing mode :: Checking coorinates ....... \n----------------------\n "

#			selected = int(auto)			 #  Get the coordinates in the form of list

#			if auto in surface_atoms:
#				print " You have selected one of the atoms on the surface"
#			else:
#				print " Warning !! check the atom you have selected, it may not be on the surface !! "


			molecule = "molecule"

			if geometry == "perpendicular":

				for i in argv[8:]:
					#print " POSCAR INDEX : ", i
					a = real_line[int(i)].split()
					b = float(a[2]) + float(distance)
					c = float(a[2]) + float(distance) + float(bond_length)

					a[2] = str(b)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[2] = str(c)
					e = ' '.join(a)
					hydrogen_line.append(e)

			elif geometry == "parallel":

				for i in argv[8:]:
					#print " POSCAR INDEX : ", i
					a = real_line[int(i)].split()
					centroid = float(a[2]) + float(distance)
					y_minus = float(a[1]) - float(bond_length)/2
					y_plus = float(a[1]) + float(bond_length)/2

					a[1] = str(y_minus)
					a[2] = str(centroid)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[1] = str(y_plus)
					a[2] = str(centroid)
					e = ' '.join(a)
					hydrogen_line.append(e)

								## for OOH ##

	elif object == "OOH":

		object_type = " O H "

		if auto == "auto":

			molecule = "molecule"

			if geometry == "perpendicular":

				for i in surface_atoms:
					a = i.split()
					b = float(a[2]) + float(distance)
					c = float(a[2]) + float(distance) + float(bond_length)
					d = c + float(bond_length)						 

					a[2] = str(b)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[2] = str(c)
					e = ' '.join(a)
					oxygen_line.append(e)


					a[2] = str(d)
					e = ' '.join(a)
					hydrogen_line.append(e)

			elif geometry == "parallel":

				for i in surface_atoms:
					a = i.split()
					centroid = float(a[2]) + float(distance)
					y_minus = float(a[1]) - float(bond_length)/2
					y_plus = float(a[1]) + float(bond_length)/2
					hydro = float(a[1]) + float(bond_length)/2 + float(bond_length)

					a[1] = str(y_minus)
					a[2] = str(centroid)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[1] = str(y_plus)
					a[2] = str(centroid)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[1] = str(hydro)
					a[2] = str(centroid)
					e = ' '.join(a)
					hydrogen_line.append(e)
				   


		elif auto != "auto":
			#print "\n-------------------\n Customizing mode :: Checking coorinates ....... \n----------------------\n "

#			selected = int(auto)			 #  Get the coordinates in the form of list

#			if auto in surface_atoms:
#				print " You have selected one of the atoms on the surface"
#			else:
#				print " Warning !! check the atom you have selected, it may not be on the surface !! "


			molecule = "molecule"

			if geometry == "perpendicular":

				for i in argv[8:]:
					#print " POSCAR INDEX : ", i
					a = real_line[int(i)].split()
					b = float(a[2]) + float(distance)
					c = float(a[2]) + float(distance) + float(bond_length)
					d = float(a[2]) + float(distance) + 2*float(bond_length)
					
					a[2] = str(b)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[2] = str(c)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[2] = str(d)
					e = ' '.join(a)
					hydrogen_line.append(e)


			elif geometry == "parallel":

				for i in argv[8:]:
					#print " POSCAR INDEX : ", i
					a = real_line[int(i)].split()
					centroid = float(a[2]) + float(distance)
					y_minus = float(a[1]) - float(bond_length)/2
					y_plus = float(a[1]) + float(bond_length)/2
					hydro = float(a[1]) + float(bond_length)/2 + float(bond_length)						

					a[1] = str(y_minus)
					a[2] = str(centroid)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[1] = str(y_plus)
					a[2] = str(centroid)
					e = ' '.join(a)
					oxygen_line.append(e)

					a[1] = str(hydro)
					a[2] = str(centroid)
					e = ' '.join(a)
					hydrogen_line.append(e)


								## for H2O ##

	elif object == "H2O":

		object_type = " O H "

		if auto == "auto":

			molecule = "molecule"

			geometry == "perpendicular"

			for i in surface_atoms:
				a = i.split()
				b = float(a[2]) + float(distance)
				c = float(a[2]) + float(distance) + 0.7071*float(bond_length)
				y_minus = float(a[1]) - 0.7071*float(bond_length)
				y_plus = float(a[1]) + 0.7071*float(bond_length)


				a[2] = str(b)
				d = ' '.join(a)
				oxygen_line.append(d)

				a[1] = str(y_minus)
				a[2] = str(c)
				d = ' '.join(a)
				hydrogen_line.append(d)

				a[1] = str(y_plus)
				a[2] = str(c)
				e = ' '.join(a)
				hydrogen_line.append(e)




		elif auto != "auto":
			#print "\n-------------------\n Customizing mode :: Checking coorinates ....... \n----------------------\n "

#			selected = int(auto)			 #  Get the coordinates in the form of list

#			if auto in surface_atoms:
#				print " You have selected one of the atoms on the surface"
#			else:
#				print " Warning !! check the atom you have selected, it may not be on the surface !! "


			molecule = "molecule"

			geometry == "perpendicular"

			for i in argv[8:]:
				 a = real_line[int(i)].split()
				 b = float(a[2]) + float(distance)
				 c = float(a[2]) + float(distance) + 0.7071*float(bond_length)
				 y_minus = float(a[1]) - 0.7071*float(bond_length)
				 y_plus = float(a[1]) + 0.7071*float(bond_length)
 
 
				 a[2] = str(b)
				 d = ' '.join(a)
				 oxygen_line.append(d)
 
				 a[1] = str(y_minus)
				 a[2] = str(c)
				 d = ' '.join(a)
				 hydrogen_line.append(d)
 
				 a[1] = str(y_plus)
				 a[2] = str(c)
				 e = ' '.join(a)
				 hydrogen_line.append(e)

								## for H2O2 ##

	elif object == "H2O2":

		object_type = " O H "

		if auto == "auto":

			molecule = "molecule"

			geometry == "perpendicular"

			for i in surface_atoms:
				a = i.split()
				b = float(a[2]) + float(distance)
				c = float(a[2]) + float(distance) + 0.7071*float(bond_length)
				y_minus = float(a[1]) - 0.7071*float(bond_length)
				y_plus = float(a[1]) + float(bond_length)
				y_plusplus = y_plus + 0.7071*float(bond_length)


				a[2] = str(b)
				d = ' '.join(a)
				oxygen_line.append(d)
			
				a[1] = str(y_plus)
				a[2] = str(b)
				d = ' '.join(a)
				oxygen_line.append(d)
			

				a[1] = str(y_minus)
				a[2] = str(c)
				d = ' '.join(a)
				hydrogen_line.append(d)

				a[1] = str(y_plusplus)
				a[2] = str(c)
				e = ' '.join(a)
				hydrogen_line.append(e)




		elif auto != "auto":
			#print "\n-------------------\n Customizing mode :: Checking coorinates ....... \n----------------------\n "

#			selected = int(auto)			 #  Get the coordinates in the form of list

#			if auto in surface_atoms:
#				print " You have selected one of the atoms on the surface"
#			else:
#				print " Warning !! check the atom you have selected, it may not be on the surface !! "


			molecule = "molecule"

			geometry == "perpendicular"

			for i in argv[8:]:
				a = real_line[int(i)].split()
				b = float(a[2]) + float(distance)
				c = float(a[2]) + float(distance) + 0.7071*float(bond_length)
				y_minus = float(a[1]) - 0.7071*float(bond_length)
				y_plus = float(a[1]) + float(bond_length)
				y_plusplus = y_plus + 0.7071*float(bond_length)


				a[2] = str(b)
				d = ' '.join(a)
				oxygen_line.append(d)

				a[1] = str(y_plus)
				a[2] = str(b)
				d = ' '.join(a)
				oxygen_line.append(d)


				a[1] = str(y_minus)
				a[2] = str(c)
				d = ' '.join(a)
				hydrogen_line.append(d)

				a[1] = str(y_plusplus)
				a[2] = str(c)
				e = ' '.join(a)
				hydrogen_line.append(e)
		


	#print "\n-------------------------\n"
	#print " Absorbate atoms\n"
	
	#if object == "H" or object == "O":
		#print "Object : ", object

		#for i in absorbate_line:
			#print i

	#else:
		#print " Oxygen : "
		#for i in oxygen_line:
			#print i
		#print "\n"  
	   
		#print " Hydrogen : "
		#for i in hydrogen_line:
			#print i


	
	#print "\n---------------------------"
	#print "Total number of absorbate atoms : ", len(absorbate_line)
	#print "\n"

	#print "\n---------------------------\n\tPOSCAR_modified\n-----------------------"
	
	compact_line[coordi_type_index] = "Cartesian" 
	compact_line[atom_species_index] =  object_type + compact_line[atom_species_index]			  
   
	if object == "H" or object == "O":
	
		compact_line[num_of_atom_index] = " %d" % len(absorbate_line) + " " + compact_line[num_of_atom_index]
		POSCAR_modified = compact_line[:coordi_start] + absorbate_line + real_line
 
	else:

		compact_line[num_of_atom_index] = " %d" % len(oxygen_line) + " "+ "%d" % len(hydrogen_line) + " " + compact_line[num_of_atom_index]
		POSCAR_modified = compact_line[:coordi_start] + oxygen_line + hydrogen_line + real_line



	#for i in POSCAR_modified:
	  #print i
#
#	print "\n---------------------------\n"
	return POSCAR_modified

#MINO surface_analyzer(file,object,molecule,geometry,distance,bond_length,auto,*argv)
fin_structure = surface_analyzer(file,object,molecule,geometry,distance,bond_length,auto,*argv)
fin_structure_stringified= "\n".join(fin_structure)


print fin_structure_stringified



 
