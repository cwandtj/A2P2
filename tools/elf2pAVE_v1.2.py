#!/usr/bin/env python
import time, sys, os, subprocess
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

version = 20171208 # 20141111 # 20130704
now_time = datetime.now()
usage    = ' Usage: %s [ELFCAR or CHGCAR] (default: ELFCAR) \n ' % sys.argv[0]
foottext = '\n Thank you\n## Yoon Su Shim (CERU, KIER) <dbstn145@gmail.com>'
print "## Plane average for grid (CHGCAR format)"
print "## Version : %s \n" % version
print now_time


# INPUT TAGS

mincut =  -0.00001
maxcut = 150.00000 
contour_levels = np.linspace(0.000, 0.5, 201) # start, stop, number of points
figure_show = False # False or True
#slice_range = [0,24] # grid range
average_axis = 'y' # If you this tag set to x, y, or z, show z-y plane, z-x plane, or y-x plane, respectively


def read_ELFCAR( filename, ispin ):

        # READ ELFCAR

        print 'Reading ... ELFCAR'

	f = open(filename)

	lines = f.readlines()
	nlines = len(lines)

        # GET ATOMIC DATA

	lattice_lines = lines[2:5]
	species_lines = lines[5].split()
	nspecies = [int(nsp) for nsp  in lines[6].split()]

        natm = 0
        for nsp in nspecies:
                natm += nsp

	position_lines = lines[9:natm+8]

	lattice = []; position = []; species = []
        species_list = {}

	i= 0
	while i< len(species_lines):
		for ns in range(int(nspecies[i])):
			species.append(species_lines[i])
                species_list['%s' % species_lines[i]] = nspecies[i]
		i+=1

        print species_list

	for line in lattice_lines:
		line = line.split()[0:3]
		line2 = []
		for li in line:
			line2.append(float(li))
		lattice.append(line2)

	for line in position_lines:
		line = line.split()[0:3]
		line2 =[]
		for li in line:
			line2.append(float(li))
		position.append(line2)


        # GET ELF GRID DATA
        
        ii = 0; numlist =[]
        while ii < len(lines):
                line = lines[ii].split()
                if len(line) == 3:
                        ini2 = ii  
                        print line
                        numlist.append(ii)
                        print ini2
                ii+=1
        print numlist

        if ispin == 1:
                elf_ini = numlist[-1]
        elif ispin == 2: 
                elf_ini = numlist[-2]
                elf_ini2 = numlist[-1]

        #elf_grid = [ int(grid) for grid in lines[natm+9].split() ]
        #elf_lines = lines[ natm+10 : natm+11+nelf_grid ]
        elf_grid  = [ int(grid) for grid in lines[elf_ini].split()]
        nelf_grid = elf_grid[0] * elf_grid[1] * elf_grid[2] / 10
        elf_lines = lines[ elf_ini+1 : elf_ini2 ]
        print elf_lines[0]
        print elf_lines[-1]

        elf_points = np.zeros( ( 2, elf_grid[0] * elf_grid[1] * elf_grid[2] ) )

        x=0; temp =[]
        while x < len( elf_lines ):
                line = elf_lines[x].split()
                x += 1
                for y in line:
                        temp.append(float(y))
        elf_points[0] = temp
        elf_points = np.array(elf_points)

        if ispin == 2:
                
                #elf_grid2  = [int(grid) for grid in lines[ natm+11+nelf_grid ].split()] 
                #elf_lines2 = lines[ natm+12+nelf_grid : natm+13+2*nelf_grid ]
                elf_grid2 = [ int(grid) for grid in lines[elf_ini2].split() ]
                elf_lines2 = lines[ elf_ini2+1: nlines]

                print elf_lines2[0]
                print elf_lines2[-1]

                
                x=0; temp=[]
                while x < len(elf_lines2):
                        line = elf_lines2[x].split()
                        x += 1
                        for y in line:
                                temp.append(float(y))
                elf_points[1] = temp
        else:
                pass

        # END READ ELFCAR

	return elf_grid, elf_points, lattice, position


def read_ISPIN():

        # READ ISPIN TAG

        hosts = subprocess.Popen(['grep','ISPIN','OUTCAR'],stdout=subprocess.PIPE)
        hosts_out = hosts.stdout.read()        
        ispin = int(hosts_out.split()[2])

        print 'ISPIN = ', ispin, '\n'

        # END READ ISPIN TAG

        return ispin


def set_min_max( mincut, maxcut, elf_points ):

        i=0
        while i < len(elf_points):

	        if elf_points[i] < mincut:

		        elf_points[i] = mincut

        	elif elf_points[i] > maxcut:

	        	elf_points[i] = maxcut

	        i+=1

        elf_points = np.array(elf_points)

        return elf_points


#def average_matrix( elf_grid, elf_points, average_axis, slice_range ):
def average_matrix( elf_grid, elf_points, average_axis ):

         # AVERAGE MATRIX       

        print 'Averaging ...'

        values_appended=[]


        if average_axis == 'y':
                cc=0
                while cc < int(elf_grid[2]):

	                aa=0
        	        target1=[]
	                while aa < int(elf_grid[0]):

#		                bb = slice_range[0]
                                bb = 0

                                target = float(elf_points[cc][0][aa])
        		        while bb < int(elf_grid[1]):
#                                while bb < slice_range[1]:
	        		        target = target + float(elf_points[cc][bb][aa])
		        	        bb=bb+1
#                                target = target/(slice_range[1]-slice_range[0])
                                target = target/bb
        		        target1.append(target)
	        	        aa=aa+1

	                values_appended.append(target1)
        	        cc=cc+1


        elif average_axis =='x':

                cc=0
                while cc < int(elf_grid[2]):

	                bb=0
	                target1=[]
        	        while bb < int(elf_grid[1]):

	        	        aa = 0

                                target = float(elf_points[cc][bb][0])
		                while aa < int(elf_grid[1]):

			                target = target + float(elf_points[cc][bb][aa])
        			        aa=aa+1

                                target = target/aa
		                target1.append(target)
		                bb=bb+1

        	        values_appended.append(target1)
        	        cc=cc+1


        elif average_axis =='z':

                bb=0
                while bb < int(elf_grid[1]):

	                aa=0
	                target1=[]
        	        while aa < int(elf_grid[0]):

	        	        cc = 0

                                target = float(elf_points[0][bb][aa])
		                while cc < int(elf_grid[2]):

			                target = target + float(elf_points[cc][bb][aa])
        			        cc=cc+1

                                target = target/cc
		                target1.append(target)
		                aa=aa+1

        	        values_appended.append(target1)
        	        bb=bb+1


        return values_appended


def plot_contour( lattice, elf_grid, values_appended, figure_show, levels, ispin, spin ):

        # GET GRID CELL

        x_step = np.linspace( 0, float(lattice[0][0]), float(elf_grid[0]) )
        y_step = np.linspace( 0, float(lattice[1][1]), float(elf_grid[1]) )
        z_step = np.linspace( 0, float(lattice[2][2]), float(elf_grid[2]) )

        # PLOT CONTOUR

        if average_axis == 'y':
                X, Z = np.meshgrid( x_step, z_step )
                cset = plt.contourf(X, Z, values_appended,levels)
        elif average_axis == 'x':
                Y, Z = np.meshgrid( y_step, z_step )
                cset = plt.contourf(Y, Z, values_appended,levels)
        elif average_axis == 'z':
                X, Y = np.meshgrid( x_step, y_step )
                cset = plt.contourf(X, Y, values_appended,levels)

        plt.colorbar(cset)#, ticks=[])

        if figure_show == False:
                pass
        else:
                print 'Plot contour...'
                plt.show()


        # SAVE RESULTS

        print 'Save data...'

        if ispin == 1:
                
                plt.savefig('ELFCAR_pAVE.png')
        else :
                if spin == 0:
                        plt.savefig('ELFCAR_pAVE_up.png')
                else:
                        plt.savefig('ELFCAR_pAVE_down.png')


# START PROGRAM

# CHECK INPUT FILE

if len(sys.argv) ==1:
        filename = 'ELFCAR'

elif len(sys.argv) == 2:
        filename = sys.argv[1]

else:
        print usage
        print foottext
        sys.exit(1)

if average_axis == 'x' or average_axis == 'X' or\
        average_axis == 'y' or average_axis == 'Y' or\
        average_axis == 'z' or average_axis == 'Z':

        print 'Averaged along "',average_axis,'"\n'

else:
        print 'ERROR: Invalid axis ! "', average_axis,'"'
        sys.exit(1)
                                        
# READ DATA

ispin =  read_ISPIN()                                                                  
elf_grid, elf_points, lattice, position = read_ELFCAR( filename, ispin )

# RUN PROGRAM

for spin in range(ispin):
        print '\nspin :', spin
        elf = elf_points[spin]
        elf = set_min_max( mincut, maxcut, elf )
        print max(elf), min(elf)
        elf = elf.reshape( ( elf_grid[2], elf_grid[1], elf_grid[0] ) )
#        values_appended = average_matrix( elf_grid, elf, average_axis, slice_range )
        values_appended = average_matrix( elf_grid, elf, average_axis )
        plot_contour( lattice, elf_grid, values_appended, figure_show, contour_levels, ispin, spin )

# END PROGRAM
