#!/usr/bin/env python
import time, sys, os, copy, subprocess
from datetime import datetime
import numpy as np
from ase.io import vasp, read, write
import NanoCore as nc

version = 20171211
now_time = datetime.now()
usage    = ' Usage: %s [POSCAR or CONTCAR] \n ' % sys.argv[0]
foottext = '\n Thank you\n## Yoon Su Shim <dbstn145@gmail.com>'
print "## Select fixed atom"
print "## Version : %s " % version
print now_time


# INPUT TAGS

atom_label = True # True or False
fix_lines = range(37, 60+1)
#fix_lines = [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]
print fix_lines


def read_POSCAR(filename):

        # READ POSCAR

	f = open(filename)

	lines = f.readlines()

	nspecies = lines[6].split()

        natms = 0     
        for na in nspecies :
                natms += int(na)
       
        if lines[7][0] == 'S' or lines[7][0] == 's' :
                start = 9
        else :
                start = 8

	position_lines = lines[start:start+natms]

	position = []
	for line in position_lines:

		line = line.split()[0:3]
		line2 =[]

		for li in line:

			line2.append(float(li))

		position.append(line2)

        # END READ POSCAR

	return lines, position


def write_POSCAR(filename, lines, position, fix_lines):
        
        # WRITE POSCAR

        f = open(filename, 'w')
        f.write(''.join(lines[0:7]))
        f.write('Selective Dynamics\nDirect\n')

        for i in range(len(position)) :

                f.write(' %19.16f %19.16f %19.16f' %(position[i][0], position[i][1], position[i][2]))

                for j in fix_lines :

                        if i == j :
                                fix = True
                                break 
                        else:
                                fix = False

                if fix == True:

                        f.write('     F   F   F\n') 

                else : 

                        f.write('     T   T   T\n')                                                       

        f.close()

        # END WRITE POSCAR


def write_atom_label(filename):
        ### atom_label.py
        ### 2017-8-9
        ### Developer: Kanghoon Yim (feihoom82@gmail.com)

        # read POSCAR
        infile = open(filename, 'rU')
        poscar = infile.readlines()
        infile.close()

        outfile = open(filename+'_atom_label', 'w')

        mag = float(poscar[1])
        atom = poscar[5].split()
        cnt = []
        sum = 0
        for i in range(len(atom)) :
        	cnt[i:] = [int(poscar[6].split()[i])]
        	sum = sum + cnt[i]
        if poscar[7].split()[0][0] == 'S' or poscar[7].split()[0][0] == 's' :
        	start = 9
        else :
	        start = 8

        for i in range(start) :
	        outfile.write(poscar[i])

        line = start
        for i in range(len(cnt)) :
	        for j in range(cnt[i]) :
		        outfile.write(poscar[line].replace('\n','')+'\t!'+atom[i]+str(j+1)+'\n')
		        line=line+1
        outfile.close()


# CHECK INPUT FILE

if len(sys.argv) == 2:
        filename = sys.argv[1]

else :
        print usage
        print foottext
        sys.exit(1)

lines, position = read_POSCAR(filename)
filename_out = filename+'_fix'
write_POSCAR(filename_out, lines, position, fix_lines)

if atom_label == True:
        write_atom_label(filename_out)
else: 
        sys.exit(1)
