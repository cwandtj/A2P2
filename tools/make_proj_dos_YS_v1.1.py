#!/usr/bin/env python

# Takes a VASP PROCAR file and produces projected densities of states
# by binning the eigenstates on to a mesh ('raw')
# and then Gaussian smearing ('smear').

# Latest updated 21/11/2017 by Yoon-Su Shim, CERU, Korea Institute of Energy Research
# Written by M. J. Wolf, Department of Chemistry, Uppsala University
# Completed 01/11/2016

import numpy as np
import copy,sys,subprocess
from datetime import datetime

# INPUT TAGS
normalize = True  # normailized by number of each atoms
print_raw = False # Prints unsmeared DOS
print_all = False # If this is set to true then it prints projection on to each and every atom individually
sigma_r = 0.1 # dispersion for Gaussian smearing
ngridpts = 1500 # number actually equals 2*ngridpts+1
min_en = None  # min energy
max_en = None  # max energy


def read_POSCAR(filename):

        # READ POSCAR

	f = open(filename)

	lines = f.readlines()
	nlines = len(lines)
	lattice_lines = lines[2:5]
	species_lines = lines[5].split()
	nspecies = lines[6].split()
	position_lines = lines[9:nlines]

	lattice = []; position = []; species = []
        species_list ={}
	
	i= 0
	while i< len(species_lines):
		for ns in range(int(nspecies[i])):
			species.append(species_lines[i])
                species_list['%s' % species_lines[i]] = nspecies[i]
		i+=1

	for line in lattice_lines:
		line = line.split()
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

        # END READ POSCAR

	return species, nspecies, lattice, position, species_lines, species_list


def read_DOSCAR():

        # READ DOSCAR

        filename = 'DOSCAR'

        with open( filename, 'r') as file:

            lines = file.readlines()
            E_F = float(lines[5].split()[3])

        print 'Fermi energy: ', E_F

        # END READ DOSCAR

        return E_F


def read_ISPIN():

        # READ ISPIN TAG

        hosts = subprocess.Popen(['grep','ISPIN','OUTCAR'],stdout=subprocess.PIPE)
        hosts_out = hosts.stdout.read()        
        ispin = int(hosts_out.split()[2])

        print 'ISPIN = ', ispin

        # END READ ISPIN TAG

        return ispin


def read_PROCAR(ispin) :

    # READ PROCAR

    filename = 'PROCAR'

    with open( filename, 'r') as file :
    
        lines = file.readlines()

        line = lines[1].split()
        NKPTS = int(line[3])
        NBANDS = int(line[7])
        NIONS = int(line[11])
        idx = lines[7].split()
        num_dat_col = len(idx)-1
        print 'Number of data column : ',num_dat_col

        eigenvals = np.zeros( [ NKPTS, NBANDS, 2 ] )
        kptweights = np.zeros( NKPTS )
        pevals = np.zeros( [NKPTS, NBANDS, 2, NIONS, num_dat_col+2] )

        spin1 = lines[2:NKPTS*(3+NBANDS*(5+NIONS))+2]

        for k in range( NKPTS ) :

                kpt = spin1[1+(3+NBANDS*(5+NIONS))*k]
                kptweights[k] = float(kpt.split()[-1])

                kpts = spin1[(NBANDS*(5+NIONS)+3)*k:(NBANDS*(5+NIONS)+3)*(k+1)]

                for b in range(NBANDS):
                        eigv = float(kpts[(5+NIONS)*b+3].split()[4])
                        eigenvals[k,b,0]=eigv
                        eigenvals[k,b,1]=eigv

                        bands = kpts[3+(5+NIONS)*b:3+(5+NIONS)*(b+1)]
                        j=0
                        for i in range(3, NIONS+3 ) :
                                line = bands[i]
                                tmp1 = [np.float( l ) for l in line.split()[1:] ]
                                tmp1.append(np.sum(tmp1[1:4]))
                                tmp1.append(np.sum(tmp1[4:9]))
                                tmp1 = np.array(tmp1)

                                pevals[k,b,0,j,:]=tmp1
                                pevals[k,b,1,j,:]=-tmp1
                                j+=1

        if ispin == 2:
                spin2= lines[NKPTS*(3+NBANDS*(5+NIONS))+3:(NKPTS*(3+NBANDS*(5+NIONS))+2)*2]

                for k in range( NKPTS ) :

                        kpts = spin2[(NBANDS*(5+NIONS)+3)*k:(NBANDS*(5+NIONS)+3)*(k+1)]

                        for b in range(NBANDS):
                        
                                eigv = float(kpts[(5+NIONS)*b+3].split()[4])
                                eigenvals[k,b,1]=eigv

                                bands = kpts[3+(5+NIONS)*b:3+(5+NIONS)*(b+1)]
                                j=0

                                for i in range(3, NIONS+3 ) :

                                        line = bands[i]
                                        tmp1 = [np.float( l ) for l in line.split()[1:] ]
                                        tmp1.append(np.sum(tmp1[1:4]))
                                        tmp1.append(np.sum(tmp1[4:9]))
                                        tmp1 = np.array(tmp1)

                                        pevals[k,b,1,j,:]=-tmp1
                                        j+=1

        idxs = {};idd = 0
        for ii in idx:
                idxs['%i' % idd] = ii
                idd+=1
        idxs['11']='p'; idxs['12']='d'
        print idxs

    # END READ PROCAR

    return [eigenvals, kptweights, pevals, idxs, num_dat_col]


def bin(eigenvals, kptweights, bandweights, ngridpts, min_en = None, max_en = None) :

    # PUT SPECTRUM ON 1D GRID

    centpt = ngridpts
    ngridpts = 2 * ngridpts + 1

    if min_en == None :

        min_en =  eigenvals.min() - 5.0

    if max_en == None :

        max_en =  eigenvals.max() + 5.0

    dE = (max_en - min_en)/ngridpts
    grid = np.multiply( np.arange( ngridpts ), dE)
    grid = np.add( grid, min_en )
    grid_data = np.zeros( [ ngridpts, 2 ] )
    NKPTS = len( kptweights )
    NBANDS = len( bandweights[0,:,0] )

    for k in range( NKPTS ) :
    
        kweight = kptweights[k]
    
        for b in range( NBANDS ) :

           gridpt = int( ( eigenvals[k,b,0] - min_en ) // dE )

           if 0 <= gridpt < ngridpts :

               bweight = bandweights[k,b,0]
               grid_data[gridpt,0] = grid_data[gridpt,0] + kweight*bweight
   
           gridpt = int( ( eigenvals[k,b,1] - min_en ) // dE )

           if 0 <= gridpt < ngridpts :

               bweight = bandweights[k,b,1]
               grid_data[gridpt,1] = grid_data[gridpt,1] + kweight*bweight

    grid_data = grid_data / dE

    # END PUT SPECTURM ON 1D GRID

    return [grid, grid_data]


def smear( grid, grid_data, sigma_r ) :

    # SMEAR DATA ON GRID

    max_en = grid.max()
    min_en = grid.min()
    dt = 1.0/( max_en - min_en )
    ngridpts = len( grid )
    centpt = ( ngridpts - 1 )/2
    centpt = centpt*dt
    
    recip_grid = np.arange( ngridpts )
    recip_grid = np.multiply( recip_grid, dt )
    
    a = 2*np.pi**2*sigma_r**2
    arg = np.multiply( -np.multiply( recip_grid - centpt, recip_grid - centpt ), a )
    gauss = np.exp( arg )
    gauss = np.fft.ifftshift( gauss )
    
    new_grid_data = copy.deepcopy( grid_data )
    
    tmp_grid_data = grid_data[:,0]
    recip_grid_data = np.fft.fft( tmp_grid_data )
    tmp_grid_data = np.fft.ifft( np.multiply( recip_grid_data, gauss ) )
    new_grid_data[:,0] = np.multiply( np.sign( np.real( tmp_grid_data ) ) , np.abs( tmp_grid_data ) )
    
    
    tmp_grid_data = grid_data[:,1]
    recip_grid_data = np.fft.fft( tmp_grid_data )
    tmp_grid_data = np.fft.ifft( np.multiply( recip_grid_data, gauss ) )
    
    new_grid_data[:,1] = np.multiply( np.sign( np.real( tmp_grid_data ) ) , np.abs( tmp_grid_data ) )

    # END SMEAR DATA ON GRID
    
    return new_grid_data


def RunAndPrint(E_F, eigenvals, kptweights, projwfn, NKPTS, NBANDS, NIONS, sel_dat_col, print_raw, ions,\
        sigma_r, ngridpts, min_en, max_en, orbital_name, normalize):  

        # RUN AND PRINT

        for i in range( len( ions ) ) :

            print 'Group', i
            group = ions[i]
            ngroup =len(group)
            print 'Number of atoms in group: ',ngroup

            for j in range( ngroup ) :

                data_on_grid = bin(eigenvals, kptweights, projwfn[:,:,:,group[j]-1,sel_dat_col-1], ngridpts, min_en, max_en)
    
                grid = data_on_grid[0]

                if j == 0 :

                    grid_data = data_on_grid[1]

                else :

                    grid_data = grid_data + data_on_grid[1]

            smear_grid_data = smear( grid, grid_data, sigma_r )

            if print_raw == True :

                with open( 'pdos_' + str(species_lines[i])+'_'+ str(sel_dat_col)+'_'+str(orbital_name) + '_raw.dat', 'w' ) as file :   
        
                    for j in range( len(grid) ) :

                            if normalize == True:

                                file.write('%15.10f %15.10f %15.10f %15.10f\n' %(grid[j], grid[j]-E_F, grid_data[j,0]/ngroup, grid_data[j,1]/ngroup))       

                            else:

                                file.write('%15.10f %15.10f %15.10f %15.10f\n' %(grid[j], grid[j]-E_F, grid_data[j,0], grid_data[j,1]))       

            with open( 'pdos_' + str(species_lines[i])+'_'+ str(sel_dat_col)+'_'+str(orbital_name) + '_smear.dat', 'w' ) as file :   
    
                for j in range( len(grid) ) :

                        if normalize == True:
    
                            file.write('%15.10f %15.10f %15.10f %15.10f\n' %(grid[j], grid[j]-E_F, smear_grid_data[j,0]/ngroup, smear_grid_data[j,1]/ngroup))

                        else:

                            file.write('%15.10f %15.10f %15.10f %15.10f\n' %(grid[j], grid[j]-E_F, smear_grid_data[j,0], smear_grid_data[j,1]))

        # END RUN AND PRINT


# START PROGRAM
start_time = datetime.now()

# CHECK INPUT 
if len(sys.argv) == 2:

	f = open(sys.argv[1])
	lines = f.readlines()

        ions = []; species_lines=[]; j=0

        for line in lines:

                group = line.split('\n')[0]

                if ',' in group:

                        group = group.split(',')

                elif '-' in group:

                        group = group.split('-')
                        group = range(int(group[0])-1,int(group[1]))

                print group

                i=0
                print '# element in group',len(group)
                while i<len(group):

                        group[i]=int(group[i])
                        i+=1    

                ions.append(group)
                j+=1
                species_lines.append('group%i' % j)

else:

        species, nspecies, lattice, positions,species_lines, species_list = read_POSCAR('CONTCAR')

        ions = []; temp1 = 0; temp2 = 0

        for element in species_lines:

                nelement = int( species_list[element] )
                temp2 += nelement
                ions.append( range( temp1,temp2 ) )
                print range( temp1,temp2 )
                temp1 = temp2

        if not len( species_lines ) == 1:

                ions.append( range( temp2 ) )
                species_lines.append('All')

print species_lines

# READ DATA
E_F   = read_DOSCAR()
ispin = read_ISPIN()
data  = read_PROCAR( ispin )

eigenvals   = data[0] # indices are (k-point, band, spin channel)
kptweights  = data[1] # just k-point weights
projwfn     = data[2] # indices are k-point, band, spin channel, ion, data column
idxs        = data[3] # indices are orbital column
num_dat_col = data[4]

NKPTS  = len( kptweights )
NBANDS = len( projwfn[0,:,0,0,0] )
NIONS  = len( projwfn[0,0,0,:,0] )

if print_all == True :

    ions = [[]]*NIONS

    for i in range(NIONS) :

        ions[i]= [i+1]

# RUN PROGRAM
for i in range( num_dat_col+2 ):

        iss = str(i+1)
        orbital_name = idxs[iss]
        sel_dat_col = i+1
        print 'orbial col.', sel_dat_col, orbital_name
        RunAndPrint( E_F, eigenvals, kptweights, projwfn, NKPTS, NBANDS, NIONS, sel_dat_col, print_raw, ions, \
                     sigma_r, ngridpts, min_en, max_en, orbital_name, normalize )

# END PROGAM
end_time = datetime.now()
elapsed_time = end_time - start_time

print 'End time:    ', end_time
print 'Elapsed time:', elapsed_time
