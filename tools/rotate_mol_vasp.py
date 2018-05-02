#################################################
## rotate_mol_vasp.py							#
## Decription: Rotate a molecule from POSCAR	#
## Developver: feihoom82@gmail.com				#
## Date: 2018-05-02								#
#################################################

### Start of Seting ###
center = 65	# Index of the center atom
ids_mol = [65, 66, 67, 68, 69]	# The center atom must be the first index.
rotate1 = 42	# Angle(degree) or initial direction vector
rotate2 = [0,0,1]	# Rotational axis  vector
# Example1) rotate1 = 45; rotate2 = [0,0,1]
# Example2) rotate1 = [1,0,0]; rotate2 = [1,1,0]
### End of Setting ###

from ase import Atoms
from ase.build import molecule
import sys

# Set first index as '0'
center = center-1
for i in range(len(ids_mol)) :
	ids_mol[i] = ids_mol[i]-1

### Read POSCAR
if len(sys.argv) < 2 :
	print 'Usage: python rotate_mol_vasp.py POSCAR'
	print 'Notice - You should set parameters in the source code !'
	sys.exit(-1)
POSCAR = open(sys.argv[1],'r').readlines()

### Read structure and elem positions from POSCAR
elem = POSCAR[5].split()
nelem = []
for i in range(len(elem)) :
	nelem[i:] = [int(POSCAR[6].split()[i])]
n_tot = sum(nelem)
if POSCAR[7].split()[0][0] in ['S','s'] :
	start = 9
else :
	start = 8
if POSCAR[start-1].split()[0][0] in ['C','c'] :
	coord_type = 'Cartesian'
elif POSCAR[start-1].split()[0][0] in ['D','d'] :
	coord_type = 'Direct'
mag = float(POSCAR[1].strip())
a = [float(x)*mag for x in POSCAR[2].split()]
b = [float(x)*mag for x in POSCAR[3].split()]
c = [float(x)*mag for x in POSCAR[4].split()]
xyz = POSCAR[start:start+n_tot]

cnt = 0
for i in range(len(elem)) :
	for j in range(nelem[i]) :
		if coord_type == 'Cartesian' :
			xyz[cnt] = [float(x) for x in xyz[cnt].split()[0:3], elem[i]]
		elif coord_type == 'Direct' :
			vec = [float(x) for x in xyz[cnt].split()[0:3]]
			xyz[cnt] = [a[0]*vec[0]+b[0]*vec[1]+c[0]*vec[2], a[1]*vec[0]+b[1]*vec[1]+c[1]*vec[2], a[2]*vec[0]+b[2]*vec[1]+c[2]*vec[2], elem[i]]
		cnt += 1

### Generate a molecule from POSCAR
chem_mol = []
pos_mol = []
for idx in ids_mol :
	chem_mol.append(xyz[idx][3])
	pos_mol.append([xyz[idx][0]-xyz[center][0],xyz[idx][1]-xyz[center][1],xyz[idx][2]-xyz[center][2]])

MOL = Atoms(''.join(chem_mol),pos_mol)
#ooh  = Atoms("OOH", [[0, 0, 0], [-1.067, -0.403, 0.796], [-0.696, -0.272, 1.706]])
#h2o  = molecule('H2O')
print 'Molecule : '+MOL.get_chemical_formula()
print '\nBefore rotation:'
print MOL.get_positions()

MOL.rotate(rotate1,rotate2,center=(0,0,0))
print '\nAfter rotation:'
print MOL.get_positions()

### Change the posision of atom in the molecule
for i in range(len(ids_mol)) :
	xyz[ids_mol[i]][0] = MOL.get_positions()[i][0]+xyz[center][0]
	xyz[ids_mol[i]][1] = MOL.get_positions()[i][1]+xyz[center][1]
	xyz[ids_mol[i]][2] = MOL.get_positions()[i][2]+xyz[center][2]

### Write POSCAR_out
out = open(sys.argv[1]+'_out','w')
out.write(''.join(POSCAR[0:start-1]))
out.write('Cartesian'+'\n')
SelDyn = ''
if start == 9 :
	SelDyn = '\tT   T   T'
for i in range(len(xyz)) :
	out.write('  %13.9F  %13.9F  %13.9F' %(xyz[i][0],xyz[i][1],xyz[i][2])+SelDyn+'\n')
out.close()

print '\nWritting '+sys.argv[1]+'_out file is done.'
