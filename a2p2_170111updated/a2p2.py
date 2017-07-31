import glob, re, os, operator, itertools, copy, collections
from math import sin, cos, acos,radians
import parser
import our_module
import ast
import sys

# change runjob.sh
#our_module.change_vasp_path()

open('compound_directories', 'w').close()
# we choose files to work on depending on arguments given to the script, eg: l, -f
cif_files = []
poscar_files = []
arguments = sys.argv
# if it's empty so we work on all files with cif and vasp extension
if len(arguments) == 1:
    cif_files = glob.glob('*.cif') # to read all files with "cif" extension in the current directory
    poscar_files = glob.glob('*.vasp') # to read all files with "vasp" extension in the current directory

# if the argument "l" is given the we get cif and vasp file names from a list

if (len(arguments) == 3) and arguments[1] == '-l':
    f_list = open(arguments[2], 'r')
    lines = f_list.readlines()
    for line in lines:
        line_list = filter(None, line.strip().split())
        for file_name in line_list:
            if file_name[-4:] == '.cif':
                cif_files.append(file_name)
            if file_name[-5:] == '.vasp':
                poscar_files.append(file_name)

if (len(arguments)>2) and arguments[1]=='-f':
    f_list = arguments[2:]
    for file_name in f_list:
        if file_name[-4:] == '.cif':
            cif_files.append(file_name)
        if file_name[-5:] == '.vasp':
            cif_files.append(file_name)


# we first convert those cif files to postcar/vasp files :) the longest part

#######################################################################################################################

for file in cif_files:
    open(file[:-4]+'.vasp', 'w').close() # a POSCAR file is created
    poscar_file = open(file[:-4]+'.vasp', 'w')
    cif_file = open(file, 'r')                  # one of the "cif" files is opened
    line = cif_file.readline()

    # variables
    chemical_formula = ''
    chemical_formula_dict = collections.OrderedDict()
    a = 0.0
    b = 0.0
    c = 0.0
    alpha = 0.0
    beta = 0.0
    gamma = 0.0
    matrix = []
    atom_site_count = 0
    number_of_isotops = 0
    nbr_of_atoms = 0
    dict = collections.OrderedDict()
    records = collections.OrderedDict() # to build up chemical formula in case we don't have one
    formulae = []
    keys = []

    while line:

        # read the chemical formula
        if our_module.chemical_formula_extract(line):
            chemical_formula = our_module.chemical_formula_extract(line)

        # read a, b, c, alpha, beta, and gamma
        if our_module.cell_length_angle_extract('cell_length_a', line):
            a = our_module.cell_length_angle_extract('cell_length_a', line)

        if our_module.cell_length_angle_extract('cell_length_b', line):
            b = our_module.cell_length_angle_extract('cell_length_b', line)
        if our_module.cell_length_angle_extract('cell_length_c', line):
            c = our_module.cell_length_angle_extract('cell_length_c', line)

        if our_module.cell_length_angle_extract('cell_angle_alpha' , line):
            alpha = radians(our_module.cell_length_angle_extract('cell_angle_alpha' , line))

        if our_module.cell_length_angle_extract('cell_angle_beta' , line):
            beta = radians(our_module.cell_length_angle_extract('cell_angle_beta' , line))

        if our_module.cell_length_angle_extract('cell_angle_gamma' , line):
            gamma = radians(our_module.cell_length_angle_extract('cell_angle_gamma' , line))

        line  = line.strip()
        if ('_atom_site_' ==  line[:11]):

            key = line.strip().split(' ')[0]
            dict[key] = []
            atom_site_count = 1
            line = cif_file.readline()
            keys.append(key)

            line = line.strip()
            while ('_atom_site_' == line[:11]):
                key = line.strip().split(' ')[0]
                keys.append(key)
                dict[key] = []
                atom_site_count += 1
                line = cif_file.readline().strip()

            line_list = filter(None, line.strip().split(' '))
            while (len(line_list) == atom_site_count):
                for j in range(atom_site_count):
                    dict[keys[j]].append(line_list[j])
                line = cif_file.readline()
                line_list = filter(None, line.strip().split(' '))
                number_of_isotops += 1

        if '_symmetry_equiv_pos_as_xyz' in line:

            tmp = cif_file.readline().strip()

            while (len(tmp) != 0):

                if tmp[0].isdigit():
                    try:
                        tmp = tmp[tmp.index("'")+1:]
                        tmp = tmp[:tmp.index("'")]
                        tmp = [i.strip() for i in tmp.split(',')]
                        if (tmp[0] != 'x') or (tmp[1] != 'y') or (tmp[2] != 'z'):
                            formulae.append(tmp)
                        tmp = cif_file.readline().strip()
                    except ValueError:
                        tmp = tmp.split(' ')[1:]
                        tmp = ' '.join(tmp)
                        tmp = tmp.strip()
                        tmp = [i.strip() for i in tmp.split(',')]
                        if (tmp[0] != 'x') or (tmp[1] != 'y') or (tmp[2] != 'z'):
                            formulae.append(tmp)
                        tmp = cif_file.readline().strip()
                else:

                    try:
                        tmp = tmp[tmp.index("'")+1:]
                        tmp = tmp[:tmp.index("'")]
                        tmp = [i.strip() for i in tmp.split(',')]
                        if (tmp[0] != 'x') or (tmp[1] != 'y') or (tmp[2] != 'z'):
                            x_cor, y_cor, z_cor = tmp_list[0], tmp_list[1], tmp_list[2]
                            x, y, z = 1,1,1
                            code_x, code_y, code_z = parser.expr(x_cor).compile(), \
                                                     parser.expr(y_cor).compile(), parser.expr(z_cor).compile()
                            if (isinstance(eval(code_x), (int, long)) and isinstance(eval(code_y), (int, long))
                                and isinstance(eval(code_x), (int, long))):
                                formulae.append(tmp)
                            tmp = cif_file.readline().strip()

                    except ValueError:

                        if not isinstance(tmp, list):
                            tmp = [i.strip() for i in tmp.split(',')]
                        if len(tmp) == 3:
                            if (tmp[0] != 'x') or (tmp[1] != 'y') or (tmp[2] != 'z'):
                                x_cor, y_cor, z_cor = tmp[0], tmp[1], tmp[2]
                                x, y, z = 1,1,1
                                code_x, code_y, code_z = parser.expr(x_cor).compile(), \
                                                         parser.expr(y_cor).compile(), parser.expr(z_cor).compile()
                                if (isinstance(eval(code_x), (int, long)) and isinstance(eval(code_y), (int, long))
                                    and isinstance(eval(code_x), (int, long))):
                                    formulae.append(tmp)

                        elif len(tmp) != 3:
                            x = cif_file.readline()
                            cif_file.seek(-2*len(x), 1)
                            break
                        tmp = cif_file.readline().strip()


        line = cif_file.readline()

    #####################################################################
    gamma_star = 0.0
    try:
        val = (cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta))
        gamma_star = acos(val)
    except ZeroDivisionError:
        print "Division by zero!"

    matrix.append([a*sin(beta), 0.0, a*cos(beta)])
    matrix.append([-b*sin(alpha)*cos(gamma_star), b*sin(alpha)*sin(gamma_star), b*cos(alpha)])
    matrix.append([0.0, 0.0, c])



    for i in range(number_of_isotops):
        x = float(dict['_atom_site_fract_x'][i])
        y = float(dict['_atom_site_fract_y'][i])
        z = float(dict['_atom_site_fract_z'][i])

        if x < 0:
            x = 1+x
        if y < 0:
            y = 1+y
        if z < 0:
            z = 1+z

        if dict['_atom_site_type_symbol'][i] in records.keys():
            records[dict['_atom_site_type_symbol'][i]].append((x,y,z))
        else:
            records.update({dict['_atom_site_type_symbol'][i]:[]})
            records[dict['_atom_site_type_symbol'][i]].append((x,y,z))

    #################################################################

    # dealing with formulae
    ###########################################################
    for formula in formulae:
        #poscar_file.write("'"+form+"'"+"\n")

        x_cor, y_cor, z_cor = formula[0].strip(), formula[1].strip(), formula[2].strip()

        if  (x_cor != 'x' or  y_cor != 'y' or  z_cor != 'z'):

            x_cor = x_cor.replace("/","*1.0/")
            y_cor = y_cor.replace("/","*1.0/")
            z_cor = z_cor.replace("/","*1.0/")

            code_x = parser.expr(x_cor).compile()
            code_y = parser.expr(y_cor).compile()
            code_z = parser.expr(z_cor).compile()

            for i in range(number_of_isotops):
                x = float (dict['_atom_site_fract_x'][i])
                y = float (dict['_atom_site_fract_y'][i])
                z = float (dict['_atom_site_fract_z'][i])

                x_n, y_n, z_n = eval(code_x), eval(code_y), eval(code_z)

                if x_n < 0.0:
                    x_n = 1.0+x_n
                elif x_n == 0:
                    x_n = 0.0
                elif x_n > 1.0:
                    x_n = x_n%1
                if y_n < 0.0:
                    y_n = 1.0+y_n
                elif y_n == 0:
                    y_n = 0.0
                elif y_n > 1.0:
                    y_n = x_n%1
                if z_n < 0.0:
                    z_n = 1.0+z_n
                elif z_n == 0:
                    z_n = 0.0
                elif z_n > 1.0:
                    z_n = z_n%1

                if dict['_atom_site_type_symbol'][i] in records.keys():
                    records[dict['_atom_site_type_symbol'][i]].append((x_n,
                                                                       y_n,
                                                                       z_n))
                else:
                    records.update({dict['_atom_site_type_symbol'][i]:[]})
                    records[dict['_atom_site_type_symbol'][i]].append((x_n,
                                                                       y_n,
                                                                       z_n))


    for key in records.keys():
        redundance_removed = list(set(records[key]))
        records[key] = redundance_removed



    for key in records.keys():
        for cordinate1 in records[key]:
            for cordinate2 in records[key]:
                cordinate1_copy = (cordinate1[0]*a,cordinate1[1]*b,cordinate1[2]*c)
                cordinate2_copy = (cordinate2[0]*a,cordinate2[1]*b,cordinate2[2]*c)
                if (our_module.calculate_distance(cordinate1_copy, cordinate2_copy) < 1) and \
                    (cordinate1_copy != cordinate2_copy):
                    records[key].remove(cordinate2)


    for elem in records.keys():
        chemical_formula_dict.update({elem:len(records[elem])})


    if chemical_formula == '':
        for elem in records.keys():
            chemical_formula += elem + str(chemical_formula_dict[elem]) + " "

    atomicity = re.findall(r'\d+', chemical_formula)
    atom_symb = ''.join(i for i in chemical_formula if not i.isdigit())
    atom_symb = atom_symb.strip().split(' ')

    # Writing to POSCAR

    #1.
    ################################################
    poscar_file.write(chemical_formula + '\n')
    poscar_file.write('1.0\n')
    ################################################

    #2.
    for i in range(len(matrix)):
        poscar_file.write("%-9f%-9f%-9f \n" %(matrix[i][0], matrix[i][1], matrix[i][2]))

    #3.
    ###########################################################
    nbr_of_atom_symb = len(atom_symb)
    for i in range(nbr_of_atom_symb):
        if (i == nbr_of_atom_symb-1):
            poscar_file.write("%-4s\n" % atom_symb[i])
        else:
            poscar_file.write("%-4s" % atom_symb[i])
    #########################################################


    #4.
    ##########################################################
    for i in range(nbr_of_atom_symb):
        if (i == nbr_of_atom_symb-1):
            poscar_file.write("%-4d\n" % int(atomicity[i]))
        else:
            poscar_file.write("%-4d" % int(atomicity[i]))
    ##########################################################

    #5.
    ##########################################################
    poscar_file.write('direct\n')
    ##########################################################

    #6
    ############################################################
    for elem in records.keys():
        for coordinate in records[elem]:
            x, y, z = coordinate[0], coordinate[1], coordinate[2]
            poscar_file.write('%-9f%-9f%-9f%-9s\n' % (x,y,z,elem))

    poscar_file.close()
    #create folders
    our_module.create_compound_folder_cif(file[:-4]+'.vasp')

#######################################################################################################################


# Then we deal with user's poscar/vasp files
for file in poscar_files:
    our_module.create_compound_folder_poscar(file)

