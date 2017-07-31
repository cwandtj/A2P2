from __future__ import division
from fractions import gcd
from sortedcontainers import SortedDict
import glob, re, operator, itertools, parser, ast, sys,os, collections, subprocess
from math import sin, cos, acos,radians,sqrt
import our_module, error_handling
import threading, time

all_elements = {'Ac':['Ac'], 'Ag':['Ag', 'Ag_GW','Ag_pv'], 'Al':['Al','Al_GW','Al_sv_GW'],
                'Am':['Am'], 'Ar':['Ar','Ar_GW'], 'As':['As', 'As_d', 'As_GW'], 'At':['At', 'At_d'],
                'Au':['Au', 'Au_GW', 'Au_pv_GW'], 'B':['B','B_GW','B_h', 'B_s'], 'Ba':['Ba_sv','Ba_sv_GW'], 'Be':['Be','Be_GW','Be_sv','Be_sv_GW'],
                'Bi':['Bi', 'Bi_d','Bi_d_GW','Bi_GW'], 'Br':['Br','Br_GW'], 'C':['C','C_GW','C_GW_new','C_h','C_s'], 'Ca':['Ca_pv','Ca_sv','Ca_sv_GW',],
                'Cd':['Cd', 'Cd_GW','Cd_pv_GW','Cd_sv_GW'], 'Ce':['Ce', 'Ce_3','Ce_GW','Ce_h'], 'Cl':['Cl', 'Cl_GW','Cl_h'],
                'Cm':['Cm'], 'Co':['Co', 'Co_GW', 'Co_pv','Co_sv','Co_sv_GW'], 'Cr':['Cr', 'Cr_pv', 'Cr_sv','Cr_sv_GW'],
                'Cs':['Cs_sv','Cs_sv_GW',], 'Cu':['Cu', 'Cu_GW','Cu_pv','Cu_pv_GW'], 'Dy':['Dy', 'Dy_3'],
                'Er':['Er','Er_2','Er_3'], 'Eu':['Eu', 'Eu_2', 'Eu_3'], 'F':['F','F_GW', 'F_GW_new','F_h','F_s'], 'Fe':['Fe','Fe_GW','Fe_pv','Fe_sv','Fe_sv_GW'],
                'Fr':['Fr_sv'],'Ga':['Ga', 'Ga_d', 'Ga_d_GW', 'Ga_GW', 'Ga_h', 'Ga_pv_GW','Ga_sv_GW'], 'Gd':['Gd', 'Gd_3'],
                'Ge':['Ge', 'Ge_d','Ge_d_GW','Ge_GW', 'Ge_h','Ge_sv_GW'], 'H':['H', 'H1.25', 'H1.33','H1.5','H1.66','H1.75','H.25','H.33','H.42','H.5','H.58','H.66','H.75','H_AE','H_GW','H_h','H_h_GW','H_s'],
                'He':['He', 'He_GW'], 'Hf':['Hf','Hf_pv','Hf_sv','Hf_sv_GW'], 'Hg':['Hg'], 'Ho':['Ho', 'Ho_3'], 'I':['I','I_GW'],
                'In':['In','In_d','In_d_GW'], 'Ir':['Ir','Ir_sv_GW'], 'K':['K_pv','K_sv','K_sv_GW'], 'Kr':['Kr','Kr_GW'],'La':['La', 'La_s'],
                'Li':['Li', 'Li_AE_GW','Li_GW','Li_sv','Li_sv_GW'],'Lu':['Lu','Lu_3'], 'Mg':['Mg','Mg_GW','Mg_pv','Mg_pv_GW','Mg_sv','Mg_sv_GW'],
                'Mn':['Mn','Mn_GW','Mn_pv','Mn_sv','Mn_sv_GW'], 'Mo':['Mo','Mo_pv','Mo_sv','Mo_sv_GW'],'N':['N','N_GW','N_GW_new','N_h','N_s','N_s_GW'],
                'Na':['Na', 'Na_pv','Na_sv','Na_sv_GW'], 'Nb':['Nb_pv','Nb_sv','Nb_sv_GW'], 'Nd':['Nd','Nd_3'], 'Ne':['Ne', 'Ne_GW','Ne_GW_soft'],
                'Ni':['Ni', 'Ni_GW', 'Ni_pv','Ni_sv_GW'], 'Np':['Np','Np_s'], 'O':['O','O_GW','O_GW_new','O_h','O_s','O_s_GW'],
                'Os':['Os', 'Os_pv', 'Os_sv_GW'], 'P':['P','P_GW','P_h'], 'Pa':['Pa','Pa_s'], 'Pb':['Pb', 'Pb_d','Pb_d_GW'], 'Pd':['Pd','Pd_GW','Pd_pv'],
                'Pm':['Pm', 'Pm_3'], 'Po':['Po','Po_d'], 'Pr':['Pr','Pr_3'], 'Pt':['Pt','Pt_GW','Pt_pv','Pt_pv_GW','Pt_sv_GW'],
                'Pu':['Pu','Pu_s'], 'Ra':['Ra_sv'],'Rb':['Rb_pv','Rb_sv','Rb_sv_GW'], 'Re':['Re','Re_pv','Re_sv_GW'],
                'Rh':['Rh','Rh_GW','Rh_pv','Rh_pv_GW','Rh_sv_GW'], 'Rn':['Rn'], 'Ru':['Ru', 'Ru_pv','Ru_pv_GW','Ru_sv','Ru_sv_GW'],
                'S':['S','S_GW','S_h'],'Sb':['Sb', 'Sb_d_GW','Sb_GW'], 'Sc':['Sc', 'Sc_sv', 'Sc_sv_GW'], 'Se':['Se', 'Se_GW'],
                'Si':['Si', 'Si_GW','Si_sv_GW'], 'Sm':['Sm','Sm_3'], 'Sn':['Sn', 'Sn_d','Sn_d_GW'],'Sr':['Sr_sv','Sr_sv_GW'],
                'Ta':['Ta', 'Ta_pv','Ta_sv_GW'], 'Tb':['Tb', 'Tb_3'], 'Tc':['Tc', 'Tc_pv', 'Tc_sv', 'Tc_sv_GW'], 'Te':['Te', 'Te_GW'],
                'Th':['Th', 'Th_s'], 'Ti':['Ti', 'Ti_pv','Ti_sv','Ti_sv_GW'],'Tl':['Tl', 'Tl_d'], 'Tm':['Tm', 'Tm_3'], 'U':['U', 'U_s'],
                'V':['V', 'V_pv', 'V_sv', 'V_sv_GW'],'W':['W', 'W_pv','W_sv_GW'], 'Xe':['Xe', 'Xe_GW'], 'Y': ['Y_sv','Y_sv_GW'],'Yb':['Yb', 'Yb_2'],
                'Zn':['Zn', 'Zn_GW', 'Zn_pv_GW','Zn_sv_GW'], 'Zr':['Zr_sv','Zr_sv_GW']
                }

preferred_psp = {'Pr': 'Pr_3', 'Ni': 'Ni_pv', 'Yb': 'Yb', 'Pd': 'Pd_GW', 'Pt': 'Pt_GW',
        'Ru': 'Ru_pv', 'Na': 'Na_pv', 'Nb': 'Nb_pv', 'Nd': 'Nd_3', 'C': 'C_GW',
        'Li': 'Li_sv_GW', 'Pb': 'Pb_d_GW', 'Y': 'Y_sv_GW', 'Tl': 'Tl_d_GW', 'Tm': 'Tm_3',
        'Rb': 'Rb_sv_GW', 'Ti': 'Ti_pv', 'Rh': 'Rh_pv', 'Tc': 'Tc_pv', 'Ta': 'Ta_pv',
        'Be': 'Be_sv_GW', 'Sm': 'Sm_3', 'Ba': 'Ba_sv_GW', 'Bi': 'Bi_d_GW', 'La': 'La_GW',
        'As': 'As_GW', 'Po': 'Po','Fe': 'Fe_pv', 'Br': 'Br_GW', 'Dy': 'Dy_3', 'Pm': 'Pm_3',
        'Hf': 'Hf_pv', 'K': 'K_sv_GW', 'At': 'At_d', 'Tb': 'Tb_3', 'Mg': 'Mg_pv_GW',
        'B': 'B_GW', 'F': 'F_GW', 'Sr': 'Sr_sv_GW','Mo': 'Mo_pv', 'Mn': 'Mn_pv', 'Lu': 'Lu_3',
        'O': 'O_GW', 'N': 'N_GW', 'Eu': 'Eu', 'Sn': 'Sn_d_GW', 'W': 'W_pv', 'V': 'V_sv_GW',
        'Sc': 'Sc_sv_GW', 'Os': 'Os_pv', 'Se': 'Se_GW','Hg': 'Hg', 'Zn': 'Zn_GW', 'Co': 'Co_GW',
        'Ag': 'Ag_GW', 'Re': 'Re_pv', 'Ca': 'Ca_sv_GW', 'Ir': 'Ir', 'Al': 'Al_GW', 'Ce': 'Ce_GW',
        'Cd': 'Cd_GW', 'Ho': 'Ho_3', 'Ge': 'Ge_d_GW','Gd': 'Gd', 'Au': 'Au_GW', 'Zr': 'Zr_sv_GW',
        'Ga': 'Ga_d_GW', 'In': 'In_d_GW', 'Cs': 'Cs_sv_GW', 'Cr': 'Cr_pv', 'Cu': 'Cu_pv', 'Er': 'Er_3'}
 



elemets_highest_orbital = {'H':'s', 'He':'s', 'Li':'s', 'Be':'s', 'B':'p', 'C':'p',
                          'N':'p','O':'p','F':'p','Ne':'p','Na':'s','Mg':'s','Al':'p',
                          'Si':'p','P':'p','S':'p','Cl':'p','Ar':'p','K':'s','Ca':'s',
                          'Sc':'d','Ti':'d','V':'d','Cr':'d','Mn':'d','Fe':'d','Co':'d',
                          'Ni':'d','Cu':'d','Zn':'d','Ga':'d','Ge':'d','As':'d','Se':'d',
                          'Br':'d','Kr':'d','Rb':'s','Sr':'s','Y':'d','Zr':'d','Nb':'d',
                          'Mo':'d','Tc':'d','Ru':'d','Rh':'d','Pd':'d','Ag':'d', 'Cd':'d',
                          'In':'d','Sn':'d','Sb':'d','Te':'d','I':'d','Xe':'d','Cs':'s',
                          'Ba':'s','La':'d','Hf':'f','Ta':'f','W':'f','Re':'f','Os':'f',
                          'Ir':'f','Pt':'f','Au':'f','Hg':'f','Tl':'f','Pb':'f','Bi':'f',
                          'Po':'f','At':'f','Rn':'f','Fr':'s','Ra':'s','Ac':'d','Ce':'f',
                          'Pr':'f','Nd':'f','Sm':'f','Eu':'f','Gd':'f','Tb':'f','Dy':'f',
                          'Ho':'f','Er':'f','Tm':'f','Yb':'f','Lu':'f','Th':'d','Pa':'f',
                          'U':'f'}

elements_uvalue = { 'Co':3.32, 'Cr':3.7, 'Fe':5.3, 'Mn':3.9, 'Mo':4.38, 'Ni':6.45, 'V':3.25, 'W':6.2}


def chemical_formula_extract(line):
    chemical_formula = ''
    if '_chemical_formula_sum' in line:
        line_parsed = filter(None, line.strip().split(' '))

        x = ' '.join(line_parsed[1:])
        y = list(x)
        try:
          y[x.index("'")] = ""
          chemical_formula = ''.join(y[1:y.index("'")])
          chemical_formula = chemical_formula.replace("'", "").strip()
        except ValueError:
            print "_chemical_formula_sum should be a string"

        return chemical_formula

def cell_length_angle_extract(string_check, line):
    a = 0.0
    if string_check in line:
        line = filter(None, line.strip().split(' '))
        string_float = line[1]

        if '(' in string_float:
            string_float = string_float[:string_float.index('(')]
            if string_float.replace('.', '').isdigit():
                a = float(string_float)

            else:
                print "Error at %s" %(string_check)
        else:
            if string_float.replace('.', '').isdigit():
                a = float(string_float)
            else:
                print "Error at %s" %(string_check)

        return a


def calculate_distance(p1, p2):
    x_d = p1[0]-p2[0]
    y_d = p1[1]-p2[1]
    z_d = p1[2]-p2[2]

    return sqrt(x_d*x_d + y_d*y_d + z_d*z_d)

def choose_best_pseudo_potential(pseudos_list):
    choosen = []
    choosen_pseudo_enmax = pseudos_list[0]
    mean = 0.0
    for i in choosen_pseudo_enmax:
        mean += i[1]

    for pseudo in pseudos_list:
        p_mean = 0.0
        for i in pseudo:
            p_mean += i[1]
        if p_mean > mean:
            mean = p_mean
            choosen_pseudo_enmax = pseudo

    for i in choosen_pseudo_enmax:
        choosen.append(i[0])
    return choosen


def poscar_to_potcar(p_file):

    all_pseudo_potentials = collections.OrderedDict()
    pseudo_potentials_enmax = collections.OrderedDict() # gonna use it
    pseudo_potentials_chosen = []
    #poscar_file.close()
    pos_file = open(p_file, 'r')

    count = 1
    line = pos_file.readline()
    chemical_elements = []
    chemical_elements_changed = []
    while line:
        if count == 6:
            chemical_elements = filter(None, line.strip().split(' '))
            break
        line = pos_file.readline()
        count += 1
    pos_file.close()

    # getting all the pseudo-potentials for atoms in the element
    for element in chemical_elements:
        all_pseudo_potentials.update({element:[]})
        for pseudo in all_elements[element]:
            all_pseudo_potentials[element].append(pseudo)




    # Reading INAP file
    PSP_LIST = False
    PSP_NAME_bool = False
    PSP_NAME = []
    PSP_ENMAX = 0
    PSP_DENMAX = 0

    inap_file = open('AP_INPUT', 'r')
    line = inap_file.readline()
    while line:
        if 'PSP_LIST' in line:
            if line[0] != "#":
                line_list = filter(None, line.strip().split(' '))
                bool_value = line_list[2].lower().title()
                bool_value = bool_value[bool_value.index('.')+1:]
                bool_value = bool_value[:bool_value.index('.')]
                try:
                    PSP_LIST = ast.literal_eval(bool_value)
                except ValueError:
                    print 'PSP_LIST must be a boolean'
        if 'PSP_NAME' in line:
            if line[0] != "#":
                PSP_NAME_bool = True
                line = line[line.index("=")+1:]
                PSP_NAME = [i.strip() for i in filter(None, line.strip().split(','))]
        if 'PSP_ENMAX' in line:
            line_list = filter(None, line.strip().split())
            PSP_ENMAX = float(line_list [2])
        if 'PSP_DENMAX' in line:
            line_list = filter(None, line.strip().split())
            PSP_DENMAX = float(line_list [2])

        line = inap_file.readline()


    # Choosing pseudo-potentials
    if PSP_LIST:
        for element in chemical_elements:
            pseudo_potentials_chosen.append(preferred_psp[element])
    elif PSP_NAME_bool:
        for i in range(len(PSP_NAME)):
            is_i_legit = False
            for element in chemical_elements:
                if PSP_NAME[i] in all_elements[element]:
                    is_i_legit = True
                    break
            if is_i_legit:
                pseudo_potentials_chosen.append(PSP_NAME[i])
            else:
                sys.exit("%s is not one of the pseudo-potentials of %s" % (PSP_NAME[i], chemical_elements[i]))
        if len(pseudo_potentials_chosen) != len(chemical_elements):
            sys.exit("Not enough pseudo-potentials")

    else:
        if (PSP_ENMAX > 0) and (PSP_DENMAX > 0):
            for i in chemical_elements:
                pseudo_potentials_enmax.update({i:[]})
                for pseudo in all_elements[i]:
                    path = os.path.abspath('')
                    path += '/POTCARS/'+pseudo+'/POTCAR'
                    open_file = open(path, 'r')
                    line = open_file.readline()
                    while line:
                        strings_line = filter(None, line.strip().split(' '))
                        if (len(strings_line)>0) and (strings_line[0] == 'ENMAX'):
                            en_max = re.findall("\d+.\d*", strings_line[2])
                            f_en_max = float(en_max[0])
                            if (f_en_max >= PSP_ENMAX):
                                pseudo_potentials_enmax[i].append((pseudo, f_en_max))
                                break
                            break


                        line = open_file.readline()
            permutation_list = []
            for i in pseudo_potentials_enmax.keys():
                permutation_list.append(pseudo_potentials_enmax[i])
            permutation_list = list(itertools.product(*permutation_list))


            for pseudos in permutation_list:
                combinations_list = itertools.combinations(pseudos,2)
                retain = True
                for i in combinations_list:
                    if abs(i[0][1]-i[1][1]) <= PSP_DENMAX:
                        retain = False
                        break

                if retain:
                    pseudo_potentials_chosen.append(pseudos)
            pseudo_potentials_chosen = our_module.choose_best_pseudo_potential(pseudo_potentials_chosen)



        else:
            for element in chemical_elements:
                pseudo_potentials_chosen.append(preferred_psp[element])


    return pseudo_potentials_chosen



# it will also copy INCAR file in this folder
def create_compound_folder_poscar(poscar):
    ids_file_r = open('ids', 'r')
    id = int (ids_file_r.readline().strip())
    ids_file_r.close()

    chemical_elements = []
    indices = []
    chemical_elements_indices = collections.OrderedDict()
    common_divisor = 1
    folder_name = ''
    sub_folder_name = ''
    pos_file = open(poscar, 'r')
    line = pos_file.readline()
    count = 1
    while line:
        if count == 6:
            chemical_elements = filter(None, line.strip().split(' '))

        if count == 7:
            indices = [int(i) for i in filter(None, line.strip().split(' '))]
        line = pos_file.readline()
        count += 1
    pos_file.close()

    common_divisor = reduce(gcd, indices)
    indices = [int(i/common_divisor) for i in indices]

    if len(indices) == len(chemical_elements):
        for i in range(len(indices)):
            chemical_elements_indices.update({chemical_elements[i]:indices[i]})
    else:
        print "Check POSCAR format"
        exit(1)

    chemical_elements_indices = SortedDict(chemical_elements_indices)


    for key in chemical_elements_indices.keys():
        folder_name += key
        sub_folder_name += str(chemical_elements_indices[key])+'_'
    sub_folder_name = sub_folder_name[:-1]


    path = os.path.abspath('')
    path += '/'+folder_name+'/'+sub_folder_name
    digit = ''
    if os.path.exists(path):
        # folder numbering
        th_folder_file = open(path+'/th_folder', 'r')
        digit = int(th_folder_file.readline().strip())
        th_folder_file.close()
        th_folder_file = open(path+'/th_folder', 'w')
        digit += 1
        if digit>9:
            digit = str(digit)
        else:
            digit = '0'+str(digit)
        th_folder_file.write(digit)
        th_folder_file.close()
        path += '/'+digit
        os.system('mkdir -p '+path)

    else:
        os.system('mkdir -p '+path)

        # create id file
        f = open(path+'/id', 'w')
        f.write(str(id))
        f.close()

        # folder numbering
        th_folder_file = open(path+'/th_folder', 'w')
        digit = '0'+str(1)
        th_folder_file.write(digit)
        th_folder_file.close()
        path += '/'+digit
        os.system('mkdir -p '+path)

        ids_file_w = open('ids', 'w')
        ids_file_w.write(str(id+1))
        ids_file_w.close()

    # creating POTCAR
    pseudos = poscar_to_potcar(poscar)
    if len(chemical_elements) == len(pseudos):
        open(path+'/POTCAR', 'w').close()
        path_2 = os.path.abspath('')
        command = 'cat '
        for elem in pseudos:
            command += path_2+'/POTCARS/'+elem+'/POTCAR '
        command += '>> '+path+'/POTCAR'
        os.system(command)
    else:
        print "Less pseudo-potentials"
        exit(1)
    os.system('cp '+ poscar + ' '+path)
    os.system('mv '+path+'/'+poscar+' '+path+'/POSCAR')

    # writing KPOINTS file
    make_kpoints(path)

    # copy INCAR_default file to this folder just created
    os.system('cp INCAR_default'+ ' '+path)
    os.system('mv '+path+'/INCAR_default'+' '+path+'/INCAR')

    # copy bash.sh file
    os.system('cp runjob.sh'+ ' '+path)

    # change INCAR file in the path
    our_module.change_incar(path, poscar)

    compounds_file = open('compound_directories', 'a')
    compounds_file.write(path+'\n')

    return path

def create_compound_folder_cif(poscar):
    path = create_compound_folder_poscar(poscar)
    os.system('mv '+ poscar + ' '+path)
    os.system('mv '+path+'/'+poscar+' '+path+'/POSCAR')


def make_kpoints(path):
    f_read = open(path+'/POSCAR', 'r')
    line = f_read.readline().strip()
    count = 1
    # useful variables
    nbr_of_atoms = 0
    a = 0
    b = 0
    c = 0
    first_line_comment = ''
    while line:
        if count == 1:
            first_line_comment = line

        if count == 3:
            line = [float(i) for i in filter(None, line.split())]
            for i in line:
                a += i*i
            a = sqrt(a)
        if count == 4:
            line = [float(i) for i in filter(None, line.split())]
            for i in line:
                b += i*i
            b = sqrt(b)
        if count == 5:
            line = [float(i) for i in filter(None, line.split())]
            for i in line:
                c += i*i
            c = sqrt(c)
        if count == 7:
            line = [int(i) for i in filter(None, line.split())]
            for i in line:
                nbr_of_atoms += i

        count += 1
        line = f_read.readline().strip()
    f_read.close()

    # then we start applying formula
    k = int(1000/nbr_of_atoms)
    # reciprocal space
    r_a, r_b, r_c = int(100/a), int(100/b), int(100/c)


    k = (k/(r_a*r_b*r_c))**(1/3.0)
    r_a, r_b, r_c = [str(i) for i in (int(r_a*k), int(r_b*k), int(r_c*k))]


    kpoint_file = open(path+'/KPOINTS', 'w')
    kpoint_file.write(first_line_comment+'\n')
    kpoint_file.write('0\n')
    kpoint_file.write('Monkhorst\n')
    kpoint_file.write(r_a+' '+ r_b+' '+r_c+'\n')
    kpoint_file.write('0 0 0')
    kpoint_file.close()


# just to be used once
def write_default_incar(file_name, new_file_directory):
    r_file = open(file_name, 'r')
    w_file = open(new_file_directory, 'w')
    lines = r_file.readlines()
    w_file.writelines(lines)

    r_file.close()
    w_file.close()

def change_incar(path, poscar):

    # change INCAR file with AP_INCAR_CHANGE tags
    ap_input_file = open('AP_INPUT', 'r')
    line = ap_input_file.readline()
    LDAU_present = False
    LDAU = True
    INCAR_change = False
    AP_INCAR_CHANGE_tags = {}
    while line:
        line = line.strip()
        if 'AP_INCAR_change' in line:
            if line[0] != "#":
                line_list = line.split('=')
                line_list = filter(None, line_list)
                bool_value = line_list[1]
                bool_value = bool_value[bool_value.index('.')+1:]
                bool_value = bool_value[:bool_value.index('.')]
                bool_value = bool_value.lower().title()
                try:
                    INCAR_change = ast.literal_eval(bool_value) # change boolean string into a boolean

                except ValueError:
                    print 'AP_INCAR_change must follow the format'

        line = ap_input_file.readline()


    if INCAR_change:
        AP_INCAR_file = open('AP_INCAR_CHANGE','r')
        lines = AP_INCAR_file.readlines()
        for line in lines:
            if len(line.split('=')) > 1:
                line_list = line.strip().split('=')
                AP_INCAR_CHANGE_tags.update({line_list[0].strip():line_list[1].strip()})
        AP_INCAR_file.close()



    # if INCAR_change:
    #     AP_INCAR_file = open('AP_INCAR_CHANGE','r')
    #     lines = AP_INCAR_file.readlines()
    #
    #     count_ap = 0
    #     for line in lines:
    #         line_list = filter(None, line.strip().split('='))
    #         if 'LDAU' in line_list:
    #             LDAU_present = True
    #             bool_value = line_list[1].lower().title()
    #             bool_value = bool_value[bool_value.index('.')+1:]
    #             bool_value = bool_value[:bool_value.index('.')]
    #             try:
    #                 LDAU = ast.literal_eval(bool_value)
    #             except ValueError:
    #                 print 'LDAU must be a boolean'
    #         # the LDAU value that is in AP_INCAR_CHANGE has high precedence, we first checked it's value before any
    #         # change related to LDAU value is made in INCAR file
    #
    #         count_ap += 1
    #
    #     AP_INCAR_file.close()




    POSCAR_file = open(poscar, 'r')
    INCAR_file = open(path+'/INCAR', 'r')

    # extracting useful info from poscar
    count_pos = 1
    elements = []
    elements_nbr = collections.OrderedDict()
    elements_nbr_list = []
    quantum_nbrs = {'s':0, 'p':1, 'd':2, 'f':3}
    line = POSCAR_file.readline()
    while line:
        if count_pos == 6:
            line = filter(None, line.strip().split(' '))
            for i in line:
                elements.append(i)
        if count_pos == 7:
            line = filter(None, line.strip().split(' '))
            c = 0
            for i in line:
                elements_nbr.update({elements[c]:i})
                c += 1
            break

        line = POSCAR_file.readline()
        count_pos += 1
    POSCAR_file.close()

    for key in elements_nbr.keys():
        elements_nbr_list.append(key)
        elements_nbr_list.append(elements_nbr[key])


    #changing INCAR file with POSCAR file
    lines = INCAR_file.readlines()
    count_inc = 0
    for line in lines:
        line_list = filter(None, line.strip().split(' '))
        if 'SYSTEM' in line_list:
            line_change = 'SYSTEM = '+ ' '.join(elements_nbr_list)
            lines[count_inc] = line_change


        if 'MAGMOM' in line_list:
            magmom = 'MAGMOM = '
            for i in elements_nbr.keys():
                magmom += elements_nbr[i]+'*0.6 '
            magmom += '\n'
            lines[count_inc] = magmom

#        if 'LMAXMIX' in line_list:
#            lmaxmix = 'LMAXMIX = '
#            min = 0
#            for elem in elements:
#                if quantum_nbrs[elemets_highest_orbital[elem]] > min:
#                    min = quantum_nbrs[elemets_highest_orbital[elem]]
#
#            lmaxmix += str(2*min)+'\n'
#            lines[count_inc] = lmaxmix

        if 'LDAU' in line_list:
            bool_value = line_list[2].lower().title()
            bool_value = bool_value[bool_value.index('.')+1:]
            bool_value = bool_value[:bool_value.index('.')]
            if not LDAU_present:
                try:
                    LDAU = ast.literal_eval(bool_value)
                except ValueError:
                    print 'LDAU must be a boolean'
            ldau = 'LDAU = .'
            ldau += str(LDAU).upper()+'.\n'
            lines[count_inc] = ldau

        if 'LDAUL' in line_list:
            ldaul = 'LDAUL = '
            for i in elements:
                ldaul += str(-1) +' '
            ldaul += '\n'
            lines[count_inc] = ldaul
        if 'LDAUU' in line_list:
            ldauu = 'LDAUU = '
            for i in elements:          
                ldauu += str(0.0) +' '
            ldauu += '\n'
            lines[count_inc] = ldauu
        if 'LDAUJ' in line_list:
            ldauj = 'LDAUJ = '
            for i in elements:
                ldauj += str(0.0)+' '
            ldauj += '\n'
            lines[count_inc] = ldauj

        if 'ENCUT' in line_list:
            encut = 'ENCUT = '
            cutmax = 0.0
            pseudos = poscar_to_potcar(poscar)
            for pseudo in pseudos:
                path2 = os.path.abspath('')
                path2 += '/POTCARS/'+pseudo+'/POTCAR'
                open_file = open(path2, 'r')
                lines2 = open_file.readlines()

                for line in lines2:
                    strings_line = filter(None, line.strip().split(' '))
                    if (len(strings_line)>0) and (strings_line[0] == 'ENMAX'):
                        en_max = re.findall("\d+.\d*", strings_line[2])
                        f_en_max = float(en_max[0])
                        if (f_en_max > cutmax):
                                cutmax = f_en_max
                        break
            encut += str(cutmax) + '\n'
            lines[count_inc] = encut


        # handling AP_INCAR_CHANGE
        if INCAR_change:
            if len(line_list)>0 and (line_list[0].strip() in AP_INCAR_CHANGE_tags.keys()):
                tag_change = line_list[0].strip() + ' = '+AP_INCAR_CHANGE_tags[line_list[0].strip()] + '\n'
                lines[count_inc] = tag_change


        count_inc += 1
    INCAR_file.close()


    INCAR_file = open(path+'/INCAR', 'w')
    INCAR_file.writelines(lines)

    INCAR_file.close()


#def change_vasp_path():
#    count = 0
#    runjob_file = open('runjob.sh', 'r')
#    lines = runjob_file.readlines()
#    runjob_file.close()
#    for line in lines:
#        line_list = line.split("=")
#        if len(line_list)>1 and (line_list[0].strip() == 'VASP_EXE'):
#            AP_INPUT_file = open('AP_INPUT','r')
#            lines2 = AP_INPUT_file.readlines()
#            AP_INPUT_file.close()
#            for line2 in lines2:
#                line_list2 = line2.split("=")
#                if len(line_list2)>1 and (line_list2[0].strip() == 'VASP_LOCATION'):
#                    vasp_loc = line_list[0].strip() + '='+line_list2[1].strip()+'\n'
#                    lines[count] =  vasp_loc
#                    break
#            break
#        count += 1

#    runjob_file_w = open('runjob.sh', 'w')
#    runjob_file_w.writelines(lines)
#    runjob_file_w.close()




