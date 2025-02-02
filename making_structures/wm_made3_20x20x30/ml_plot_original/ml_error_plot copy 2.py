import numpy as np
from ase.io import read
import re

A = 0.529177249
O = 15.99491**(1/2)
H = 1.00782**(1/2)
E = 27.2114079527
#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240
def modify_list(lst, factor_O, factor_H):
    result = [lst[i] * (factor_O if i < 3 else factor_H) for i in range(9)]
    return result
def internal_coordinates(zmat_file):
    distance_1_list = []
    distance_2_list = []
    angle_1_list = []
    with open(zmat_file, 'r') as zmat_file:
        lines = zmat_file.readlines()
    block = []
    for line in lines:
        if line.strip() == 'O':
            if block:
                if len(block) == 3:
                    distance_1 = float(block[1].split()[2])
                    distance_2 = float(block[2].split()[2])
                    angle_1 = float(block[2].split()[4])
                    distance_1_list.append(distance_1)
                    distance_2_list.append(distance_2)
                    angle_1_list.append(angle_1)
                block = []  # Reset the block
        block.append(line.strip())
    # Process the last block (if it exists)
    if len(block) == 3:
        distance_1 = float(block[1].split()[2])
        distance_2 = float(block[2].split()[2])
        angle_1 = float(block[2].split()[4])
        distance_1_list.append(distance_1)
        distance_2_list.append(distance_2)
        angle_1_list.append(angle_1)
    return distance_1_list, distance_2_list, angle_1_list
with open( '../data/raw.c1.int.log' , 'r') as f:
    lines = f.readlines()
    #get R_0 from the raw file and convert tot mass**(1/2) * bohr unit
    atoms = [line.split()[2:5] for line in lines[2:5]]
    init_coordinate = np.array(atoms, dtype=float).flatten()
    init_coordinate_bohr = init_coordinate
    init_coordinate_amu = modify_list(init_coordinate_bohr, O, H)
    #R0_amu_vec는 mass**(1/2) * bohr unit
    R0_amu_vec= np.array(init_coordinate_amu)
    # #init_coordniate_ang는 angstrom unit
    init_coordinate_ang = init_coordinate * A

    temp_list = [list(map(float, re.findall(r'-?\d+\.\d+', lines[i]))) for i in range(38, 54)]
    #force_3Nx3N은 non-mass weighted hessian matrix (bohr unit)
    force_3Nx3N = np.mat([
        temp_list[0] + [temp_list[6][0]] + [temp_list[7][0]] + [temp_list[8][0]],
        temp_list[1] + [temp_list[6][1]] + [temp_list[7][1]] + [temp_list[8][1]],
        temp_list[2] + [temp_list[6][2]] + [temp_list[7][2]] + [temp_list[8][2]],
        temp_list[3] + [temp_list[6][3]] + [temp_list[7][3]] + [temp_list[8][3]],
        temp_list[4] + [temp_list[6][4]] + [temp_list[7][4]] + [temp_list[8][4]],
        temp_list[5] + [temp_list[6][5]] + [temp_list[7][5]] + [temp_list[8][5]],
        temp_list[6] + temp_list[13],
        temp_list[7] + temp_list[14],
        temp_list[8] + temp_list[15]
        ])
    
    #get Q_i from the raw file and covert to mass*(1/2)*bohr unit
    Q_bohr = [[] for _ in range(3)]
    for i in range(23, 32):
        line_data = lines[i].split()
        Q_bohr[0].append(float(line_data[4] if len(line_data) > 6 else line_data[2]))
        Q_bohr[1].append(float(line_data[5] if len(line_data) > 6 else line_data[3]))
        Q_bohr[2].append(float(line_data[6] if len(line_data) > 6 else line_data[4]))
    Q_bohr_vec = np.array(Q_bohr)
    #Q_vec는 mass**(1/2)*bohr unit
    Q_vec = np.array([modify_list(Q_bohr[i], O, H) for i in range(3)])
def normal(xyz_file, R0_amu_vec=R0_amu_vec, Q_vec=Q_vec):
    IN = read(xyz_file ,index = ':')
    #extracting the energies from wm_.xyz input file, E_input in eV
    E_input_eV=[i.get_total_energy() for i in IN]
    #extracting the positions from wm_.xyz input file, R_ang in ang
    R_ang=[np.reshape(np.array(i.get_positions()),9) for i in IN]
    #convert positions from wm_.xyz to have mass**(1/2) * bohr unit
    R_ang_amu_vec = [np.concatenate((r[:3]*(1/A)*O, r[3:9]*(1/A)*H)) for r in R_ang]
    #displacement vector, R-R_0, in amu unit
    dis_vec = [r -R0_amu_vec for r in R_ang_amu_vec]
    #get the coefficients for each vibration modes for each coordinates
    c_1 = np.array([np.dot(vec, Q_vec[0]) for vec in dis_vec])
    c_2 = np.array([np.dot(vec, Q_vec[1]) for vec in dis_vec])
    c_3 = np.array([np.dot(vec, Q_vec[2]) for vec in dis_vec])
    return c_1, c_2, c_3


def quip(word):
    energies=[]
    for i in range(1,6):
        quip_file=f'./data/quip_test/{word}/quip_test_{i}.xyz'
        e_temp=[i.get_total_energy() for i in read(quip_file ,index = ':')]
        energies.append(e_temp)
    energies_2=[]
    for i in range(1,6):
        quip_file=f'./data/quip_test_0104/{word}/quip_test_{i}.xyz'
        e_temp_2=[i.get_total_energy() for i in read(quip_file ,index = ':')]
        energies_2.append(e_temp_2)

    e_test=[i.get_total_energy() for i in read(f'./data/{word}_test.xyz', index=':')]

    energies=np.array(energies)
    energies_2=np.array(energies_2)

    e_quip=np.mean(energies, axis=0)
    e_quip_2=np.mean(energies_2, axis=0)

    ml_delta = [1000*x for x in e_test]
    ml_error_abs = [abs(1000*(a-b)) for a,b in zip(e_test,e_quip)]
    ml_error_2_abs = [abs(1000*(a-b)) for a,b in zip(e_test,e_quip_2)]
    ml_error_diff = [a-b for a,b in zip(ml_error_2_abs, ml_error_abs)]

    ml_pot=[x*1000 for x in e_quip]
    ml_pot_2=[x*1000 for x in e_quip_2]

    return ml_delta, ml_pot, ml_pot_2, ml_error_abs, ml_error_2_abs, ml_error_diff

c_1,c_2,c_3 = normal('../made_3.xyz')
dis_1, dis_2, ang_1 = internal_coordinates(zmat_file='../data/wm_combine.zmat')

delta_har, pot_har_full, pot_har_104, error_har_abs_full, error_har_abs_104, har_error_diff = quip('har')
delta_dir, pot_dir_full, pot_dir_104, error_dir_abs_full, error_dir_abs_104, dir_error_diff = quip('dir')
delta_nor, pot_nor_full, pot_nor_104, error_nor_abs_full, error_nor_abs_104, nor_error_diff = quip('nor')
delta_int, pot_int_full, pot_int_104, error_int_abs_full, error_int_abs_104, int_error_diff = quip('int')

data = np.column_stack((c_1, c_2, c_3,
                        dis_1, dis_2, ang_1,
                        delta_dir, delta_har, delta_nor, delta_int,
                        pot_dir_full, pot_har_full, pot_nor_full, pot_int_full,
                        pot_dir_104, pot_har_104, pot_nor_104, pot_int_104,
                        error_dir_abs_full, error_dir_abs_104, dir_error_diff,
                        error_har_abs_full, error_har_abs_104, har_error_diff,
                        error_nor_abs_full, error_nor_abs_104, nor_error_diff,
                        error_int_abs_full, error_int_abs_104, int_error_diff))

header = '''c_1 c_2 c_3 dis_1 dis_2 ang_1 delta_dir delta_har delta_nor delta_int pot_dir_full pot_har_full pot_nor_full pot_int_full pot_dir_104 pot_har_104 pot_nor_104 pot_int_104 error_dir_abs_full error_dir_abs_104 dir_error_diff error_har_abs_full error_har_abs_104 har_error_diff error_nor_abs_full error_nor_abs_104 nor_error_diff error_int_abs_full error_int_abs_104 int_error_diff'''
# print(error_dir_abs_104[:10])
# np.savetxt(f'original_data.txt', data, delimiter=' ', header=header, comments='', fmt='%s')