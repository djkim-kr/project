import numpy as np
from ase.io import read

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240
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

def quip(word):
    energies=[]
    for i in range(1,6):
        quip_file=f'./data/quip_test/{word}/quip_test_{i}.xyz'
        e_temp=[i.get_total_energy() for i in read(quip_file ,index = ':')]
        energies.append(e_temp)
    energies_2=[]
    for i in range(1,6):
        quip_file=f'./data/quip_test_60/{word}/quip_test_{i}.xyz'
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


#made_5
c_1_list = np.linspace(-0.8, 0.8, 17)
c_2_list = np.linspace(-0.6,0.6, 13)
c_3_list = np.linspace(-0.8, 0.8, 17)
combined_list = []
for x in c_1_list:
    for y in c_2_list:
        for z in c_3_list:
            combined_list.append([x, y, z])

c_1 = [x[0] for x in combined_list]
c_2 = [x[1] for x in combined_list]
c_3 = [x[2] for x in combined_list]

dis_1, dis_2, ang_1 = internal_coordinates(zmat_file='../data/wm_combine5.zmat')

delta_dir, pot_dir_full, pot_dir_104, error_dir_abs_full, error_dir_abs_104, dir_error_diff = quip('dir')
delta_har, pot_har_full, pot_har_104, error_har_abs_full, error_har_abs_104, har_error_diff = quip('har')
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

np.savetxt(f'made_data.txt', data, delimiter=' ', header=header, comments='', fmt='%s')