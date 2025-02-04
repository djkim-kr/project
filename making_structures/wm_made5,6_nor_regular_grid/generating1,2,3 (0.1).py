import random
import numpy as np
from ase.io import read
import os

A = 0.529177249
O = 15.9994**(1/2)
H = 1.00797**(1/2)
E = 27.2114079527


def make_float(x):
    x_float = float(x[0,0])
    return x_float
def R_format(position):
    symbols = ['O', 'H', 'H']
    output_str = ''
    for i, (symbol, pos) in enumerate(zip(symbols, position)):
        output_str += "{:2s}   {:14.8f}   {:14.8f}   {:14.8f}".format(
            symbol, pos[0], pos[1], pos[2]
        )
        if i < len(symbols) -1:
            output_str += '\n'
    return output_str
def modify_list1(lst, factor_O, factor_H):
    result = lst.copy()
    # Multiply factors for the components
    result[0] *= factor_O
    result[1] *= factor_O
    result[2] *= factor_O
    result[3] *= factor_H
    result[4] *= factor_H
    result[5] *= factor_H
    result[6] *= factor_H
    result[7] *= factor_H
    result[8] *= factor_H

    return result
def reshape_3X3(matrix_1X9):
    reshaped_matrix = matrix_1X9.reshape((3,3))
    with_atom = 'O: ' + str(reshaped_matrix[0]) + '\nH1: ' + str(reshaped_matrix[1]) + '\nH2: ' + str(reshaped_matrix[2])
    return with_atom

with open( 'data/raw.c1.int.log' , 'r') as f:
    lines = f.readlines()
    
    row3 = lines[2].split()
    row4 = lines[3].split()
    row5 = lines[4].split()

    O1 = [float(row3[2]), float(row3[3]), float(row3[4])]    
    H1 = [float(row4[2]), float(row4[3]), float(row4[4])]
    H2 = [float(row5[2]), float(row5[3]), float(row5[4])]

    init_coordinate = O1 + H1 + H2
    init_coordinate_bohr = np.mat(init_coordinate)

    # init_coordinate_amu = []
    # for i in range(3):
    #     init_coordinate_amu.append(init_coordinate[i]*O)
    # for i in range(3,9):
    #     init_coordinate_amu.append(init_coordinate[i]*H)
    # R0_amu_vec= np.array(init_coordinate_amu)
    # print(R0_amu_vec)

    # init_coordinate_ang = []
    # for i in range(len(init_coordinate)):
    #     init_coordinate_ang.append(init_coordinate[i]*A)
    # init_coordinate_ang_mat = np.array(init_coordinate_ang)

    #get Q_i from the raw file and covert to mass**(1/2) *bohr unit
    Q_1_bohr = []
    Q_2_bohr = []
    Q_3_bohr = []
    for i in range(23,32):
        if len(lines[i].split()) > 6: 
            Q_1_bohr.append(float(float(lines[i].split()[4])))
            Q_2_bohr.append(float(float(lines[i].split()[5])))
            Q_3_bohr.append(float(float(lines[i].split()[6])))
        else:
            Q_1_bohr.append(float(float(lines[i].split()[2])))
            Q_2_bohr.append(float(float(lines[i].split()[3])))
            Q_3_bohr.append(float(float(lines[i].split()[4])))
    Q_1_bohr_vec = np.array(Q_1_bohr)
    Q_2_bohr_vec = np.array(Q_2_bohr)
    Q_3_bohr_vec = np.array(Q_3_bohr)

    # print(Q_1_bohr_vec)
    # print(Q_2_bohr_vec)
    # print(Q_3_bohr_vec)

    # Q_1_vec = np.array(modify_list1(Q_1_bohr, O, H))
    # Q_2_vec = np.array(modify_list1(Q_2_bohr, O, H))
    # Q_3_vec = np.array(modify_list1(Q_3_bohr, O, H))

#made_5
c_1_list = np.linspace(-0.8, 0.8, 17)
c_2_list = np.linspace(-0.6,0.6, 13)
c_3_list = np.linspace(-0.8, 0.8, 17)
#made_6
c_1_list = np.linspace(-0.8, 0.8, 10)
c_2_list = np.linspace(-0.6,0.6, 8)
c_3_list = np.linspace(-0.8, 0.8, 10)


combined_list = []
for x in c_1_list:
    for y in c_2_list:
        for z in c_3_list:
            combined_list.append([x, y, z])
            

temp =[]
for i in range(len(combined_list)):
    temp_1 = (init_coordinate_bohr 
                + (combined_list[i][0] * Q_1_bohr_vec)
                + (combined_list[i][1] * Q_2_bohr_vec)
                + (combined_list[i][2] * Q_3_bohr_vec))
    temp_1_ang = temp_1 * A
    temp.append(temp_1_ang)

# print(temp[0])

# print(type(np.array(only_1[0].reshape((3,3)))))

file_log = "data/made_6.xyz"
if not os.path.exists("data"):
    os.makedirs("data")
with open(file_log, 'w') as f:
    for i in range(len(temp)):
        f.write('3\n')
        f.write(f'Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={np.round(combined_list[i][0],3)}/ {np.round(combined_list[i][1],3)}/ {np.round(combined_list[i][2],3)} pbc="T T T"\n')
        f.write('{}\n'.format(R_format(np.array(temp[i].reshape((3,3))))))        
print('Results written to {}'.format(file_log))