import re
import numpy as np
import os
from ase.io import read


energy_boundary=0.05
#file에서 정보 가져오기
filename = str(input('enter the xyz file to read : '))
# filename = 'wm/wm_test.xyz'
a = read(filename ,index = ':')
E_anh_eV=[]
for i in a:
    E_anh_eV.append(i.get_total_energy())

tempx=[]
for i in a:
    tempx.append(i.get_positions())
tempxx = [np.array(d) for d in tempx]
R_ang = [np.reshape(b, 9) for b in tempxx]


#matrix형식으로 되어 있는 DE 를 float으로 바꿔주는 함수
def make_float(x):
    x_float = float(x[0,0])
    return x_float
#1x9로 되어있는 함수를 3x3으로 바꾸는 것
def reshape_3X3(matrix_1X9):
    reshaped_matrix = matrix_1X9.reshape((3,3))
    with_atom = 'O: ' + str(reshaped_matrix[0]) + '\nH1: ' + str(reshaped_matrix[1]) + '\nH2: ' + str(reshaped_matrix[2])
    return with_atom
#atom의 coordinate를 특정 순서와 형식에 맞게 표현하는 것
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


#writing the results in the log file in data directory in the right format for quippy
base_filename = os.path.basename(filename)
file_log = f"wm_data/{energy_boundary}_{base_filename}"

if not os.path.exists("wm_data"):
    os.makedirs("wm_data")

with open(file_log, 'w') as f:

    for i in range(len(R_ang)):
        if E_anh_eV[i] < energy_boundary:
            f.write('3\n')
            f.write('Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={} pbc="T T T"\n'.format(E_anh_eV[i]))
            f.write('{}\n'.format(R_format(R_ang[i].reshape((3,3)))))

# Print a message indicating the file has been written
print('Results written to {}'.format(file_log))