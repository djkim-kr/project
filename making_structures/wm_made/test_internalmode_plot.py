import numpy as np
from ase.io import read
import os
import re
from scipy import interpolate
from scipy.interpolate import interp1d 
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt
#1 Hartree = 27.2114079527eV
#1 Bohr = 0.529177249 Angstrom
A = 0.529177249
O = 15.99491**(1/2)
H = 1.00782**(1/2)
E = 27.2114079527

word = str(input('enter word : '))
# word = 'train'  
energy_boundary=50

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
def cubic_spline_interpolator(x, y):
    # Create a CubicSpline interpolation
    cs = CubicSpline(x, y)

    def interpolate(x_interp):
        # Interpolate the y-value at the specified x_interp point
        y_interp = cs(x_interp)
        return y_interp

    return interpolate

##origin zmat
## O
## H      1      0.95883
## H      1      0.95883    2    104.34240



##getting informations from data points to make interpolation.
source_dis = './data/dis_1.xyz' 
source_ang = './data/ang_1.xyz'
a = read(source_dis ,index = ':')
b = read(source_ang ,index = ':')

E_abs_eV_dis=[]
E_abs_eV_ang=[]
for i in a:
    E_abs_eV_dis.append(i.get_total_energy()*E)
for i in b:
    E_abs_eV_ang.append(i.get_total_energy()*E)
E_eV_dis = []
E_eV_ang = []
for i in range(len(E_abs_eV_dis)):
    E_eV_dis.append(E_abs_eV_dis[i] - E_abs_eV_dis[0])
    E_eV_ang.append(E_abs_eV_ang[i] - E_abs_eV_ang[0])

##origin data values and number of added data points
distance_1 = [0.95883]
angle_1 = [104.34240]
num_values = 14
distance_values = np.linspace(0.75, 1.30, num_values)
angle_values = np.linspace(67.5, 147.0, num_values)
distance_1.extend(distance_values)
angle_1.extend(angle_values)
sorted_indices_dis = np.argsort(distance_1)
sorted_indices_ang = np.argsort(angle_1)

dis_coeff = np.array(distance_1)[sorted_indices_dis]
ang_coeff = np.array(angle_1)[sorted_indices_ang]
E_ev_dis_array = np.array(E_eV_dis)[sorted_indices_dis]
E_ev_ang_array = np.array(E_eV_ang)[sorted_indices_ang]

# print(dis_coeff)
# print(ang_coeff)
# print(E_ev_dis_array)
# print(E_ev_ang_array)

int_func_dis = cubic_spline_interpolator(dis_coeff, E_ev_dis_array)
int_func_ang = cubic_spline_interpolator(ang_coeff, E_ev_ang_array)

zmat_file = f'./data/wm_{word}.zmat'
# zmat_file = f'./data/origin.zmat'

##get dis, ang info. from wm_{}.zmat
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
# print(distance_1_list[0])
# print(distance_2_list[0])
# print(angle_1_list[0])

##interpolate the energy from dis, ang info. from each coord.
e_dis_1 = []
e_dis_2 = []
e_ang_1 = []
e_sum = []
for i in range(len(distance_1_list)):
    e_dis_1.append(int_func_dis(distance_1_list[i]))
    e_dis_2.append(int_func_dis(distance_2_list[i]))
    e_ang_1.append(int_func_ang(angle_1_list[i]))
    e_sum.append(int_func_dis(distance_1_list[i])+int_func_dis(distance_2_list[i])+int_func_ang(angle_1_list[i]))
    
##extracting infos from .xyz file
xyz_file = f'./data/wm_{word}.xyz'
# xyz_file = f'./data/origin.xyz'
##extracting the energies from wm_.xyz input file, E_input in eV
IN = read(xyz_file ,index = ':')
E_input_eV=[]
for i in IN:
    E_input_eV.append(i.get_total_energy())
##extracting the positions from wm_.xyz input file, R_ang in ang
tempx=[]
for i in IN:
    tempx.append(i.get_positions())
tempxx = [np.array(b) for b in tempx]
R_ang = [np.reshape(b, 9) for b in tempxx]


## get Delta_ML. input energy - e_sum 과 같은 값.
Delta_ML = []
for i in range(len(R_ang)):
    Delta_ML.append('{:.10f}'.format(float(E_input_eV[i]-e_sum[i])))

################################################################
#Delta_ML을 xyz file에 기존 coordinate와 함께 입히기. for QUIP, int_{}.xyz file 로 만듦.
# base_filename = os.path.basename(xyz_file)
# file_log = "data/int_{}".format(base_filename.split("_", 1)[-1])
# if not os.path.exists("data"):
#     os.makedirs("data")
# with open(file_log, 'w') as f:
#     for i in range(len(R_ang)):
#         if E_input_eV[i] < energy_boundary:
#             f.write('3\n')
#             f.write('Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={} pbc="T T T"\n'.format(Delta_ML[i]))
#             f.write('{}\n'.format(R_format(R_ang[i].reshape((3,3)))))

# print('Results written to {}'.format(file_log))

x_interp = np.linspace(0.75,1.30, 100)  # Interpolate between 1 and 10, with 100 points

# Interpolate the y-values at the specified points
y_interp = int_func_dis(x_interp)
# Plot the original data and the interpolated curve
plt.figure(figsize=(8, 6))
plt.scatter(dis_coeff, E_ev_dis_array, label='Original Data', color='red', marker='o')
plt.plot(distance_1_list, e_dis_1,'.b' , label='interpolated data')
plt.plot(x_interp, y_interp, label='CubicSpline Interpolation', color='black')
plt.xlabel('distance(Å)')
plt.ylabel('E_eV')
plt.legend()
plt.title('CubicSpline Interpolation_dis_1')
plt.grid(True)
plt.show()

# Plot the original data and the interpolated curve
plt.figure(figsize=(8, 6))
plt.scatter(dis_coeff, E_ev_dis_array, label='Original Data', color='red', marker='o')
plt.plot(distance_2_list, e_dis_2,'.b' , label='interpolated data')
plt.plot(x_interp, y_interp, label='CubicSpline Interpolation', color='black')
plt.xlabel('distance(Å)')
plt.ylabel('E_eV')
plt.legend()
plt.title('CubicSpline Interpolation_dis_2')
plt.grid(True)
plt.show()

x_interp_2 = np.linspace(67.5,147.0, 100) 
y_interp_2 = int_func_ang(x_interp_2)
# Plot the original data and the interpolated curve
plt.figure(figsize=(8, 6))
plt.scatter(ang_coeff, E_ev_ang_array, label='Original Data', color='red', marker='o')
plt.plot(angle_1_list, e_ang_1,'.b' , label='interpolated data')
plt.plot(x_interp_2, y_interp_2, label='CubicSpline Interpolation', color='black')
plt.xlabel('angle')
plt.ylabel('E_eV')
plt.legend()
plt.title('CubicSpline Interpolation_ang')
plt.grid(True)
plt.show()