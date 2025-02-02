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
def har_E(dis_vec, force_3Nx3N):
    cQT = np.array(dis_vec)
    cQ = np.transpose(cQT)
    pro1 = np.dot(force_3Nx3N, cQ)
    E_har_hartree = (1/2) * np.dot(cQT, pro1)
    E_har_eV = E_har_hartree * E
    A = make_float(E_har_eV)
    X = '{:.10f}'.format(A)
    return(str(X))

filename_1 = str(input('enter the wm_xyz file to read : '))
energy_boundary = 50
# filename_1 = './data/wm_train.xyz'

#finding coefficients
with open( './data/raw.c1.int.log' , 'r') as f:
    lines = f.readlines()
#get R_0 from the raw file and convert tot mass**(1/2) * bohr unit
    row3 = lines[2].split()
    row4 = lines[3].split()
    row5 = lines[4].split()
    O1 = [float(row3[2]), float(row3[3]), float(row3[4])]    
    H1 = [float(row4[2]), float(row4[3]), float(row4[4])]
    H2 = [float(row5[2]), float(row5[3]), float(row5[4])]
    init_coordinate = O1 + H1 + H2
    init_coordinate_bohr = np.array(init_coordinate)
    #R0_amu_vec는 mass**(1/2) * bohr unit
    init_coordinate_amu = []
    for i in range(3):
        init_coordinate_amu.append(init_coordinate[i]*O)
    for i in range(3,9):
        init_coordinate_amu.append(init_coordinate[i]*H)
    R0_amu_vec= np.array(init_coordinate_amu)

    # #init_coordniate_ang는 angstrom unit
    init_coordinate_ang = []
    for i in range(len(init_coordinate)):
        init_coordinate_ang.append(init_coordinate[i]*A)
    init_coordinate_ang_mat = np.array(init_coordinate_ang)

    for i in range(38,54):
        read_row = lines[i]
        numbers = re.findall(r'-?\d+\.\d+', read_row)
        
        new_list_name = f"temp_{i}"
        globals()[new_list_name] = []
        for m in range(len(numbers)): 
            globals()[new_list_name].append(float(numbers[m]))

    row1 = temp_38 + temp_44[0:3]
    row2 = temp_39 + temp_45[0:3]
    row3 = temp_40 + temp_46[0:3]
    row4 = temp_41 + temp_44[3:6]
    row5 = temp_42 + temp_45[3:6]
    row6 = temp_43 + temp_46[3:6]
    row7 = temp_44 + temp_51
    row8 = temp_45 + temp_52
    row9 = temp_46 + temp_53
    #force_3Nx3N은 non-mass weighted hessian matrix (bohr unit)
    force_3Nx3N = np.mat([row1,row2,row3,row4,row5,row6,row7,row8,row9])

#get Q_i from the raw file and covert to mass*(1/2)*bohr unit
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
    #Q_1 vec는 mass**(1/2)*bohr unit
    Q_1_vec = np.array(modify_list1(Q_1_bohr, O, H))
    Q_2_vec = np.array(modify_list1(Q_2_bohr, O, H))
    Q_3_vec = np.array(modify_list1(Q_3_bohr, O, H))

#extracting the energies from wm_.xyz input file, E_input in eV
IN = read(filename_1 ,index = ':')
E_input_eV=[]
for i in IN:
    E_input_eV.append(i.get_total_energy())
#extracting the positions from wm_.xyz input file, R_ang in ang
tempx=[]
for i in IN:
    tempx.append(i.get_positions())
tempxx = [np.array(b) for b in tempx]
R_ang = [np.reshape(b, 9) for b in tempxx]
#convert positions from wm_.xyz to have mass**(1/2) * bohr unit
R_ang_amu_vec=[]
for i in range(len(R_ang)):
    temp_a = R_ang[i][0:3]*(1/A)*O
    temp_b = R_ang[i][3:9]*(1/A)*H
    R_ang_amu_vec.append(np.concatenate((temp_a,temp_b)))
#convert positions from wm_.xyz to bohr unit
R_bohr_vec = []
for i in range(len(R_ang)):
    temp_a = R_ang[i][0:3]*(1/A)
    temp_b = R_ang[i][3:9]*(1/A)
    R_bohr_vec.append(np.concatenate((temp_a,temp_b)))
#displacement vector, R-R_0, in amu unit
dis_vec=[]
for i in range(len(R_ang)):
    dis_vec.append(R_ang_amu_vec[i]-R0_amu_vec)
#displacement vector, R-R_0, in bohr unit)
dis_vec_nonmass = []
for i in range(len(R_ang)):
    dis_vec_nonmass.append(R_bohr_vec[i]-init_coordinate_bohr)
#get the coefficients for each vibration modes for each coordinates
c_1= []
c_2= []
c_3= []
for i in range(len(dis_vec)):
    c_1.append(np.dot(dis_vec[i],Q_1_vec))
    c_2.append(np.dot(dis_vec[i],Q_2_vec))
    c_3.append(np.dot(dis_vec[i],Q_3_vec))

#############################################################
#get the interpolation function for each vibration mode
source_1 = './data/wm_1_0.1.xyz' 
source_2 = './data/wm_2_0.1.xyz' 
source_3 = './data/wm_3_0.1.xyz' 
a = read(source_1 ,index = ':')
b = read(source_2 ,index = ':')
c = read(source_3 ,index = ':')

E_abs_eV_1=[]
E_abs_eV_2=[]
E_abs_eV_3=[]
for i in a:
    E_abs_eV_1.append(i.get_total_energy()*E)
for i in b:
    E_abs_eV_2.append(i.get_total_energy()*E)
for i in c:
    E_abs_eV_3.append(i.get_total_energy()*E)
E_eV_1 = []
E_eV_2 = []
E_eV_3 = []
# for i in range(len(E_abs_eV_1)):
#     E_eV_1.append(E_abs_eV_1[i] - E_abs_eV_1[5])
#     E_eV_2.append(E_abs_eV_2[i] - E_abs_eV_2[5])
#     E_eV_3.append(E_abs_eV_3[i] - E_abs_eV_3[5])  
# E_ev_1_array = np.array(E_eV_1)
# E_ev_2_array = np.array(E_eV_2)
# E_ev_3_array = np.array(E_eV_3)

# x_coeff = np.array([-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
# int_func_1 = cubic_spline_interpolator(x_coeff, E_ev_1_array)
# int_func_2 = cubic_spline_interpolator(x_coeff, E_ev_2_array)
# int_func_3 = cubic_spline_interpolator(x_coeff, E_ev_3_array)
for i in range(len(E_abs_eV_1)):
    E_eV_1.append(E_abs_eV_1[i] - E_abs_eV_1[9])
    E_eV_2.append(E_abs_eV_2[i] - E_abs_eV_2[9])
    E_eV_3.append(E_abs_eV_3[i] - E_abs_eV_3[9])  
E_ev_1_array = np.array(E_eV_1)
E_ev_2_array = np.array(E_eV_2)
E_ev_3_array = np.array(E_eV_3)

x_coeff = np.array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
int_func_1 = cubic_spline_interpolator(x_coeff, E_ev_1_array)
int_func_2 = cubic_spline_interpolator(x_coeff, E_ev_2_array)
int_func_3 = cubic_spline_interpolator(x_coeff, E_ev_3_array)

#get the energy for each vibraton mode for each coordinates
e_1 = []
e_2 = []
e_3 = []
e_sum = []
for i in range(len(R_ang)):
    e_1.append(int_func_1(c_1[i]))
    e_2.append(int_func_2(c_2[i]))
    e_3.append(int_func_3(c_3[i]))
    e_sum.append(int_func_1(c_1[i])+int_func_2(c_2[i])+int_func_3(c_3[i]))

#get \Delta_ML. input energy - e_sum 과 같은 값.
Delta_ML = []
for i in range(len(R_ang)):
    Delta_ML.append('{:.10f}'.format(float(E_input_eV[i]-e_sum[i])))

##############################################################
#cacluate E_1,2,3 with harmonic way so that to compare. 
#e_har_1,2,3은 1,2,3 vector에 대한 계산값, 123은 1,2,3 더한 후 계산한 값, har은 1~9 모두 포함된 값.
e_har_1 = []
e_har_2 = []
e_har_3 = []
e_har_123 = []
e_har = []
for i in range(len(R_ang)):
    e_har_1.append(float(har_E(np.reshape(c_1[i]*Q_1_bohr_vec, (1,9)), force_3Nx3N)))
    e_har_2.append(float(har_E(np.reshape(c_1[i]*Q_2_bohr_vec, (1,9)), force_3Nx3N)))
    e_har_3.append(float(har_E(np.reshape(c_1[i]*Q_3_bohr_vec, (1,9)), force_3Nx3N)))
    e_har_123.append(float(har_E(np.reshape(c_1[i]*Q_1_bohr_vec +c_2[i]*Q_2_bohr_vec +c_3[i]*Q_3_bohr_vec, (1,9)), force_3Nx3N)))
    e_har.append(float(har_E(np.reshape(dis_vec_nonmass[i], (1,9)), force_3Nx3N)))
Delta_ML_har = []
for i in range(len(R_ang)):
    Delta_ML_har.append('{:.10f}'.format(float(E_input_eV[i]-e_har[i])))

################################################################
#Delta_ML을 xyz file에 기존 coordinate와 함께 입히기. for QUIP, DE_{}.xyz file 로 만듦.
base_filename = os.path.basename(filename_1)

file_log = "data/nor_{}".format(base_filename.split("_", 1)[-1])
if not os.path.exists("data"):
    os.makedirs("data")
with open(file_log, 'w') as f:
    for i in range(len(R_ang)):
        if E_input_eV[i] < energy_boundary:
            f.write('3\n')
            f.write('Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={} pbc="T T T"\n'.format(Delta_ML[i]))
            f.write('{}\n'.format(R_format(R_ang[i].reshape((3,3)))))

print('Results written to {}'.format(file_log))

#기존 harmonic way에 대한 계산값을 만들기 input - har. for quip, har_{}.xyz file로 만듦.
file_log_2 = "data/har_{}".format(base_filename.split("_", 1)[-1])
with open(file_log_2, 'w') as f:
    for i in range(len(R_ang)):
        if E_input_eV[i] < energy_boundary:                
            f.write('3\n')
            f.write('Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={} pbc="T T T"\n'.format(Delta_ML_har[i]))
            f.write('{}\n'.format(R_format(R_ang[i].reshape((3,3)))))

print('Results written to {}'.format(file_log_2))
