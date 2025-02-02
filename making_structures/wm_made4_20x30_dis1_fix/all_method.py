import numpy as np
from ase.io import read, write
import os
import re
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit

A = 0.529177249
O = 15.99491**(1/2)
H = 1.00782**(1/2)
E = 27.2114079527

def modify_list(lst, factor_O, factor_H):
    result = [lst[i] * (factor_O if i < 3 else factor_H) for i in range(9)]
    return result

def har_E(dis_vec, force_3Nx3N):
    cQT = np.array(dis_vec)
    cQ = np.transpose(cQT)
    pro1 = np.dot(force_3Nx3N, cQ)
    E_har_hartree = (1/2) * np.dot(cQT, pro1)
    E_har_eV = E_har_hartree * E
    A = float(E_har_eV[0,0])
    X = f'{A:.10f}'
    return(str(X))

def cubic_spline_interpolator(x, y):
    # Create a CubicSpline interpolation
    cs = CubicSpline(x, y)
    def interpolate(x_interp):
        # Interpolate the y-value at the specified x_interp point
        y_interp = cs(x_interp)
        return y_interp
    return interpolate

def int_fn(source_xyz):
    a = read(source_xyz, index = ':')
    E_abs_eV = [i.get_total_energy()*E for i in a]
    E_eV = np.array([energy - E_abs_eV[9] for energy in E_abs_eV])
    x_coeff = np. array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    int_func = cubic_spline_interpolator(x_coeff, E_eV)
    return int_func

def morse(r, D_e, a, r_e):
    return (D_e * (np.exp(-2*a*(r-r_e))-2*np.exp(-a*(r-r_e))) + D_e)

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

###########################################################
################read the info from log file################
with open( './data/raw.c1.int.log' , 'r') as f:
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

###########################################################
###########################Method##########################
def direct(xyz_file):
    IN = read(xyz_file ,index = ':')
    E_input_eV=[i.get_total_energy() for i in IN]
    return E_input_eV

def harmonic(xyz_file, init_coordinate_bohr=init_coordinate_bohr, force_3Nx3N=force_3Nx3N):
    IN = read(xyz_file ,index = ':')
    E_input_eV=[i.get_total_energy() for i in IN]
    #extracting the positions from wm_.xyz input file, R_ang in ang
    R_ang=[np.reshape(np.array(i.get_positions()),9) for i in IN]
    #convert positions from wm_.xyz to bohr unit
    R_bohr_vec = [np.concatenate((r[:3]*(1/A), r[3:9]*(1/A))) for r in R_ang]
    #displacement vector, R-R_0, in bohr unit)
    dis_vec_nonmass = [r -init_coordinate_bohr for r in R_bohr_vec]
    e_har = [float(har_E(np.reshape(vec, (1,9)), force_3Nx3N)) for vec in dis_vec_nonmass]
    Delta_ML = ['{:.10f}'.format(float(input_energy - harmonic_energy)) 
                for input_energy, harmonic_energy in zip(E_input_eV, e_har)]
    return Delta_ML

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
    c_vec = np.array([[np.dot(vec, Q) for vec in dis_vec] for Q in Q_vec])
    int_funcs = [int_fn(f'./data/wm_{i}_0.1.xyz') for i in range(1, 4)]
    e_sum = np.sum([int_func(c) for int_func, c in zip(int_funcs, c_vec)], axis=0)
    Delta_ML = ['{:.10f}'.format(input_energy - normal_energy) 
                for input_energy, normal_energy in zip(E_input_eV, e_sum)]
    return Delta_ML

def internal(xyz_file, zmat_file):
    ##getting informations from data points for interpolation.
    a = read('./data/dis_1(2).xyz', index = ':')
    b = read('./data/ang_1(2).xyz', index = ':')
    E_eV_dis=[i.get_total_energy() for i in a]
    E_eV_ang=[i.get_total_energy() for i in b]
    ##origin data values and number of added data points
    distance_1 = [0.95883]
    angle_1 = [104.34240]
    num_values = 15
    distance_values = np.linspace(0.70, 1.40, num_values)
    angle_values = np.linspace(60.0, 150.0, num_values)
    distance_1.extend(distance_values)
    angle_1.extend(angle_values)
    dis_coeff = np.array(distance_1)[np.argsort(distance_1)]
    ang_coeff = np.array(angle_1)[np.argsort(angle_1)]
    E_ev_dis_array = np.array(E_eV_dis)[np.argsort(distance_1)]
    E_ev_ang_array = np.array(E_eV_ang)[np.argsort(angle_1)]
    ##cubicspline fitting
    int_func_dis = cubic_spline_interpolator(dis_coeff, E_ev_dis_array)
    int_func_ang = cubic_spline_interpolator(ang_coeff, E_ev_ang_array)
    #get dis_1, dis_2, ang_1 from zmat_file
    dis_1, dis_2, ang_1 = internal_coordinates(zmat_file)
    #get internal energy from cubicspline
    e_dis_1 = [int_func_dis(a) for a in dis_1]
    e_dis_2 = [int_func_dis(b) for b in dis_2]
    e_ang_1 = [int_func_ang(c) for c in ang_1]
    e_sum = [sum(values) for values in zip(e_dis_1, e_dis_2, e_ang_1)]

    IN = read(xyz_file ,index = ':')
    E_input_eV=[i.get_total_energy() for i in IN]
    Delta_ML = ['{:.10f}'.format(float(input_energy - internal_energy)) 
            for input_energy, internal_energy in zip(E_input_eV, e_sum)]
    return Delta_ML

def internal_morse(xyz_file, zmat_file):
    ##getting informations from data points for interpolation.
    a = read('./data/dis_1.xyz', index = ':')
    b = read('./data/ang_1.xyz', index = ':')
    E_dis=[i.get_total_energy()*E for i in a]
    E_ang=[i.get_total_energy()*E for i in b]
    E_eV_dis=[x - E_dis[0] for x in E_dis]
    E_eV_ang=[x - E_ang[0] for x in E_ang]
    ##origin data values and number of added data points
    distance_1 = [0.95883]
    angle_1 = [104.34240]
    num_values = 14
    distance_values = np.linspace(0.75, 1.30, num_values)
    angle_values = np.linspace(67.5, 147.0, num_values)
    distance_1.extend(distance_values)
    angle_1.extend(angle_values)
    dis_coeff = np.array(distance_1)[np.argsort(distance_1)]
    ang_coeff = np.array(angle_1)[np.argsort(angle_1)]
    E_ev_dis_array = np.array(E_eV_dis)[np.argsort(distance_1)]
    E_ev_ang_array = np.array(E_eV_ang)[np.argsort(angle_1)]
    ##morse potential fitting
    p0_dis=[1,1,distance_1[0]]
    p0_ang=[1,1,angle_1[0]]
    popt_dis, pcov_dis = curve_fit(morse, dis_coeff, E_ev_dis_array, p0=p0_dis)
    popt_ang, pcov_ang = curve_fit(morse, ang_coeff, E_ev_ang_array, p0=p0_ang)
    #get dis_1, dis_2, ang_1 from zmat_file
    dis_1, dis_2, ang_1 = internal_coordinates(zmat_file)
    #get internal energy from morse potential
    e_dis_1_mo = morse(dis_1, popt_dis[0], popt_dis[1], popt_dis[2])
    e_dis_2_mo = morse(dis_2, popt_dis[0], popt_dis[1], popt_dis[2])
    e_ang_1_mo = morse(ang_1, popt_ang[0], popt_ang[1], popt_ang[2])
    e_sum_mo = [sum(values) for values in zip(e_dis_1_mo,e_dis_2_mo,e_ang_1_mo)]

    IN = read(xyz_file ,index = ':')
    E_input_eV=[i.get_total_energy() for i in IN]
    Delta_ML_morse = ['{:.10f}'.format(float(input_energy - internal_energy)) 
            for input_energy, internal_energy in zip(E_input_eV, e_sum_mo)]
    return Delta_ML_morse

def delta_write(xyz_file, new_xyz_file, E_delta):
    data=[]
    with open(xyz_file, 'r') as file:
        num=0
        for line in file:
            if 'energy=' in line:
                line = line.replace(line.split('energy=')[1].split()[0], str(E_delta[num]))
                num += 1
            data.append(line)
        if num -1 != len(E_delta):
            ValueError('Length of Energy is not the same')
    with open(new_xyz_file, 'w') as f:
        f.writelines(data)
    print(f'{new_xyz_file} is successfully written')

def naming(xyz_file, method):
    base_filename = os.path.basename(xyz_file)
    file_log = f'data/{method[:3]}_{base_filename.split("_", 1)[-1]}'
    return file_log

def erase_force_info(xyz_file, new_xyz):
    a = read(xyz_file, index=':')
    for i in range(len(a)):
        a_copy = a[i].copy()
        del a_copy.arrays['forces']
        a[i] = a_copy
    write(new_xyz, a)
    print(f'{xyz_file}: force eliminated \n{new_xyz}: created')

########################################################
#########################excute#########################
word = 'combine'
xyz_file = f'./data/wm_{word}.xyz'
zmat_file = f'./data/wm_{word}.zmat'

erase_force_info('./made_4.xyz', xyz_file)

delta_write(xyz_file, naming(xyz_file, method='direct'), direct(xyz_file))
delta_write(xyz_file, naming(xyz_file, method='harmonic'), harmonic(xyz_file))
delta_write(xyz_file, naming(xyz_file, method='normal'), normal(xyz_file))
delta_write(xyz_file, naming(xyz_file, method='internal'), internal(xyz_file, zmat_file))
delta_write(xyz_file, naming(xyz_file, method='inmorse'), internal_morse(xyz_file, zmat_file))

