import numpy as np
from ase.io import read
import os
import re
from scipy.interpolate import interp2d
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

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
def reshape_3X3(matrix_1X9):
    reshaped_matrix = matrix_1X9.reshape((3,3))
    with_atom = 'O: ' + str(reshaped_matrix[0]) + '\nH1: ' + str(reshaped_matrix[1]) + '\nH2: ' + str(reshaped_matrix[2])
    return with_atom
def co_combine(co_a, co_b):
    return [f'{a} / {b}' for a in co_a for b in co_b]
def cubic_spline_interpolator(x, y):
    # Create a CubicSpline interpolation
    cs = CubicSpline(x, y)

    def interpolate(x_interp):
        # Interpolate the y-value at the specified x_interp point
        y_interp = cs(x_interp)
        return y_interp

    return interpolate
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

    init_coordinate_amu = []
    for i in range(3):
        init_coordinate_amu.append(init_coordinate[i]*O)
    for i in range(3,9):
        init_coordinate_amu.append(init_coordinate[i]*H)
    R0_amu_vec= np.array(init_coordinate_amu)
    

    init_coordinate_ang = []
    for i in range(len(init_coordinate)):
        init_coordinate_ang.append(init_coordinate[i]*A)
    init_coordinate_ang_mat = np.array(init_coordinate_ang)

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

    Q_1_vec = np.array(modify_list1(Q_1_bohr, O, H))
    Q_2_vec = np.array(modify_list1(Q_2_bohr, O, H))
    Q_3_vec = np.array(modify_list1(Q_3_bohr, O, H))    

def plot_2B(source):
    source_read = read(source, index = ':')
    
    numbers = source.replace('.xyz', '').split('_')[-1].split('-')
    a = float(numbers[0])
    b = float(numbers[1])
    stra=numbers[0]
    strb=numbers[1]
    co_set = co_combine(co_lists[a], co_lists[b])
    co_a = np.array([float(item.split('/')[0]) for item in co_set])
    co_b = np.array([float(item.split('/')[-1]) for item in co_set])
    points = np.column_stack((co_a, co_b))

    E_eV_array = np.array([i.get_total_energy() for i in source_read])

    grid_x, grid_y = np.mgrid[min(co_a):max(co_a):complex(0, 100), 
                            min(co_b):max(co_b):complex(0, 100)]

    Z_new = griddata(points, E_eV_array, (grid_x, grid_y), method='cubic')

    plt.figure()
    plt.imshow(Z_new,
            extent=(min(co_lists[a])-0.1, max(co_lists[a])+0.1 ,min(co_lists[b])-0.1 ,max(co_lists[b])+0.1),
            origin='lower')
    plt.scatter(co_a, co_b, c=E_eV_array, edgecolor='k')  # 원본 데이터 포인트 표시
    plt.title(f'2B_{stra}_{strb}_interpolation')
    plt.xlabel(f'Coeff_{stra}')
    plt.ylabel(f'Coeff_{strb}')
    plt.colorbar(label='E(eV)')
    plt.savefig(f"./figure/2B_{stra}_{strb}_plot.png", dpi=300)
    print(f"./figure/2B_{stra}_{strb}_plot.png has been saved")
    
    # 3D 플롯 생성
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # 3D 표면 플롯 생성
    surf = ax.plot_surface(grid_x, grid_y, Z_new, cmap='viridis', edgecolor = 'none')
    # 원본 데이터 포인트 표시
    ax.scatter(co_a, co_b, E_eV_array, color='r', marker='o')

    ax.set_xlabel(f'Coeff_{stra}')
    ax.set_ylabel(f'Coeff_{strb}')
    ax.set_zlabel('E(eV)')
    ax.set_title(f'3D Surface of 2B_{stra}_{strb}_interpolation Interpolation')
    # 컬러바 추가
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.savefig(f"./figure/2B_{stra}_{strb}_3D.png", dpi=300)
    print(f"./figure/2B_{stra}_{strb}_3D.png has been saved")
def int_fn(source_xyz):
    a = read(source_xyz, index = ':')
    E_abs_eV = [i.get_total_energy()*E for i in a]
    E_eV = np.array([energy - E_abs_eV[9] for energy in E_abs_eV])
    x_coeff = np. array([-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8])
    int_func = cubic_spline_interpolator(x_coeff, E_eV)
    return int_func
def int_2B(word):
    source = f'./data/wm_{word}.xyz'
    source_read = read(source, index = ':')
    numbers = source.replace('.xyz', '').split('_')[-1].split('-')
    a = float(numbers[0])
    b = float(numbers[1])
    co_set = co_combine(co_lists[a], co_lists[b])
    co_a = np.array([float(item.split('/')[0]) for item in co_set])
    co_b = np.array([float(item.split('/')[-1]) for item in co_set])
    points = np.column_stack((co_a, co_b))
    E_eV_array = np.array([i.get_total_energy() for i in source_read])
    return {'points':points, 'energy':E_eV_array}


co_list_1 = np.round(np.arange(-0.9, 0.9, 0.1).tolist(),2)
co_list_2 = np.round(np.arange(-0.5, 0.8, 0.1).tolist(),2)
co_list_3 = np.round(np.arange(-0.6, 0.8, 0.1).tolist(),2)
co_lists = {1: co_list_1,
            2: co_list_2,
            3: co_list_3}


# xyzfile ='./data/wm_train.xyz'
def energy_xyz(xyzfile):
    IN = read(xyzfile ,index = ':')
    #extracting the energies from wm_.xyz input file, E_input in eV
    E_input_eV=[i.get_total_energy() for i in IN]
    #extracting the positions from wm_.xyz input file, R_ang in ang
    R_ang=[np.reshape(np.array(i.get_positions()),9) for i in IN]
    #convert positions from wm_.xyz to have mass**(1/2) * bohr unit
    R_ang_amu_vec = [np.concatenate((r[:3]*(1/A)*O, r[3:9]*(1/A)*H)) for r in R_ang]
    #displacement vector, R-R_0, in amu unit
    dis_vec = [r -R0_amu_vec for r in R_ang_amu_vec]
    #get the coefficients for each vibration modes for each coordinates
    c_1 = np.array([np.dot(vec, Q_1_vec) for vec in dis_vec])
    c_2 = np.array([np.dot(vec, Q_2_vec) for vec in dis_vec])
    c_3 = np.array([np.dot(vec, Q_3_vec) for vec in dis_vec])

    int_func_1 = int_fn('./data/wm_1_0.1.xyz')
    int_func_2 = int_fn('./data/wm_2_0.1.xyz')
    int_func_3 = int_fn('./data/wm_3_0.1.xyz')

    e_1 = [int_func_1(c) for c in c_1]
    e_2 = [int_func_2(c) for c in c_2]
    e_3 = [int_func_3(c) for c in c_3]
    e_sum = [sum(values) for values in zip(e_1, e_2, e_3)]

    co_12 = np.column_stack((c_1,c_2))
    co_13 = np.column_stack((c_1,c_3))
    co_23 = np.column_stack((c_2,c_3))

    E_12 = griddata(int_2B('1-2')['points'], int_2B('1-2')['energy'], 
                    co_12 , method='cubic')
    E_13 = griddata(int_2B('1-3')['points'], int_2B('1-3')['energy'], 
                    co_13 , method='cubic')
    E_23 = griddata(int_2B('2-3')['points'], int_2B('2-3')['energy'], 
                    co_23 , method='cubic')
    return {'E_1':e_1, 'E_2':e_2, 'E_3':e_3, 'co_12':co_12, 'c_1':c_1, 'c_2':c_2,
            'E_12':E_12,'E_13':E_13,'E_23':E_23,
            'E_SUM':e_sum, 'E_input':E_input_eV, 'R_ang':R_ang}
    
def delta_write(xyzfile, boundary, E_input_eV, R_ang, *args):
    E_sum = [sum(values) for values in zip(*args)]
    Delta_ML = ['{:.10f}'.format(float(E_input_eV[i]-E_sum[i])) for i in range(len(E_input_eV))]
    base_filename = os.path.basename(xyzfile)
    tale = base_filename.split("_", 1)[-1]
    file_log =f"data/{boundary}/nor_13_{tale}"
    if not os.path.exists(f"data/{boundary}"):
        os.makedirs(f'data/{boundary}')
    with open(file_log, 'w')as f:
        for i in range(len(E_input_eV)):
            if E_input_eV[i] < boundary:
                f.write('3\n')
                f.write('Lattice="15.0 0.0 0.0 0.0 15.0 0.0 0.0 0.0 15.0" Properties=species:S:1:pos:R:3 energy={} pbc="T T T"\n'.format(Delta_ML[i]))
                f.write('{}\n'.format(R_format(R_ang[i].reshape((3,3)))))
    print('Results written to {}'.format(file_log))

boundary = 50
file = './data/dir_train.xyz'
energy_data = energy_xyz(file)
delta_write(file, boundary,
            energy_data['E_input'],
            energy_data['R_ang'],
            energy_data['E_13'],
            energy_data['E_2'])
            # energy_data['E_13'],
            # energy_data['E_23'],
            # [-x for x in energy_data['E_SUM']],)

file = './data/dir_test.xyz'
energy_data = energy_xyz(file)
delta_write(file, boundary,
            energy_data['E_input'],
            energy_data['R_ang'],
            energy_data['E_13'],
            energy_data['E_2'])
            # energy_data['E_13'],
            # energy_data['E_23'],
            # [-x for x in energy_data['E_SUM']],)

file = './data/dir_valid.xyz'
energy_data = energy_xyz(file)
delta_write(file, boundary,
            energy_data['E_input'],
            energy_data['R_ang'],
            energy_data['E_13'],
            energy_data['E_2'])
            # energy_data['E_13'],
            # energy_data['E_23'],
            # [-x for x in energy_data['E_SUM']],)
