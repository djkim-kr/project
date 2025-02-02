import numpy as np
from ase.io import read

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240

def quip(word):
    energies=[]
    for i in range(1,6):
        quip_file=f'./data/quip_test/{word}/quip_test_{i}.xyz'
        e_temp=[i.get_total_energy() for i in read(quip_file ,index = ':')]
        energies.append(e_temp)
    e_test=[i.get_total_energy() for i in read(f'./data/{word}_test.xyz', index=':')]
    energies=np.array(energies)
    e_quip=np.mean(energies, axis=0)
    ml_delta = [1000*x for x in e_test]
    ml_error = [1000*(a-b) for a,b in zip(e_test,e_quip)]
    ml_pot=[x*1000 for x in e_quip]
    return ml_delta, ml_error, ml_pot

distance_list = np.linspace(0.70, 1.40, 20)
angle_list = np.linspace(60, 150, 30)

combined_list = []
for angle_value in angle_list:
    for distance_value1 in distance_list:
        for distance_value2 in distance_list:
            combined_list.append([distance_value1, distance_value2, angle_value])
dis_1 = [x[0] for x in combined_list]
dis_2 = [x[1] for x in combined_list]
ang_1 = [x[2] for x in combined_list]

delta_har, ml_error_har, ml_pot_har = quip('har')
delta_int, ml_error_int, ml_pot_int = quip('int')
delta_nor, ml_error_nor, ml_pot_nor = quip('nor')
delta_dir, ml_error_dir, ml_pot_dir = quip('har')

data = np.column_stack((dis_1, dis_2, ang_1,
                        delta_dir, delta_har, delta_nor, delta_int,
                        ml_error_dir, ml_error_har, ml_error_nor, ml_error_int,
                        ml_pot_dir, ml_pot_har, ml_pot_nor, ml_pot_int))

header = 'dis_1 dis_2 ang_1 delta_dir, delta_har, delta_nor delta_int ml_error_dir ml_error_har ml_error_nor ml_error_int ml_pot_dir ml_pot_har ml_pot_nor ml_pot_int'
np.savetxt(f'made_data.txt', data, delimiter=' ', header=header, comments='', fmt='%s')