import numpy as np
import os
import subprocess

distance_1 = [0.70, 1.40]
# distance_2 = [0.75, 1.30]
angle_1 = [60, 150.0]

# Create 13 evenly spaced values between the minimum and maximum values
num_values = 15
distance_1_min = min(distance_1)
distance_1_max = max(distance_1)
distance_1_values = np.linspace(distance_1_min, distance_1_max, num_values)

angle_1_min = min(angle_1)
angle_1_max = max(angle_1)
angle_1_values = np.linspace(angle_1_min, angle_1_max, num_values)

distance_1 = []
# distance_2 = []
angle_1 = []

distance_1.extend(distance_1_values)
# distance_2.extend(distance_2_values)
angle_1.extend(angle_1_values)

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240

output_folder = './data/e_training_2/'
dis1_xyz = './data/e_training_2/dis_1.xyz'
# dis2_xyz = './data/e_training/dis_2.xyz'
ang1_xyz = './data/e_training_2/ang_1.xyz'

os.makedirs(output_folder, exist_ok=True)
for i, val in enumerate(distance_1):
    new_list = []
    new_list.append("O")
    new_list.append(f"H      1      {val:.5f}")
    new_list.append(f"H      1      0.95883    2    104.34240")
    formatted_output = "\n".join(new_list)
    output_file_path = os.path.join(output_folder, f'dis_1_{i+1}.zmat')
    with open(output_file_path, 'w') as output_file:
        output_file.write(formatted_output)
    command = f'python3 gc.py -zmat {output_file_path} >> {dis1_xyz}'
    subprocess.run(command, shell=True, cwd=os.getcwd())
    
for i, val in enumerate(angle_1):
    new_list = []
    new_list.append("O")
    new_list.append(f"H      1      0.95883")
    new_list.append(f"H      1      0.95883    2    {angle_1[i]:.5f}")
    formatted_output = "\n".join(new_list)
    output_file_path = os.path.join(output_folder, f'ang_1_{i+1}.zmat')
    with open(output_file_path, 'w') as output_file:
        output_file.write(formatted_output)
    command = f'python3 gc.py -zmat {output_file_path} >> {ang1_xyz}'
    subprocess.run(command, shell=True, cwd=os.getcwd())   
    
print('Conversion complete.')