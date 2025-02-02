import numpy as np
import os
import subprocess

distance_list = np.linspace(0.70, 1.40, 20)
angle_list = np.linspace(60, 150, 30)

#len(combined_list)=12000 = 20x20x30
combined_list = []

for angle_value in angle_list:
    for distance_value1 in distance_list:
            combined_list.append([distance_value1, angle_value])

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240

output_folder = './data/zmat'
made_3_xyz = './data/made_4.xyz'

os.makedirs(output_folder, exist_ok=True)
for i in range(len(combined_list)):
    new_list = []
    new_list.append("O")
    new_list.append(f"H      1      0.95883")
    new_list.append(f"H      1      {combined_list[i][0]:.5f}    2    {combined_list[i][1]:.5f}")
    formatted_output = "\n".join(new_list)
    output_file_path = os.path.join(output_folder, f'made_{i+1}.zmat')
    with open(output_file_path, 'w') as output_file:
        output_file.write(formatted_output)
    command = f'python3 gc.py -zmat {output_file_path} >> {made_3_xyz}'
    subprocess.run(command, shell=True, cwd=os.getcwd())
     
    
# print('Conversion complete.')