import os
import subprocess


word = 'final'  

output_folder = f'./data/wm_{word}'
os.makedirs(output_folder, exist_ok=True)

# Read the XYZ file
input_xyz_file = f'./data/wm_{word}.xyz'
with open(input_xyz_file, 'r') as input_file:
    lines = input_file.readlines()

# Split the lines into groups of 5
group_size = 5
for i in range(0, len(lines), group_size):
    group = lines[i:i + group_size]
    # Generate the output file name
    output_file_name = os.path.join(output_folder, f'wm_{word}_{i // group_size + 1}.xyz')
    # Write the group of lines to the output file
    with open(output_file_name, 'w') as output_file:
        output_file.writelines(group)
print(f'Split {len(lines)} lines into {len(os.listdir(output_folder))} files in the {output_folder} folder.')


final_zmat = f'./data/wm_{word}.zmat'
# Process the split XYZ files and append to the final ZMAT file
for i in range(len(os.listdir(output_folder))):
    number = i + 1
    filename = f'wm_{word}_{number}.xyz'
    xyz_file = os.path.join(output_folder, filename)
    command = f'python3 gc.py -xyz {xyz_file} >> {final_zmat}'
    # Execute the command in the terminal
    subprocess.run(command, shell=True, cwd=os.getcwd())

print('Conversion complete.')
