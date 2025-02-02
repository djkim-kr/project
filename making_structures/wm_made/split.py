import os

input_file = str(input('enter xyz file: '))
base_name =os.path.splitext(input_file)[0]
output_files = [f'{base_name}_train.xyz', f'{base_name}_valid.xyz', f'{base_name}_test.xyz']
line_ranges = [(1,12000),(12001,13500),(13501,15000)]
with open(input_file, 'r') as file:
    lines = file.readlines()

for i, (start, end) in enumerate(line_ranges):
    output_file = output_files[i]
    with open (output_file, 'w') as output:
        for line in lines[start - 1:end]: 
            output.write(line)

print('File split into three pieces successfully.')