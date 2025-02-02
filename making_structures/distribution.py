import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# 임의의 데이터 생성
# data = np.random.normal(0, 1, 1000)

# # 히스토그램 그리기
# plt.hist(data, bins=30, alpha=0.5, color='blue', edgecolor='black')
# plt.title('Histogram of Data')
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.show()

# # 밀도 플롯 그리기
# sns.kdeplot(data, fill=True)
# plt.title('Density Plot of Data')
# plt.xlabel('Value')
# plt.show()



def get_list(zmat_file):
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
    distance = distance_1_list + distance_2_list
    return np.array(distance), np.array(angle_1_list)

def density_plt(list, xlabel, set):
    sns.kdeplot(list, label=f'{set}: {np.round(list.min(),2)}~{np.round(list.max(),2)}', fill=True, alpha=0.4)
    plt.xlabel(xlabel)
    
zmat_file = './wm_made/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'made')

zmat_file = './wm_made2/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'made_2')

zmat_file = '../internal mode/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'original')

zmat_file = './wm_made3_20x20x30/data/wm_combine.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'made_3')

zmat_file = './wm_made5,6_nor_regular_grid/data/wm_combine5.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'made_5')

zmat_file = './wm_made5,6_nor_regular_grid/data/wm_combine6.zmat'
distance, angle = get_list(zmat_file)
density_plt(distance, 'distance(Å)', 'made_6')


plt.title('Distribution of distance')
plt.legend()
plt.savefig('Distance distribution.png', dpi=300)
plt.show()


zmat_file = './wm_made/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'made')

zmat_file = './wm_made2/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'made_2')

zmat_file = '../internal mode/data/wm_train.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'original')

zmat_file = './wm_made3_20x20x30/data/wm_combine.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'made_3')

zmat_file = './wm_made5,6_nor_regular_grid/data/wm_combine5.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'made_5')

zmat_file = './wm_made5,6_nor_regular_grid/data/wm_combine6.zmat'
distance, angle = get_list(zmat_file)
density_plt(angle, 'Angle(°)', 'made_6')


plt.title('Distribution of Angle')
plt.legend()
plt.savefig('Angle distribution.png', dpi=300)
plt.show()

