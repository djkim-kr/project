import numpy as np
from ase.io import read
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240
eq_dis = 0.95883
eq_ang = 104.34240
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

# word='nor'
# file = f'./data/{word}_combine.xyz'

# # zmat_1='./data/wm_original.zmat'
# # org_dis_1, org_dis_2, org_ang = internal_coordinates(zmat_1)
# # zmat_2='./data/wm_made.zmat'
# # made_dis_1, made_dis_2, made_ang = internal_coordinates(zmat_1)

# distance_list = np.linspace(0.70, 1.40, 20)
# angle_list = np.linspace(60, 150, 30)
# combined_list = []

# for angle_value in angle_list:
#     for distance_value1 in distance_list:
#             combined_list.append([distance_value1, angle_value])

# dis = np.array([combined_list[i][0] for i in range(len(combined_list))])
# ang = np.array([combined_list[i][1] for i in range(len(combined_list))])
# source = read(file, index=':')
# delta_energy = np.array([i.get_total_energy() for i in source])
# z= delta_energy

# xi, yi = np.linspace(dis.min(), dis.max(), 100), np.linspace(ang.min(), ang.max(), 100)
# xi, yi = np.meshgrid(xi, yi)

# Z = griddata(np.c_[dis, ang], z, (xi, yi), method='cubic')
# fig = plt.figure(figsize=(12,8))
# fig.set_facecolor('white')
# ax=fig.add_subplot()
# contour1=ax.contour(xi, yi, Z, levels=10, colors='k', linewidths=1, linestyles='--')
# contour2 = ax.contourf(xi, yi, Z, levels=256, cmap='jet')

# ax.clabel(contour1, contour1.levels, inline=True)
# fig.colorbar(contour2, shrink=0.5, label='E(eV)')
# ax.scatter(dis, ang, color='k', alpha=0.5, s=5)
# plt.title(f'{word}_delta_contour')
# plt.xlabel(f'distance(Å)')
# plt.ylabel(f'angle(°)')
# # plt.savefig(f'./figure/{word}_delta_contour.png',dpi=300)
# plt.show()

def generate_contour_plot(word, ax, vmin, vmax):
    file = f'./data/{word}_combine.xyz'
    zmat_1 = './data/wm_original.zmat'
    org_dis_1, org_dis_2, org_ang = internal_coordinates(zmat_1)
    zmat_2 = './data/wm_made.zmat'
    made_dis_1, made_dis_2, made_ang = internal_coordinates(zmat_2)

    distance_list = np.linspace(0.70, 1.40, 20)
    angle_list = np.linspace(60, 150, 30)
    combined_list = []

    for angle_value in angle_list:
        for distance_value1 in distance_list:
            combined_list.append([distance_value1, angle_value])

    dis = np.array([combined_list[i][0] for i in range(len(combined_list))])
    ang = np.array([combined_list[i][1] for i in range(len(combined_list))])
    source = read(file, index=':')
    delta_energy = np.array([i.get_total_energy() for i in source])
    z = delta_energy

    xi, yi = np.linspace(dis.min(), dis.max(), 100), np.linspace(ang.min(), ang.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    Z = griddata(np.c_[dis, ang], z, (xi, yi), method='cubic')
    contour1 = ax.contour(xi, yi, Z, levels=50, colors='k', linewidths=1, linestyles='--')
    contour2 = ax.contourf(xi, yi, Z, levels=256, cmap='jet'
                           , vmin=vmin, vmax=vmax)
    
    ax.scatterplot(x=made_dis_2, y=made_ang, s=5, color='white', label='made', alpha=0.5)

    ax.axhline(y=eq_ang, color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=eq_dis, color='red', linestyle='--', alpha=0.5)

    ax.clabel(contour1, contour1.levels, inline=True)
    # fig.colorbar(contour2, shrink=0.5, label='E(eV)')
    ax.scatter(dis, ang, color='k', alpha=0.5, s=5)
    # eq_dis = 0.95883
    # eq_ang = 104.34240  
    # ax.set_xlim(eq_dis - 0.2, eq_dis + 0.2)
    # ax.set_ylim(eq_ang - 20, eq_ang + 20)
    ax.set_title(f'{word}_delta_contour')
    ax.set_xlabel('distance(Å)')
    ax.set_ylabel('angle(°)')

# Create the figure and axes
fig, axs = plt.subplots(2, 2, figsize=(18, 12), constrained_layout=True)
fig.set_facecolor('white')

# Generate the plots
words = ['dir', 'nor', 'int', 'har']
vmin= -0.1
vmax= 0.1
positions = [(0, 0), (0, 1), (1, 0), (1, 1)]
for word, pos in zip(words, positions):
    ax = axs[pos[0], pos[1]]
    generate_contour_plot(word, ax, vmin, vmax)

## 컬러바 생성
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
colormapping = cm.ScalarMappable(norm=norm, cmap='jet')
fig.colorbar( colormapping, ax=axs, fraction=0.05, pad=0.1, label='E(eV)')

# plt.savefig(f'./figure/100meV_contour_more.png',dpi=300)
plt.show()