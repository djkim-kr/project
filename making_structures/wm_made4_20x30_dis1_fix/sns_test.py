import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from ase.io import read
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.cm as cm
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

#origin point
# O
# H      1      0.95883
# H      1      0.95883    2    104.34240
eq_dis = 0.95883
eq_ang = 104.34240

zmat_1='./data/wm_original.zmat'
org_dis_1, org_dis_2, org_ang = internal_coordinates(zmat_1)
org_dis_2 = [org_dis_2[i] for i in range(len(org_dis_2)) 
             if eq_dis -0.1 < org_dis_1[i] < eq_dis + 0.1]
org_ang = [org_ang[i] for i in range(len(org_ang)) 
             if eq_dis -0.1 < org_dis_1[i] < eq_dis + 0.1]

zmat_2='./data/wm_made.zmat'
made_dis_1, made_dis_2, made_ang = internal_coordinates(zmat_2)
made_dis_2 = [made_dis_2[i] for i in range(len(made_dis_2)) 
             if eq_dis -0.1 < made_dis_1[i] < eq_dis + 0.1]
made_ang = [made_ang[i] for i in range(len(made_ang)) 
             if eq_dis -0.1 < made_dis_1[i] < eq_dis + 0.1]

# # 플롯 생성
# fig = plt.figure(figsize=(10, 8))
# grid = plt.GridSpec(4, 4,) 
#                     # hspace=0.4, wspace=0.4)
# # 축 설정
# main_ax = fig.add_subplot(grid[1:4, 0:3])
# y_hist = fig.add_subplot(grid[1:4, 3], sharey=main_ax)
# x_hist = fig.add_subplot(grid[0, 0:3], sharex=main_ax)

# # 메인 플롯
# sns.scatterplot(x=made_dis_2, y=made_ang, marker="+", s=10, color='orange', label='made', ax=main_ax)
# # sns.kdeplot(x=made_dis_2, y=made_ang, cmap="Oranges", fill=True, thresh=0, alpha=0.5, ax=main_ax)

# sns.scatterplot(x=org_dis_2, y=org_ang, marker="+", s=10, color='blue', label='original', ax=main_ax)
# # sns.kdeplot(x=org_dis_2, y=org_ang, cmap="Blues", fill=True, thresh=0, alpha=0.5, ax=main_ax)

# # x축 히스토그램
# sns.histplot(made_dis_2, kde=True, color='orange', ax=x_hist, alpha=0.5)
# sns.histplot(org_dis_2, kde=True, color='blue', ax=x_hist, alpha=0.5)
# x_hist.axis('off')

# # y축 히스토그램
# sns.histplot(y=made_ang, kde=True, color='orange', ax=y_hist, alpha=0.5 )
# sns.histplot(y=org_ang, kde=True, color='blue', ax=y_hist, alpha=0.5 )
# y_hist.axis('off')

# # 레전드
# main_ax.legend()
# # plt.show()

words = ['nor', 'int','har','dir']
for word in words:
    # word='nor'
    file = f'./data/{word}_combine.xyz'

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
    z= delta_energy

    xi, yi = np.linspace(dis.min(), dis.max(), 100), np.linspace(ang.min(), ang.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    Z = griddata(np.c_[dis, ang], z, (xi, yi), method='cubic')
    fig = plt.figure(figsize=(12,8))
    fig.set_facecolor('white')
    grid = plt.GridSpec(4, 5) 
    # 축 설정
    ax = fig.add_subplot(grid[1:4, 0:3])
    x_hist = fig.add_subplot(grid[0, 0:3], sharex=ax)
    y_hist = fig.add_subplot(grid[1:4, 3], sharey=ax)
    colorbar = fig.add_subplot(grid[1:4, 4])
    colorbar.axis('off')

    contour1=ax.contour(xi, yi, Z, levels=10, colors='k', linewidths=1, linestyles='--')
    contour2 = ax.contourf(xi, yi, Z, levels=256, cmap='jet')
                            # , vmin=-2.0, vmax=2.0)

    ax.clabel(contour1, contour1.levels, inline=True)
    # ax.scatter(dis, ang, color='k', alpha=0.5, s=5)

    ## 컬러바 생성
    # norm = mpl.colors.Normalize(vmin=-2.0, vmax=2.0)
    # colormapping = cm.ScalarMappable(norm=norm, cmap='jet')
    # fig.colorbar( colormapping, ax=colorbar, fraction=0.05, pad=0.1, location = 'left',label='E(eV)')

    fig.colorbar( contour2, ax=colorbar, location = 'left', shrink=0.5, label='E(eV)')

    ax.axhline(y=eq_ang, color='red', linestyle='--', alpha=0.5)
    ax.axvline(x=eq_dis, color='red', linestyle='--', alpha=0.5)
    # 메인 플롯
    sns.scatterplot(x=made_dis_2, y=made_ang, s=5, color='green', label='made', ax=ax, alpha=0.7)
    sns.kdeplot(x=made_dis_2, y=made_ang, cmap="Greens", fill=False, thresh=0, ax=ax)
    sns.scatterplot(x=org_dis_2, y=org_ang, s=5, color='blue', label='original', ax=ax, alpha=0.7)
    sns.kdeplot(x=org_dis_2, y=org_ang, cmap="Blues", fill=False, thresh=0, ax=ax)

    # x축 히스토그램
    sns.histplot(made_dis_2, kde=True, color='green', ax=x_hist, alpha=0.5)
    sns.histplot(org_dis_2, kde=True, color='blue', ax=x_hist, alpha=0.5)
    x_hist.axis('off')

    # y축 히스토그램
    sns.histplot(y=made_ang, kde=True, color='green', ax=y_hist, alpha=0.5 )
    sns.histplot(y=org_ang, kde=True, color='blue', ax=y_hist, alpha=0.5 )
    y_hist.axis('off')

    # 레전드
    ax.legend()
    ax.set_xlim(0.70, 1.40)
    ax.set_ylim(60, 150)
    ax.set_title(f'{word}_delta_contour')
    ax.set_xlabel(f'distance(Å)')
    ax.set_ylabel(f'angle(°)')
    plt.savefig(f'./figure/{word}_distribution_contour.png',dpi=300)
    plt.tight_layout()
    plt.show()