import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, pi, sqrt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal



# velocity: 0.2580697580112788
# total num droplets: 2.2788731170907947e+20

def calcNumTotalDroplets(core_radius, droplet_radius):
    core_volume = (4/3) * pi * (core_radius**3)
    droplet_volume = (4/3) * pi * (droplet_radius**3)
    num_droplets = core_volume / droplet_volume
    return num_droplets

def iterReverseD(obj_concs, cell_concs, index=0, iterReverseDList=[]):

    if index < len(list(obj_concs)):
        obj = list(obj_concs)[index]
        cell_concs_range = list(cell_concs)[0:index + 1]
        avg_cell_concs_range = sum(cell_concs_range) / (len(cell_concs_range))
        # print(cell_concs[index], avg_cell_concs_range)
        avg_D = obj / avg_cell_concs_range
        iterReverseDList.append(avg_D)
        return iterReverseD(obj_concs=obj_concs, cell_concs=cell_concs, index=(index + 1), iterReverseDList=iterReverseDList)
    else:
        return iterReverseDList

def calcDiffusionLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def meltLengthWidth(diff_length, droplet_radius):
    length_width = (2 * droplet_radius) + (2 * diff_length)
    return length_width

def recalcConcentration(predicted_d, original_moles_silicate, original_moles_metal, volume_mesh, radius_object):
    volume_obj = (4 / 3) * pi * (radius_object ** 3)

    original_conc_silicate = original_moles_silicate / volume_mesh
    original_conc_metal = original_moles_metal / volume_obj
    concs_mesh = [original_conc_silicate]
    concs_objs = [original_conc_metal]
    moles_mesh = [original_moles_silicate]
    moles_objs = [original_moles_metal]
    verify_D = []
    for index, d in enumerate(list(predicted_d)):
        old_moles_obj = moles_objs[index - 1]
        old_moles_cell = moles_mesh[index - 1]

        adj_matrix = (old_moles_cell /
                      (1 + (3 * volume_mesh * (
                                  (4 * pi * (radius_object ** 3) * d) ** (-1)))))
        adj_object = (old_moles_obj /
                      (1 + (4 * pi * (radius_object ** 3) * d) * (
                                  (3 ** (-1)) * (volume_mesh**(-1)))))
        adj_moles = adj_matrix - adj_object

        # adjust the moles of the element in the object and matrix, respectively
        new_moles_obj = old_moles_obj + adj_moles
        new_moles_cell = old_moles_cell - adj_moles
        new_obj_conc = new_moles_obj / volume_obj
        new_mesh_conc = new_moles_cell / volume_mesh
        check_D = new_obj_conc / new_mesh_conc
        moles_objs.append(new_moles_obj)
        moles_mesh.append(new_moles_cell)
        concs_mesh.append(new_mesh_conc)
        concs_objs.append(new_obj_conc)
        verify_D.append(check_D)

    return concs_mesh[1:], concs_objs[1:], moles_mesh[1:], moles_objs[1:], verify_D[1:]

df_fO2_08 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-0.8.csv")
df_fO2_110 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-1.10.csv")
df_fO2_225 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-2.25.csv")
df_fO2_245 = pd.read_csv("avg_eucrite_models/eucrites_avg_D_coeffs_fO2_-2.45.csv")

D_fO2_08 = list(df_fO2_08['D'])
depth_fO2_08 = list(df_fO2_08['z-depth'] / 1000)
cell_temperatures_fO2_08 = df_fO2_08['cell_temperature']
cell_pressures_fO2_08 = df_fO2_08['cell_pressure']
obj_conc_fO2_08 = df_fO2_08['object_conc']
cell_conc_fO2_08 = df_fO2_08['cell_conc']
cell_moles_fO2_08 = list(df_fO2_08['cell_moles'])
object_moles_fO2_08 = list(df_fO2_08['object_moles'])
reverse_D_fO2_08 = iterReverseD(obj_concs=obj_conc_fO2_08, cell_concs=cell_conc_fO2_08, index=0, iterReverseDList=[])

D_fO2_110 = list(df_fO2_110['D'])
depth_fO2_110 = list(df_fO2_110['z-depth'] / 1000)
cell_temperatures_fO2_110 = df_fO2_110['cell_temperature']
cell_pressures_fO2_110 = df_fO2_110['cell_pressure']
obj_conc_fO2_110 = df_fO2_110['object_conc']
cell_conc_fO2_110 = df_fO2_110['cell_conc']
cell_moles_fO2_110 = list(df_fO2_110['cell_moles'])
object_moles_fO2_110 = list(df_fO2_110['object_moles'])
reverse_D_fO2_110 = iterReverseD(obj_concs=obj_conc_fO2_110, cell_concs=cell_conc_fO2_110, index=0, iterReverseDList=[])

D_fO2_225 = list(df_fO2_225['D'])
depth_fO2_225 = list(df_fO2_225['z-depth'] / 1000)
cell_temperatures_fO2_225 = df_fO2_225['cell_temperature']
cell_pressures_fO2_225 = df_fO2_225['cell_pressure']
obj_conc_fO2_225 = df_fO2_225['object_conc']
cell_conc_fO2_225 = df_fO2_225['cell_conc']
cell_moles_fO2_225 = list(df_fO2_225['cell_moles'])
object_moles_fO2_225 = list(df_fO2_225['object_moles'])
reverse_D_fO2_225 = iterReverseD(obj_concs=obj_conc_fO2_225, cell_concs=cell_conc_fO2_225, index=0, iterReverseDList=[])

D_fO2_245 = list(df_fO2_245['D'])
depth_fO2_245 = list(df_fO2_245['z-depth'] / 1000)
cell_temperatures_fO2_245 = df_fO2_245['cell_temperature']
cell_pressures_fO2_245 = df_fO2_245['cell_pressure']
obj_conc_fO2_245 = df_fO2_245['object_conc']
cell_conc_fO2_245 = df_fO2_245['cell_conc']
cell_moles_fO2_245 = list(df_fO2_245['cell_moles'])
object_moles_fO2_245 = list(df_fO2_245['object_moles'])
reverse_D_fO2_245 = iterReverseD(obj_concs=obj_conc_fO2_245, cell_concs=cell_conc_fO2_245, index=0, iterReverseDList=[])

num_dropelts_vesta = calcNumTotalDroplets(core_radius=113 * 1000, droplet_radius=0.0185)




original_moles_silicate = 0
original_moles_metal = 0
droplet_radius = 0.0185

diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.2580697580112788, droplet_radius=droplet_radius)
length = meltLengthWidth(droplet_radius=droplet_radius, diff_length=diff_length)
width = length
volume_mesh = 400 * length * width

concs_mesh_fO2_08, concs_objs_fO2_08, moles_mesh_fO2_08, moles_objs_fO2_08, verify_D_fO2_08 = recalcConcentration(predicted_d=D_fO2_08,
                              original_moles_silicate=0.034937625, original_moles_metal=0, volume_mesh=volume_mesh, radius_object=0.0185)
concs_mesh_fO2_110, concs_objs_fO2_110, moles_mesh_fO2_110, moles_objs_fO2_110, verify_D_fO2_110 = recalcConcentration(predicted_d=D_fO2_110,
                              original_moles_silicate=0.034937625, original_moles_metal=0, volume_mesh=volume_mesh, radius_object=0.0185)
concs_mesh_fO2_225, concs_objs_fO2_225, moles_mesh_fO2_225, moles_objs_fO2_225, verify_D_fO2_225 = recalcConcentration(predicted_d=D_fO2_225,
                              original_moles_silicate=0.034937625, original_moles_metal=0, volume_mesh=volume_mesh, radius_object=0.0185)
concs_mesh_fO2_245, concs_objs_fO2_245, moles_mesh_fO2_245, moles_objs_fO2_245, verify_D_fO2_245 = recalcConcentration(predicted_d=D_fO2_245,
                              original_moles_silicate=0.034937625, original_moles_metal=0, volume_mesh=volume_mesh, radius_object=0.0185)

new_reverse_D_fO2_08 = iterReverseD(obj_concs=concs_objs_fO2_08, cell_concs=concs_mesh_fO2_08, index=0, iterReverseDList=[])
new_reverse_D_fO2_110 = iterReverseD(obj_concs=concs_objs_fO2_110, cell_concs=concs_mesh_fO2_110, index=0, iterReverseDList=[])
new_reverse_D_fO2_225 = iterReverseD(obj_concs=concs_objs_fO2_225, cell_concs=concs_mesh_fO2_225, index=0, iterReverseDList=[])
new_reverse_D_fO2_245 = iterReverseD(obj_concs=concs_objs_fO2_245, cell_concs=concs_mesh_fO2_245, index=0, iterReverseDList=[])

new_df = pd.DataFrame({'depth': depth_fO2_08, 'old_obj_conc': obj_conc_fO2_08, 'new_obj_conc': concs_objs_fO2_08, 'old_mesh_conc': cell_conc_fO2_08, 'new_mesh_conc': concs_mesh_fO2_08})
new_df.to_csv("08_check.csv")


fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.plot(depth_fO2_08, D_fO2_08, label='fO2 = IW-0.8', color='red')
ax1.plot(depth_fO2_08, new_reverse_D_fO2_08, label='Reversed fO2 = IW-0.8', color='red', linestyle="--")
ax1.plot(depth_fO2_110, D_fO2_110, label='fO2 = IW-1.10', color='blue')
ax1.plot(depth_fO2_110, new_reverse_D_fO2_110, label='Reversed fO2 = IW-1.10', color='blue', linestyle="--")
ax1_2 = fig1.add_subplot(212)
ax1_2.plot(depth_fO2_08, [j - i for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)],
           label="Reversed - Distributed (fO2 = IW-0.8)", linestyle="--", color='red')
ax1_2.plot(depth_fO2_110, [j - i for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)],
           label="Reversed - Distributed (fO2 = IW-1.10)", linestyle="--", color='blue')
ax1.grid()
ax1.fill_between(depth_fO2_08, D_fO2_08, D_fO2_110, color='red', alpha=0.2)
ax1_2.fill_between(depth_fO2_225, [j - i for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)],
                   [j - i for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)], color='red', alpha=0.2)
ax1.set_title("Iron Meteorite fO2 Model")
ax1.legend(loc='center right')
ax1_2.legend(loc='upper left')
ax1_2.set_xlabel("Depth (km)")
ax1.set_ylabel("Bulk D ('Reversed')")
ax1_2.set_ylabel("Reversed - Distributed")
ax1_2.grid()


fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
ax2.plot(depth_fO2_225, D_fO2_225, label='fO2 = IW-2.25', color='red')
ax2.plot(depth_fO2_225, new_reverse_D_fO2_225, label='Reversed fO2 = IW-2.25', color='red', linestyle="--")
ax2.plot(depth_fO2_245, D_fO2_245, label='fO2 = IW-2.45', color='blue')
ax2.plot(depth_fO2_245, new_reverse_D_fO2_245, label='Reversed fO2 = IW-2.45', color='blue', linestyle="--")
ax2_2 = fig2.add_subplot(212)
ax2_2.plot(depth_fO2_225, [j - i for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)],
           label="Reversed - Distributed (fO2 = IW-2.25)", linestyle="--", color='red')
ax2_2.plot(depth_fO2_245, [j - i for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)],
           label="Reversed - Distributed (fO2 = IW-2.45)", linestyle="--", color='blue')
ax2.grid()
ax2.fill_between(depth_fO2_225, D_fO2_225, D_fO2_245, color='red', alpha=0.2)
ax2_2.fill_between(depth_fO2_225, [j - i for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)],
                   [j - i for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)], color='red', alpha=0.2)
ax2.legend(loc='center right')
ax2_2.legend(loc='upper left')
ax2.set_title("Vesta fO2 Model")
ax2_2.set_xlabel("Depth (km)")
ax2.set_ylabel("Bulk D ('Reversed')")
ax2_2.set_ylabel("Reversed - Distributed")
ax2_2.grid()


fig3 = plt.figure()
ax3 = fig3.add_subplot(211)
ax3.plot(depth_fO2_08, moles_mesh_fO2_08, linewidth=2.0, label="Moles in Mesh Cells (fO2: IW-0.8)")
ax3.plot(depth_fO2_110, moles_mesh_fO2_110, linewidth=2.0, label="Moles in Mesh Cells (fO2: IW-1.10)")
ax3_2 = fig3.add_subplot(212)
ax3_2.plot(depth_fO2_225, moles_mesh_fO2_225, linewidth=2.0, label="Moles in Mesh Cells (fO2: IW-2.25)")
ax3_2.plot(depth_fO2_245, moles_mesh_fO2_245, linewidth=2.0, label="Moles in Mesh Cells (fO2: IW-2.45)")
ax3.grid()
ax3_2.grid()
ax3_2.set_xlabel("Depth (km)")
ax3.set_ylabel("Moles in Mesh Cell")
ax3_2.set_ylabel("Moles in Mesh Cell")
ax3.set_title("Mole Evolution in Mesh Cells Along Droplet Path of Descent")
ax3.fill_between(depth_fO2_08, moles_mesh_fO2_08, moles_mesh_fO2_110, color='red', alpha=0.2, label="Oxidizing Model")
ax3_2.fill_between(depth_fO2_08, moles_mesh_fO2_08, moles_mesh_fO2_110, color='blue', alpha=0.2, label="Oxidizing Model")
ax3.legend(loc='upper left')
ax3_2.legend(loc='upper left')

fig4 = plt.figure()
ax5 = fig4.add_subplot(211)
ax5.plot(depth_fO2_08, moles_objs_fO2_08, linewidth=2.0, label="Moles in Metal Droplets (fO2: IW-0.8)")
ax5.plot(depth_fO2_110, moles_objs_fO2_110, linewidth=2.0, label="Moles in Metal Droplets (fO2: IW-1.10)")
ax5_2 = fig4.add_subplot(212)
ax5_2.plot(depth_fO2_225, moles_objs_fO2_225, linewidth=2.0, label="Moles in Metal Droplets (fO2: IW-2.25)")
ax5_2.plot(depth_fO2_245, moles_objs_fO2_245, linewidth=2.0, label="Moles in Metal Droplets (fO2: IW-2.45)")
ax5.grid()
ax5_2.grid()
ax5_2.set_xlabel("Depth (km)")
ax5.set_ylabel("Moles in Metal Droplet")
ax5_2.set_ylabel("Moles in Metal Droplet")
ax5.set_title("Mole Evolution in Metal Droplet Along Droplet Path of Descent")
ax5.fill_between(depth_fO2_08, moles_objs_fO2_08, moles_objs_fO2_110, color='red', alpha=0.2, label="Reducing Model")
ax5_2.fill_between(depth_fO2_225, moles_objs_fO2_225, moles_objs_fO2_245, color='blue', alpha=0.2, label="Reducing Model")
ax5.legend(loc='upper left')
ax5_2.legend(loc='upper left')


fig5 = plt.figure()
ax5 = fig5.add_subplot(211)
ax5_2 = fig5.add_subplot(212)
ax5.plot(depth_fO2_08, [((j - i) / i) * 100 for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)])
ax5.plot(depth_fO2_110, [((j - i) / i) * 100 for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)])
ax5_2.plot(depth_fO2_225, [((j - i) / i) * 100 for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)])
ax5_2.plot(depth_fO2_245, [((j - i) / i) * 100 for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)])
ax5.grid()
ax5_2.grid()
ax5_2.set_xlabel("Depth (km)")
ax5.set_ylabel("Relative Change (%)")
ax5_2.set_ylabel("Relative Change (%)")
ax5.set_title("Relative Change in Averaged vs. Absolute D Over Depth (relative to absolute)")
ax5.fill_between(depth_fO2_08, [((j - i) / i) * 100 for i, j in zip(D_fO2_08, new_reverse_D_fO2_08)],
                 [((j - i) / i) * 100 for i, j in zip(D_fO2_110, new_reverse_D_fO2_110)], color='red', alpha=0.2, label="Reducing Model")
ax5_2.fill_between(depth_fO2_225, [((j - i) / i) * 100 for i, j in zip(D_fO2_225, new_reverse_D_fO2_225)],
                   [((j - i) / i) * 100 for i, j in zip(D_fO2_245, new_reverse_D_fO2_245)], color='blue', alpha=0.2, label="Reducing Model")
ax5.legend(loc='lower left')
ax5_2.legend(loc='lower left')

# a = [moles_objs_fO2_08[-1] * num_dropelts_vesta, moles_objs_fO2_110[-1] * num_dropelts_vesta, moles_objs_fO2_225[-1] * num_dropelts_vesta, moles_objs_fO2_245[-1] * num_dropelts_vesta]
# for i in a:
#     print('%.4E' % Decimal(str(i)))

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(depth_fO2_08, concs_mesh_fO2_08, linewidth=2.0, label=)

plt.show()