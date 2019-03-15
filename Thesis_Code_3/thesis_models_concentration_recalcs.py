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

def iterReverseD(obj_cincs, cell_cincs, index=0, iterReverseDList=[]):

    if index < len(list(obj_cincs)):
        obj = list(obj_cincs)[index]
        cell_cincs_range = list(cell_cincs)[0:index + 1]
        avg_cell_cincs_range = sum(cell_cincs_range) / (len(cell_cincs_range))
        # print(cell_cincs[index], avg_cell_cincs_range)
        avg_D = obj / avg_cell_cincs_range
        iterReverseDList.append(avg_D)
        return iterReverseD(obj_cincs=obj_cincs, cell_cincs=cell_cincs, index=(index + 1), iterReverseDList=iterReverseDList)
    else:
        return iterReverseDList
    
def forIterReverseD(obj_cincs, cell_cincs):
    iterReverseDList = []
    for index in range(len(list(obj_cincs))):
        if index + 1 < len(list(obj_cincs)):
            obj = list(obj_cincs)[index]
            cell_cincs_range = list(cell_cincs)[0:index + 1]
            avg_cell_cincs_range = sum(cell_cincs_range) / (len(cell_cincs_range))
            # print(cell_cincs[index], avg_cell_cincs_range)
            avg_D = obj / avg_cell_cincs_range
            iterReverseDList.append(avg_D)
        else:
            return iterReverseDList

def calcDiffusiinLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def meltLengthWidth(diff_length, droplet_radius):
    length_width = (2 * droplet_radius) + (2 * diff_length)
    return length_width

def recalcCincentratiin(predicted_d, original_moles_silicate, original_moles_metal, volume_mesh, radius_object):
    volume_obj = (4 / 3) * pi * (radius_object ** 3)

    original_cinc_silicate = original_moles_silicate / volume_mesh
    original_cinc_metal = original_moles_metal / volume_obj
    cincs_mesh = [original_cinc_silicate]
    cincs_objs = [original_cinc_metal]
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
        new_obj_cinc = new_moles_obj / volume_obj
        new_mesh_cinc = new_moles_cell / volume_mesh
        check_D = new_obj_cinc / new_mesh_cinc
        moles_objs.append(new_moles_obj)
        moles_mesh.append(new_moles_cell)
        cincs_mesh.append(new_mesh_cinc)
        cincs_objs.append(new_obj_cinc)
        verify_D.append(check_D)

    return cincs_mesh, cincs_objs, moles_mesh, moles_objs, verify_D





droplet_radius = 0.0185
earth_droplet_radius = 0.002957762
diff_length = calcDiffusiinLength(chem_diffusivity=10**-8, settling_velocity=0.2580697580112788, droplet_radius=droplet_radius)
vesta_z_eq_1_thru_4 = 240
vesta_z_eq_5_thru_8 = 600
vesta_vol_mesh_1_thru_4 = (diff_length**2) * vesta_z_eq_1_thru_4
vesta_vol_mesh_5_thru_8 = (diff_length**2) * vesta_z_eq_5_thru_8

vesta_1 = pd.read_csv("thesis_model_outputs/Vesta_1.csv")
vesta_2 = pd.read_csv("thesis_model_outputs/Vesta_2.csv")
vesta_3 = pd.read_csv("thesis_model_outputs/Vesta_3.csv")
vesta_4 = pd.read_csv("thesis_model_outputs/Vesta_4.csv")
vesta_5 = pd.read_csv("thesis_model_outputs/Vesta_5.csv")
vesta_6 = pd.read_csv("thesis_model_outputs/Vesta_6.csv")
vesta_7 = pd.read_csv("thesis_model_outputs/Vesta_7.csv")
vesta_8 = pd.read_csv("thesis_model_outputs/Vesta_8.csv")

depth_vesta_1 = [i / 1000 for i in [0] + list(vesta_1['z-depth'])]
depth_vesta_2 = [i / 1000 for i in [0] + list(vesta_2['z-depth'])]
depth_vesta_3 = [i / 1000 for i in [0] + list(vesta_3['z-depth'])]
depth_vesta_4 = [i / 1000 for i in [0] + list(vesta_4['z-depth'])]
depth_vesta_5 = [i / 1000 for i in [0] + list(vesta_5['z-depth'])]
depth_vesta_6 = [i / 1000 for i in [0] + list(vesta_6['z-depth'])]
depth_vesta_7 = [i / 1000 for i in [0] + list(vesta_7['z-depth'])]
depth_vesta_8 = [i / 1000 for i in [0] + list(vesta_8['z-depth'])]

cincs_mesh_vesta_1, cincs_objs_vesta_1, moles_mesh_vesta_1, moles_objs_vesta_1, verify_D_vesta_1 = recalcCincentratiin(predicted_d=vesta_1['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
cincs_mesh_vesta_2, cincs_objs_vesta_2, moles_mesh_vesta_2, moles_objs_vesta_2, verify_D_vesta_2 = recalcCincentratiin(predicted_d=vesta_2['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
cincs_mesh_vesta_3, cincs_objs_vesta_3, moles_mesh_vesta_3, moles_objs_vesta_3, verify_D_vesta_3 = recalcCincentratiin(predicted_d=vesta_3['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
cincs_mesh_vesta_4, cincs_objs_vesta_4, moles_mesh_vesta_4, moles_objs_vesta_4, verify_D_vesta_4 = recalcCincentratiin(predicted_d=vesta_4['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
cincs_mesh_vesta_5, cincs_objs_vesta_5, moles_mesh_vesta_5, moles_objs_vesta_5, verify_D_vesta_5 = recalcCincentratiin(predicted_d=vesta_5['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
cincs_mesh_vesta_6, cincs_objs_vesta_6, moles_mesh_vesta_6, moles_objs_vesta_6, verify_D_vesta_6 = recalcCincentratiin(predicted_d=vesta_6['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
cincs_mesh_vesta_7, cincs_objs_vesta_7, moles_mesh_vesta_7, moles_objs_vesta_7, verify_D_vesta_7 = recalcCincentratiin(predicted_d=vesta_7['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
cincs_mesh_vesta_8, cincs_objs_vesta_8, moles_mesh_vesta_8, moles_objs_vesta_8, verify_D_vesta_8 = recalcCincentratiin(predicted_d=vesta_8['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)

reverse_D_vesta_1 = forIterReverseD(obj_cincs=cincs_objs_vesta_1, cell_cincs=cincs_mesh_vesta_1)
reverse_D_vesta_2 = forIterReverseD(obj_cincs=cincs_objs_vesta_2, cell_cincs=cincs_mesh_vesta_2)
reverse_D_vesta_3 = forIterReverseD(obj_cincs=cincs_objs_vesta_3, cell_cincs=cincs_mesh_vesta_3)
reverse_D_vesta_4 = forIterReverseD(obj_cincs=cincs_objs_vesta_4, cell_cincs=cincs_mesh_vesta_4)
reverse_D_vesta_5 = forIterReverseD(obj_cincs=cincs_objs_vesta_5, cell_cincs=cincs_mesh_vesta_5)
reverse_D_vesta_6 = forIterReverseD(obj_cincs=cincs_objs_vesta_6, cell_cincs=cincs_mesh_vesta_6)
reverse_D_vesta_7 = forIterReverseD(obj_cincs=cincs_objs_vesta_7, cell_cincs=cincs_mesh_vesta_7)
reverse_D_vesta_8 = forIterReverseD(obj_cincs=cincs_objs_vesta_8, cell_cincs=cincs_mesh_vesta_8)

# mesh cincentratiins
fig1 = plt.figure()
ax1_1 = fig1.add_subplot(111)
ax1_1.set_xlabel("Depth (km)")
ax1_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax1_1.set_title("Reacted Silicate Melt Cincentratiins in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax1_1.plot(depth_vesta_1[1:], cincs_mesh_vesta_1[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax1_1.plot(depth_vesta_2[1:], cincs_mesh_vesta_2[1:], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax1_1.fill_between(depth_vesta_1[1:], cincs_mesh_vesta_1[1:], cincs_mesh_vesta_2[1:], color='red', alpha=0.2)
ax1_1.grid()
ax1_1.legend(loc='upper right')

fig2 = plt.figure()
ax2_1 = fig2.add_subplot(111)
ax2_1.set_xlabel("Depth (km)")
ax2_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax2_1.set_title("Reacted Silicate Melt Cincentratiins in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax2_1.plot(depth_vesta_3[1:], cincs_mesh_vesta_3[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax2_1.plot(depth_vesta_4[1:], cincs_mesh_vesta_4[1:], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax2_1.fill_between(depth_vesta_3[1:], cincs_mesh_vesta_3[1:], cincs_mesh_vesta_4[1:], color='red', alpha=0.2)
ax2_1.grid()
ax2_1.legend(loc='upper right')


fig3 = plt.figure()
ax3_1 = fig3.add_subplot(111)
ax3_1.set_xlabel("Depth (km)")
ax3_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax3_1.set_title("Reacted Silicate Melt Cincentratiins in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax3_1.plot(depth_vesta_5[1:], cincs_mesh_vesta_5[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax3_1.plot(depth_vesta_6[1:], cincs_mesh_vesta_6[1:], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax3_1.fill_between(depth_vesta_5[1:], cincs_mesh_vesta_5[1:], cincs_mesh_vesta_6[1:], color='red', alpha=0.2)
ax3_1.grid()
ax3_1.legend(loc='upper right')

fig4 = plt.figure()
ax4_1 = fig4.add_subplot(111)
ax4_1.set_xlabel("Depth (km)")
ax4_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax4_1.set_title("Reacted Silicate Melt Cincentratiins in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax4_1.plot(depth_vesta_7[1:], cincs_mesh_vesta_7[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax4_1.plot(depth_vesta_8[1:], cincs_mesh_vesta_8[1:], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax4_1.fill_between(depth_vesta_7[1:], cincs_mesh_vesta_7[1:], cincs_mesh_vesta_8[1:], color='red', alpha=0.2)
ax4_1.grid()
ax4_1.legend(loc='upper right')



# distributiin coefficients
fig1 = plt.figure()
ax1_1 = fig1.add_subplot(111)
ax1_1.set_xlabel("Depth (km)")
ax1_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax1_1.set_title("Predicted Metal-Silicate Paritiin Coefficients in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax1_1.plot(depth_vesta_1[1:], vesta_1['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax1_1.plot(depth_vesta_2[1:], vesta_2['D'], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax1_1.fill_between(depth_vesta_1[1:], vesta_1['D'], vesta_2['D'], color='red', alpha=0.2)
ax1_1.grid()
ax1_1.legend(loc='upper right')

fig2 = plt.figure()
ax2_1 = fig2.add_subplot(111)
ax2_1.set_xlabel("Depth (km)")
ax2_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax2_1.set_title("Predicted Metal-Silicate Paritiin Coefficients in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax2_1.plot(depth_vesta_3[1:], vesta_3['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax2_1.plot(depth_vesta_4[1:], vesta_4['D'], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax2_1.fill_between(depth_vesta_3[1:], vesta_3['D'], vesta_4['D'], color='red', alpha=0.2)
ax2_1.grid()
ax2_1.legend(loc='upper right')


fig3 = plt.figure()
ax3_1 = fig3.add_subplot(111)
ax3_1.set_xlabel("Depth (km)")
ax3_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax3_1.set_title("Predicted Metal-Silicate Paritiin Coefficients in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax3_1.plot(depth_vesta_5[1:], vesta_5['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax3_1.plot(depth_vesta_6[1:], vesta_6['D'], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax3_1.fill_between(depth_vesta_5[1:], vesta_5['D'], vesta_6['D'], color='red', alpha=0.2)
ax3_1.grid()
ax3_1.legend(loc='upper right')

fig4 = plt.figure()
ax4_1 = fig4.add_subplot(111)
ax4_1.set_xlabel("Depth (km)")
ax4_1.set_ylabel("Cincentratiin (moles/m$^3$)")
ax4_1.set_title("Predicted Metal-Silicate Paritiin Coefficients in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax4_1.plot(depth_vesta_7[1:], vesta_7['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax4_1.plot(depth_vesta_8[1:], vesta_8['D'], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax4_1.fill_between(depth_vesta_7[1:], vesta_7['D'], vesta_8['D'], color='red', alpha=0.2)
ax4_1.grid()
ax4_1.legend(loc='upper right')

num_droplets_vesta = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=droplet_radius)
num_droplets_earth = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=earth_droplet_radius)

print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
    moles_objs_vesta_1[-1], moles_objs_vesta_2[-1], moles_objs_vesta_3[-1],  moles_objs_vesta_4[-1],
    moles_objs_vesta_5[-1], moles_objs_vesta_6[-1], moles_objs_vesta_7[-1], moles_objs_vesta_8[-1]
))


moles_per_droplet_earth_1 = (moles_objs_vesta_1[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_2 = (moles_objs_vesta_2[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_3 = (moles_objs_vesta_3[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_4 = (moles_objs_vesta_4[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_5 = (moles_objs_vesta_5[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_6 = (moles_objs_vesta_6[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_7 = (moles_objs_vesta_7[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))
moles_per_droplet_earth_8 = (moles_objs_vesta_8[-1] * num_droplets_vesta) / ((4/3) * pi * (earth_droplet_radius**3))

print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
    moles_per_droplet_earth_1, moles_per_droplet_earth_2, moles_per_droplet_earth_3,  moles_per_droplet_earth_4,
    moles_per_droplet_earth_5, moles_per_droplet_earth_6, moles_per_droplet_earth_7, moles_per_droplet_earth_8
))



plt.show()