import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, pi, sqrt, log10
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal
import bisect

plt.rcParams.update({'font.size': 16})


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
    
def forIterReverseD(obj_concs, cell_concs, initial_moles_obj, initial_moles_mesh, volume_mesh, radius_droplet):
    volume_obj = (4/3) * pi * (radius_droplet**3)
    initial_conc_mesh = initial_moles_mesh / volume_mesh
    initial_conc_obj = initial_moles_obj / volume_obj
    initial_D = initial_conc_obj / initial_conc_mesh
    iterReverseDList = [initial_D]
    for index in range(len(list(obj_concs))):
        if index + 1 < len(list(obj_concs)):
            obj = list(obj_concs)[index]
            cell_concs_range = list(cell_concs)[0:index + 1]
            avg_cell_concs_range = sum(cell_concs_range) / (len(cell_concs_range))
            # print(cell_concs[index], avg_cell_concs_range)
            avg_D = obj / avg_cell_concs_range
            iterReverseDList.append(avg_D)
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
    verify_D = [concs_objs[0] / concs_mesh[0]]
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

    return concs_mesh, concs_objs, moles_mesh, moles_objs, verify_D



magma_ocean_concs_182w = pd.read_csv("182w_mantle.csv")
magma_ocean_concs_184w = pd.read_csv("184w_mantle.csv")


droplet_radius = 0.0185
earth_droplet_radius = 0.002957762
diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.2580697580112788, droplet_radius=droplet_radius)
earth_diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.646064514, droplet_radius=earth_droplet_radius)
vesta_z_eq_1_thru_4 = 240
vesta_z_eq_5_thru_8 = 600
vesta_vol_mesh_1_thru_4 = ((2 * (droplet_radius + diff_length))**2) * vesta_z_eq_1_thru_4
vesta_vol_mesh_5_thru_8 = ((2 * (droplet_radius + diff_length))**2) * vesta_z_eq_5_thru_8
earth_z_eq_1_thru_4 = 25
earth_z_eq_5_thru_8 = 65
earth_vol_mesh_1_thru_4 = ((2 * (earth_droplet_radius + earth_diff_length))**2) * earth_z_eq_1_thru_4
earth_vol_mesh_5_thru_8 = ((2 * (earth_droplet_radius + earth_diff_length))**2) * earth_z_eq_5_thru_8

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



concs_mesh_vesta_1, concs_objs_vesta_1, moles_mesh_vesta_1, moles_objs_vesta_1, verify_D_vesta_1 = recalcConcentration(predicted_d=vesta_1['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_2, concs_objs_vesta_2, moles_mesh_vesta_2, moles_objs_vesta_2, verify_D_vesta_2 = recalcConcentration(predicted_d=vesta_2['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_3, concs_objs_vesta_3, moles_mesh_vesta_3, moles_objs_vesta_3, verify_D_vesta_3 = recalcConcentration(predicted_d=vesta_3['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_4, concs_objs_vesta_4, moles_mesh_vesta_4, moles_objs_vesta_4, verify_D_vesta_4 = recalcConcentration(predicted_d=vesta_4['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_5, concs_objs_vesta_5, moles_mesh_vesta_5, moles_objs_vesta_5, verify_D_vesta_5 = recalcConcentration(predicted_d=vesta_5['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_6, concs_objs_vesta_6, moles_mesh_vesta_6, moles_objs_vesta_6, verify_D_vesta_6 = recalcConcentration(predicted_d=vesta_6['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_7, concs_objs_vesta_7, moles_mesh_vesta_7, moles_objs_vesta_7, verify_D_vesta_7 = recalcConcentration(predicted_d=vesta_7['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_8, concs_objs_vesta_8, moles_mesh_vesta_8, moles_objs_vesta_8, verify_D_vesta_8 = recalcConcentration(predicted_d=vesta_8['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)

# reverse_D_vesta_1 = forIterReverseD(obj_concs=concs_objs_vesta_1, cell_concs=concs_mesh_vesta_1)
# reverse_D_vesta_2 = forIterReverseD(obj_concs=concs_objs_vesta_2, cell_concs=concs_mesh_vesta_2)
# reverse_D_vesta_3 = forIterReverseD(obj_concs=concs_objs_vesta_3, cell_concs=concs_mesh_vesta_3)
# reverse_D_vesta_4 = forIterReverseD(obj_concs=concs_objs_vesta_4, cell_concs=concs_mesh_vesta_4)
# reverse_D_vesta_5 = forIterReverseD(obj_concs=concs_objs_vesta_5, cell_concs=concs_mesh_vesta_5)
# reverse_D_vesta_6 = forIterReverseD(obj_concs=concs_objs_vesta_6, cell_concs=concs_mesh_vesta_6)
# reverse_D_vesta_7 = forIterReverseD(obj_concs=concs_objs_vesta_7, cell_concs=concs_mesh_vesta_7)
# reverse_D_vesta_8 = forIterReverseD(obj_concs=concs_objs_vesta_8, cell_concs=concs_mesh_vesta_8)

# mesh concentrations
fig1 = plt.figure()
ax1_1 = fig1.add_subplot(111)
ax1_1.set_xlabel("Depth (km)")
ax1_1.set_ylabel("Concentration (moles/m$^3$)")
ax1_1.set_title("Reacted Silicate Melt Concentrations in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax1_1.plot(depth_vesta_1[1:], concs_mesh_vesta_1[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax1_1.plot(depth_vesta_2[1:], concs_mesh_vesta_2[1:], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax1_1.fill_between(depth_vesta_1[1:], concs_mesh_vesta_1[1:], concs_mesh_vesta_2[1:], color='red', alpha=0.2)
ax1_1.grid()
ax1_1.legend(loc='upper right')

fig2 = plt.figure()
ax2_1 = fig2.add_subplot(111)
ax2_1.set_xlabel("Depth (km)")
ax2_1.set_ylabel("Concentration (moles/m$^3$)")
ax2_1.set_title("Reacted Silicate Melt Concentrations in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax2_1.plot(depth_vesta_3[1:], concs_mesh_vesta_3[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax2_1.plot(depth_vesta_4[1:], concs_mesh_vesta_4[1:], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax2_1.fill_between(depth_vesta_3[1:], concs_mesh_vesta_3[1:], concs_mesh_vesta_4[1:], color='red', alpha=0.2)
ax2_1.grid()
ax2_1.legend(loc='upper right')


fig3 = plt.figure()
ax3_1 = fig3.add_subplot(111)
ax3_1.set_xlabel("Depth (km)")
ax3_1.set_ylabel("Concentration (moles/m$^3$)")
ax3_1.set_title("Reacted Silicate Melt Concentrations in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax3_1.plot(depth_vesta_5[1:], concs_mesh_vesta_5[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax3_1.plot(depth_vesta_6[1:], concs_mesh_vesta_6[1:], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax3_1.fill_between(depth_vesta_5[1:], concs_mesh_vesta_5[1:], concs_mesh_vesta_6[1:], color='red', alpha=0.2)
ax3_1.grid()
ax3_1.legend(loc='upper right')

fig4 = plt.figure()
ax4_1 = fig4.add_subplot(111)
ax4_1.set_xlabel("Depth (km)")
ax4_1.set_ylabel("Concentration (moles/m$^3$)")
ax4_1.set_title("Reacted Silicate Melt Concentrations in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax4_1.plot(depth_vesta_7[1:], concs_mesh_vesta_7[1:], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax4_1.plot(depth_vesta_8[1:], concs_mesh_vesta_8[1:], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax4_1.fill_between(depth_vesta_7[1:], concs_mesh_vesta_7[1:], concs_mesh_vesta_8[1:], color='red', alpha=0.2)
ax4_1.grid()
ax4_1.legend(loc='upper right')



# partition coefficients
fig9 = plt.figure()
ax9_1 = fig9.add_subplot(111)
ax9_1.set_xlabel("Depth (km)")
ax9_1.set_ylabel("Concentration (moles/m$^3$)")
ax9_1.set_title("Predicted Metal-Silicate Partition Coefficients in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax9_1.plot(depth_vesta_1[1:], vesta_1['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax9_1.plot(depth_vesta_2[1:], vesta_2['D'], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax9_1.fill_between(depth_vesta_1[1:], vesta_1['D'], vesta_2['D'], color='red', alpha=0.2)
ax9_1.grid()
ax9_1.legend(loc='upper right')

fig10 = plt.figure()
ax10_1 = fig10.add_subplot(111)
ax10_1.set_xlabel("Depth (km)")
ax10_1.set_ylabel("Concentration (moles/m$^3$)")
ax10_1.set_title("Predicted Metal-Silicate Partition Coefficients in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax10_1.plot(depth_vesta_3[1:], vesta_3['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax10_1.plot(depth_vesta_4[1:], vesta_4['D'], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax10_1.fill_between(depth_vesta_3[1:], vesta_3['D'], vesta_4['D'], color='red', alpha=0.2)
ax10_1.grid()
ax10_1.legend(loc='upper right')


fig11 = plt.figure()
ax11_1 = fig11.add_subplot(111)
ax11_1.set_xlabel("Depth (km)")
ax11_1.set_ylabel("Concentration (moles/m$^3$)")
ax11_1.set_title("Predicted Metal-Silicate Partition Coefficients in a Vestian Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax11_1.plot(depth_vesta_5[1:], vesta_5['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax11_1.plot(depth_vesta_6[1:], vesta_6['D'], linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax11_1.fill_between(depth_vesta_5[1:], vesta_5['D'], vesta_6['D'], color='red', alpha=0.2)
ax11_1.grid()
ax11_1.legend(loc='upper right')

fig12 = plt.figure()
ax12_1 = fig12.add_subplot(111)
ax12_1.set_xlabel("Depth (km)")
ax12_1.set_ylabel("Concentration (moles/m$^3$)")
ax12_1.set_title("Predicted Metal-Silicate Partition Coefficients in a Vestian Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax12_1.plot(depth_vesta_7[1:], vesta_7['D'], linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax12_1.plot(depth_vesta_8[1:], vesta_8['D'], linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax12_1.fill_between(depth_vesta_7[1:], vesta_7['D'], vesta_8['D'], color='red', alpha=0.2)
ax12_1.grid()
ax12_1.legend(loc='upper right')

num_droplets_vesta = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=droplet_radius)
num_droplets_earth = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=earth_droplet_radius)

print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
    moles_objs_vesta_1[-1], moles_objs_vesta_2[-1], moles_objs_vesta_3[-1],  moles_objs_vesta_4[-1],
    moles_objs_vesta_5[-1], moles_objs_vesta_6[-1], moles_objs_vesta_7[-1], moles_objs_vesta_8[-1]
))


moles_per_droplet_earth_1 = (moles_objs_vesta_1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_2 = (moles_objs_vesta_2[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_3 = (moles_objs_vesta_3[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_4 = (moles_objs_vesta_4[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_5 = (moles_objs_vesta_5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_6 = (moles_objs_vesta_6[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_7 = (moles_objs_vesta_7[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_8 = (moles_objs_vesta_8[-1] * num_droplets_vesta) / num_droplets_earth

print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
    moles_per_droplet_earth_1, moles_per_droplet_earth_2, moles_per_droplet_earth_3,  moles_per_droplet_earth_4,
    moles_per_droplet_earth_5, moles_per_droplet_earth_6, moles_per_droplet_earth_7, moles_per_droplet_earth_8
))

print("Num Droplets Vesta: {}\nNum Droplets Earth: {}".format(num_droplets_vesta, num_droplets_earth))


earth_1 = pd.read_csv("thesis_model_outputs/Earth_1.csv")
earth_2 = pd.read_csv("thesis_model_outputs/Earth_2.csv")
earth_3 = pd.read_csv("thesis_model_outputs/Earth_3.csv")
earth_4 = pd.read_csv("thesis_model_outputs/Earth_4.csv")
earth_5 = pd.read_csv("thesis_model_outputs/Earth_5.csv")
earth_6 = pd.read_csv("thesis_model_outputs/Earth_6.csv")
earth_7 = pd.read_csv("thesis_model_outputs/Earth_7.csv")
earth_8 = pd.read_csv("thesis_model_outputs/Earth_8.csv")

depth_earth_1 = [i / 1000 for i in [0] + list(earth_1['z-depth'])]
depth_earth_2 = [i / 1000 for i in [0] + list(earth_2['z-depth'])]
depth_earth_3 = [i / 1000 for i in [0] + list(earth_3['z-depth'])]
depth_earth_4 = [i / 1000 for i in [0] + list(earth_4['z-depth'])]
depth_earth_5 = [i / 1000 for i in [0] + list(earth_5['z-depth'])]
depth_earth_6 = [i / 1000 for i in [0] + list(earth_6['z-depth'])]
depth_earth_7 = [i / 1000 for i in [0] + list(earth_7['z-depth'])]
depth_earth_8 = [i / 1000 for i in [0] + list(earth_8['z-depth'])]

concs_mesh_earth_1, concs_objs_earth_1, moles_mesh_earth_1, moles_objs_earth_1, verify_D_earth_1 = recalcConcentration(predicted_d=earth_1['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_1, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
concs_mesh_earth_2, concs_objs_earth_2, moles_mesh_earth_2, moles_objs_earth_2, verify_D_earth_2 = recalcConcentration(predicted_d=earth_2['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_2, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
concs_mesh_earth_3, concs_objs_earth_3, moles_mesh_earth_3, moles_objs_earth_3, verify_D_earth_3 = recalcConcentration(predicted_d=earth_3['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_3, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
concs_mesh_earth_4, concs_objs_earth_4, moles_mesh_earth_4, moles_objs_earth_4, verify_D_earth_4 = recalcConcentration(predicted_d=earth_4['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_4, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
concs_mesh_earth_5, concs_objs_earth_5, moles_mesh_earth_5, moles_objs_earth_5, verify_D_earth_5 = recalcConcentration(predicted_d=earth_5['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_5, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
concs_mesh_earth_6, concs_objs_earth_6, moles_mesh_earth_6, moles_objs_earth_6, verify_D_earth_6 = recalcConcentration(predicted_d=earth_6['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_6, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
concs_mesh_earth_7, concs_objs_earth_7, moles_mesh_earth_7, moles_objs_earth_7, verify_D_earth_7 = recalcConcentration(predicted_d=earth_7['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_7, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
concs_mesh_earth_8, concs_objs_earth_8, moles_mesh_earth_8, moles_objs_earth_8, verify_D_earth_8 = recalcConcentration(predicted_d=earth_8['D'],
                              original_moles_silicate=0.27950089725326804, original_moles_metal=moles_per_droplet_earth_8, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)

# reverse_D_earth_1 = forIterReverseD(obj_concs=concs_objs_earth_1, cell_concs=concs_mesh_earth_1)
# reverse_D_earth_2 = forIterReverseD(obj_concs=concs_objs_earth_2, cell_concs=concs_mesh_earth_2)
# reverse_D_earth_3 = forIterReverseD(obj_concs=concs_objs_earth_3, cell_concs=concs_mesh_earth_3)
# reverse_D_earth_4 = forIterReverseD(obj_concs=concs_objs_earth_4, cell_concs=concs_mesh_earth_4)
# reverse_D_earth_5 = forIterReverseD(obj_concs=concs_objs_earth_5, cell_concs=concs_mesh_earth_5)
# reverse_D_earth_6 = forIterReverseD(obj_concs=concs_objs_earth_6, cell_concs=concs_mesh_earth_6)
# reverse_D_earth_7 = forIterReverseD(obj_concs=concs_objs_earth_7, cell_concs=concs_mesh_earth_7)
# reverse_D_earth_8 = forIterReverseD(obj_concs=concs_objs_earth_8, cell_concs=concs_mesh_earth_8)

# mesh concentrations
fig5 = plt.figure()
ax5_1 = fig5.add_subplot(111)
ax5_1.set_xlabel("Depth (km)")
ax5_1.set_ylabel("Concentration (moles/m$^3$)")
ax5_1.set_title("Reacted Silicate Melt Concentrations in an Earth Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax5_1.plot(depth_earth_1, concs_mesh_earth_1, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax5_1.plot(depth_earth_2, concs_mesh_earth_2, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax5_1.fill_between(depth_earth_1, concs_mesh_earth_1, concs_mesh_earth_2, color='red', alpha=0.2)
ax5_1.grid()
ax5_1.legend(loc='upper right')

fig6 = plt.figure()
ax6_1 = fig6.add_subplot(111)
ax6_1.set_xlabel("Depth (km)")
ax6_1.set_ylabel("Concentration (moles/m$^3$)")
ax6_1.set_title("Reacted Silicate Melt Concentrations in an Earth Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax6_1.plot(depth_earth_3, concs_mesh_earth_3, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax6_1.plot(depth_earth_4, concs_mesh_earth_4, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax6_1.fill_between(depth_earth_3, concs_mesh_earth_3, concs_mesh_earth_4, color='red', alpha=0.2)
ax6_1.grid()
ax6_1.legend(loc='upper right')


fig7 = plt.figure()
ax7_1 = fig7.add_subplot(111)
ax7_1.set_xlabel("Depth (km)")
ax7_1.set_ylabel("Concentration (moles/m$^3$)")
ax7_1.set_title("Reacted Silicate Melt Concentrations in an Earth Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax7_1.plot(depth_earth_5, concs_mesh_earth_5, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax7_1.plot(depth_earth_6, concs_mesh_earth_6, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax7_1.fill_between(depth_earth_5, concs_mesh_earth_5, concs_mesh_earth_6, color='red', alpha=0.2)
ax7_1.grid()
ax7_1.legend(loc='upper right')

fig8 = plt.figure()
ax8_1 = fig8.add_subplot(111)
ax8_1.set_xlabel("Depth (km)")
ax8_1.set_ylabel("Concentration (moles/m$^3$)")
ax8_1.set_title("Reacted Silicate Melt Concentrations in an Earth Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax8_1.plot(depth_earth_7, concs_mesh_earth_7, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax8_1.plot(depth_earth_8, concs_mesh_earth_8, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax8_1.fill_between(depth_earth_7, concs_mesh_earth_7, concs_mesh_earth_8, color='red', alpha=0.2)
ax8_1.grid()
ax8_1.legend(loc='upper right')

# partition coefficients
fig13 = plt.figure()
ax13_1 = fig13.add_subplot(111)
ax13_1.set_xlabel("Depth (km)")
ax13_1.set_ylabel("D")
ax13_1.set_title("Predicted Metal-Silicate Partition Coefficients in an Earth Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax13_1.plot(depth_earth_1, verify_D_earth_1, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax13_1.plot(depth_earth_2, verify_D_earth_2, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax13_1.fill_between(depth_earth_1, verify_D_earth_1, verify_D_earth_2, color='red', alpha=0.2)
ax13_1.grid()
ax13_1.legend(loc='upper right')

fig14 = plt.figure()
ax14_1 = fig14.add_subplot(111)
ax14_1.set_xlabel("Depth (km)")
ax14_1.set_ylabel("D")
ax14_1.set_title("Predicted Metal-Silicate Partition Coefficients in an Earth Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax14_1.plot(depth_earth_3, verify_D_earth_3, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax14_1.plot(depth_earth_4, verify_D_earth_4, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax14_1.fill_between(depth_earth_3, verify_D_earth_3, verify_D_earth_4, color='red', alpha=0.2)
ax14_1.grid()
ax14_1.legend(loc='upper right')


fig15 = plt.figure()
ax15_1 = fig15.add_subplot(111)
ax15_1.set_xlabel("Depth (km)")
ax15_1.set_ylabel("D")
ax15_1.set_title("Predicted Metal-Silicate Partition Coefficients in an Earth Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax15_1.plot(depth_earth_5, verify_D_earth_5, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax15_1.plot(depth_earth_6, verify_D_earth_6, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax15_1.fill_between(depth_earth_5, verify_D_earth_5, verify_D_earth_6, color='red', alpha=0.2)
ax15_1.grid()
ax15_1.legend(loc='upper right')

fig16 = plt.figure()
ax16_1 = fig16.add_subplot(111)
ax16_1.set_xlabel("Depth (km)")
ax16_1.set_ylabel("D")
ax16_1.set_title("Predicted Metal-Silicate Partition Coefficients in an Earth Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax16_1.plot(depth_earth_7, verify_D_earth_7, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax16_1.plot(depth_earth_8, verify_D_earth_8, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax16_1.fill_between(depth_earth_7, verify_D_earth_7, verify_D_earth_8, color='red', alpha=0.2)
ax16_1.grid()
ax16_1.legend(loc='upper right')

fig17 = plt.figure()
ax17_1 = fig17.add_subplot(111)
ax17_1.set_xlabel("Depth (km)")
ax17_1.set_ylabel("Concentration (moles/m$^3$)")
ax17_1.set_title("Reacted Silicate Melt Concentrations in a Vesta Magma Ocean (Oxidizing Model, $\eta$=10$^{-3.5}$)")
ax17_1.plot(depth_vesta_1, concs_mesh_vesta_1, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax17_1.plot(depth_vesta_2, concs_mesh_vesta_2, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax17_1.fill_between(depth_vesta_1, concs_mesh_vesta_1, concs_mesh_vesta_2, color='red', alpha=0.2)
ax17_1.grid()
ax17_1.legend(loc='upper right')

fig18 = plt.figure()
ax18_1 = fig18.add_subplot(111)
ax18_1.set_xlabel("Depth (km)")
ax18_1.set_ylabel("Concentration (moles/m$^3$)")
ax18_1.set_title("Reacted Silicate Melt Concentrations in a Vesta Magma Ocean (Reducing Model, $\eta$=10$^{-3.5}$)")
ax18_1.plot(depth_vesta_3, concs_mesh_vesta_3, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax18_1.plot(depth_vesta_4, concs_mesh_vesta_4, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax18_1.fill_between(depth_vesta_3, concs_mesh_vesta_3, concs_mesh_vesta_4, color='red', alpha=0.2)
ax18_1.grid()
ax18_1.legend(loc='upper right')


fig19 = plt.figure()
ax19_1 = fig19.add_subplot(111)
ax19_1.set_xlabel("Depth (km)")
ax19_1.set_ylabel("Concentration (moles/m$^3$)")
ax19_1.set_title("Reacted Silicate Melt Concentrations in a Vesta Magma Ocean (Oxidizing Model, $\eta$=10$^{-1.0}$)")
ax19_1.plot(depth_vesta_5, concs_mesh_vesta_5, linewidth=2.0, color='blue', label='fO$_2$=IW-0.8')
ax19_1.plot(depth_vesta_6, concs_mesh_vesta_6, linewidth=2.0, color='green', label='fO$_2$=IW-1.10')
ax19_1.fill_between(depth_vesta_5, concs_mesh_vesta_5, concs_mesh_vesta_6, color='red', alpha=0.2)
ax19_1.grid()
ax19_1.legend(loc='upper right')

fig20 = plt.figure()
ax20_1 = fig20.add_subplot(111)
ax20_1.set_xlabel("Depth (km)")
ax20_1.set_ylabel("Concentration (moles/m$^3$)")
ax20_1.set_title("Reacted Silicate Melt Concentrations in a Vesta Magma Ocean (Reducing Model, $\eta$=10$^{-1.0}$)")
ax20_1.plot(depth_vesta_7, concs_mesh_vesta_7, linewidth=2.0, color='blue', label='fO$_2$=IW-2.25')
ax20_1.plot(depth_vesta_8, concs_mesh_vesta_8, linewidth=2.0, color='green', label='fO$_2$=IW-2.45')
ax20_1.fill_between(depth_vesta_7, concs_mesh_vesta_7, concs_mesh_vesta_8, color='red', alpha=0.2)
ax20_1.grid()
ax20_1.legend(loc='upper right')


ratio_1 = [i / j for i, j in zip(list(vesta_3['D']), list(vesta_1['D']))]
ratio_2 = [i / j for i, j in zip(list(vesta_4['D']), list(vesta_2['D']))]
print(sum(ratio_1) / len(ratio_1))
print(sum(ratio_2) / len(ratio_2))

ratio_3 = [i / j for i, j in zip(concs_mesh_vesta_1, concs_mesh_vesta_5)]
ratio_4 = [i / j for i, j in zip(concs_mesh_vesta_2, concs_mesh_vesta_6)]
print(sum(ratio_3) / len(ratio_3))
print(sum(ratio_4) / len(ratio_4))


bulk_D_vesta_1 = forIterReverseD(cell_concs=concs_mesh_vesta_1, obj_concs=concs_objs_vesta_1,
                                 initial_moles_obj=0, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_1_thru_4, radius_droplet=droplet_radius)
bulk_D_vesta_2 = forIterReverseD(cell_concs=concs_mesh_vesta_2, obj_concs=concs_objs_vesta_2,
                                 initial_moles_obj=0, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_1_thru_4, radius_droplet=droplet_radius)
bulk_D_vesta_3 = forIterReverseD(cell_concs=concs_mesh_vesta_3, obj_concs=concs_objs_vesta_3,
                                 initial_moles_obj=0, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_1_thru_4, radius_droplet=droplet_radius)
bulk_D_vesta_4 = forIterReverseD(cell_concs=concs_mesh_vesta_4, obj_concs=concs_objs_vesta_4,
                                 initial_moles_obj=0, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_1_thru_4, radius_droplet=droplet_radius)
bulk_D_vesta_5 = forIterReverseD(cell_concs=concs_mesh_vesta_5, obj_concs=concs_objs_vesta_5,
                                 initial_moles_obj=0, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_5_thru_8, radius_droplet=droplet_radius)
bulk_D_vesta_6 = forIterReverseD(cell_concs=concs_mesh_vesta_6, obj_concs=concs_objs_vesta_6, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_5_thru_8, radius_droplet=droplet_radius, initial_moles_obj=0)
bulk_D_vesta_7 = forIterReverseD(cell_concs=concs_mesh_vesta_7, obj_concs=concs_objs_vesta_7, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_5_thru_8, radius_droplet=droplet_radius, initial_moles_obj=0)
bulk_D_vesta_8 = forIterReverseD(cell_concs=concs_mesh_vesta_8, obj_concs=concs_objs_vesta_8, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=vesta_vol_mesh_5_thru_8, radius_droplet=droplet_radius, initial_moles_obj=0)
bulk_D_earth_1 = forIterReverseD(cell_concs=concs_mesh_earth_1, obj_concs=concs_objs_earth_1, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_1_thru_4, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_1)
bulk_D_earth_2 = forIterReverseD(cell_concs=concs_mesh_earth_2, obj_concs=concs_objs_earth_2, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_1_thru_4, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_2)
bulk_D_earth_3 = forIterReverseD(cell_concs=concs_mesh_earth_3, obj_concs=concs_objs_earth_3, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_1_thru_4, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_3)
bulk_D_earth_4 = forIterReverseD(cell_concs=concs_mesh_earth_4, obj_concs=concs_objs_earth_4, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_5_thru_8, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_4)
bulk_D_earth_5 = forIterReverseD(cell_concs=concs_mesh_earth_5, obj_concs=concs_objs_earth_5, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_5_thru_8, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_5)
bulk_D_earth_6 = forIterReverseD(cell_concs=concs_mesh_earth_6, obj_concs=concs_objs_earth_6, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_5_thru_8, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_6)
bulk_D_earth_7 = forIterReverseD(cell_concs=concs_mesh_earth_7, obj_concs=concs_objs_earth_7, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_5_thru_8, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_7)
bulk_D_earth_8 = forIterReverseD(cell_concs=concs_mesh_earth_8, obj_concs=concs_objs_earth_8, initial_moles_mesh=0.27950089725326804,
                                 volume_mesh=earth_vol_mesh_5_thru_8, radius_droplet=earth_droplet_radius, initial_moles_obj=moles_per_droplet_earth_8)

fig21 = plt.figure()
fig22 = plt.figure()
ax21 = fig21.add_subplot(111)
ax22 = fig22.add_subplot(111)
ax21.plot(depth_vesta_1, bulk_D_vesta_1, linewidth=2.0, label='Vesta Model 1')
ax21.plot(depth_vesta_2, bulk_D_vesta_2, linewidth=2.0, label='Vesta Model 2')
ax21.plot(depth_vesta_3, bulk_D_vesta_3, linewidth=2.0, label='Vesta Model 3')
ax21.plot(depth_vesta_4, bulk_D_vesta_4, linewidth=2.0, label='Vesta Model 4')
ax21.plot(depth_vesta_5, bulk_D_vesta_5, linewidth=2.0, label='Vesta Model 5')
ax21.plot(depth_vesta_6, bulk_D_vesta_6, linewidth=2.0, label='Vesta Model 6')
ax21.plot(depth_vesta_7, bulk_D_vesta_7, linewidth=2.0, label='Vesta Model 7')
ax21.plot(depth_vesta_8, bulk_D_vesta_8, linewidth=2.0, label='Vesta Model 8')
ax22.plot(depth_earth_1, bulk_D_earth_1, linewidth=2.0, label='Earth Model 1')
ax22.plot(depth_earth_2, bulk_D_earth_2, linewidth=2.0, label='Earth Model 2')
ax22.plot(depth_earth_3, bulk_D_earth_3, linewidth=2.0, label='Earth Model 3')
ax22.plot(depth_earth_4, bulk_D_earth_4, linewidth=2.0, label='Earth Model 4')
ax22.plot(depth_earth_5, bulk_D_earth_5, linewidth=2.0, label='Earth Model 5')
ax22.plot(depth_earth_6, bulk_D_earth_6, linewidth=2.0, label='Earth Model 6')
ax22.plot(depth_earth_7, bulk_D_earth_7, linewidth=2.0, label='Earth Model 7')
ax22.plot(depth_earth_8, bulk_D_earth_8, linewidth=2.0, label='Earth Model 8')
ax22.axhline(40, linewidth=2.0, linestyle="--", color='red', label='Observed Earth Bulk D')
ax21.set_title("Bulk D For Modeled Vesta")
ax21.set_xlabel("Depth (km)")
ax21.set_ylabel("Bulk D")
ax22.set_title("Bulk D For Modeled Earth")
ax22.set_xlabel("Depth (km)")
ax22.set_ylabel("Bulk D")
ax21.grid()
ax22.grid()
ax21.legend(loc='upper right')
ax22.legend(loc='upper right')

fig23 = plt.figure()
fig24 = plt.figure()
ax23 = fig23.add_subplot(111)
ax24 = fig24.add_subplot(111)
ax23.plot(depth_vesta_1, bulk_D_vesta_1, linewidth=2.0, label='Vesta Model 1')
ax23.plot(depth_vesta_2, bulk_D_vesta_2, linewidth=2.0, label='Vesta Model 2')
ax23.plot(depth_vesta_3, bulk_D_vesta_3, linewidth=2.0, label='Vesta Model 3')
ax23.plot(depth_vesta_4, bulk_D_vesta_4, linewidth=2.0, label='Vesta Model 4')
ax23.plot(depth_vesta_5, bulk_D_vesta_5, linewidth=2.0, label='Vesta Model 5')
ax23.plot(depth_vesta_6, bulk_D_vesta_6, linewidth=2.0, label='Vesta Model 6')
ax23.plot(depth_vesta_7, bulk_D_vesta_7, linewidth=2.0, label='Vesta Model 7')
ax23.plot(depth_vesta_8, bulk_D_vesta_8, linewidth=2.0, label='Vesta Model 8')
ax24.plot(depth_earth_1, bulk_D_earth_1, linewidth=2.0, label='Earth Model 1')
ax24.plot(depth_earth_2, bulk_D_earth_2, linewidth=2.0, label='Earth Model 2')
ax24.plot(depth_earth_3, bulk_D_earth_3, linewidth=2.0, label='Earth Model 3')
ax24.plot(depth_earth_4, bulk_D_earth_4, linewidth=2.0, label='Earth Model 4')
ax24.plot(depth_earth_5, bulk_D_earth_5, linewidth=2.0, label='Earth Model 5')
ax24.plot(depth_earth_6, bulk_D_earth_6, linewidth=2.0, label='Earth Model 6')
ax24.plot(depth_earth_7, bulk_D_earth_7, linewidth=2.0, label='Earth Model 7')
ax24.plot(depth_earth_8, bulk_D_earth_8, linewidth=2.0, label='Earth Model 8')
ax24.axhline(40, linewidth=2.0, linestyle="--", color='red', label='Observed Earth Bulk log(D)')
ax23.set_title("Bulk log(D) For Modeled Vesta")
ax23.set_xlabel("Depth (km)")
ax23.set_ylabel("Bulk log(D)")
ax24.set_title("Bulk log(D) For Modeled Earth")
ax24.set_xlabel("Depth (km)")
ax24.set_ylabel("Bulk log(D)")
ax23.set_yscale('log')
ax24.set_yscale('log')
ax23.grid()
ax24.grid()
ax23.legend(loc='lower right')
ax24.legend(loc='lower right')

fig25 = plt.figure()
ax25 = fig25.add_subplot(111)
ax25.plot(depth_earth_1[2:], bulk_D_earth_1[2:], linewidth=2.0, label='Earth Model 1')
ax25.plot(depth_earth_2[2:], bulk_D_earth_2[2:], linewidth=2.0, label='Earth Model 2')
ax25.plot(depth_earth_3[2:], bulk_D_earth_3[2:], linewidth=2.0, label='Earth Model 3')
ax25.plot(depth_earth_4[2:], bulk_D_earth_4[2:], linewidth=2.0, label='Earth Model 4')
ax25.plot(depth_earth_5[2:], bulk_D_earth_5[2:], linewidth=2.0, label='Earth Model 5')
ax25.plot(depth_earth_6[2:], bulk_D_earth_6[2:], linewidth=2.0, label='Earth Model 6')
ax25.plot(depth_earth_7[2:], bulk_D_earth_7[2:], linewidth=2.0, label='Earth Model 7')
ax25.plot(depth_earth_8[2:], bulk_D_earth_8[2:], linewidth=2.0, label='Earth Model 8')
ax25.axhline(40, linewidth=2.0, linestyle="--", color='red', label='Observed Earth Bulk log D')
ax25.set_title("Bulk D For Modeled Earth w/o Initial Disequilibrium")
ax25.set_xlabel("Depth (km)")
ax25.set_ylabel("Bulk D)")
ax25.grid()
ax25.legend(loc='lower right')

fig26 = plt.figure()
ax26 = fig26.add_subplot(111)
ax26.plot(depth_earth_1[2:], bulk_D_earth_1[2:], linewidth=2.0, label='Earth Model 1')
ax26.plot(depth_earth_2[2:], bulk_D_earth_2[2:], linewidth=2.0, label='Earth Model 2')
ax26.plot(depth_earth_3[2:], bulk_D_earth_3[2:], linewidth=2.0, label='Earth Model 3')
ax26.plot(depth_earth_4[2:], bulk_D_earth_4[2:], linewidth=2.0, label='Earth Model 4')
ax26.plot(depth_earth_5[2:], bulk_D_earth_5[2:], linewidth=2.0, label='Earth Model 5')
ax26.plot(depth_earth_6[2:], bulk_D_earth_6[2:], linewidth=2.0, label='Earth Model 6')
ax26.plot(depth_earth_7[2:], bulk_D_earth_7[2:], linewidth=2.0, label='Earth Model 7')
ax26.plot(depth_earth_8[2:], bulk_D_earth_8[2:], linewidth=2.0, label='Earth Model 8')
ax26.axhline(40, linewidth=2.0, linestyle="--", color='red', label='Observed Earth Bulk log(D)')
ax26.set_title("Bulk log(D) For Modeled Earth w/o Initial Disequilibrium")
ax26.set_xlabel("Depth (km)")
ax26.set_ylabel("Bulk log(D)")
ax26.set_yscale('log')
ax26.grid()
ax26.legend(loc='lower right')


plt.show()