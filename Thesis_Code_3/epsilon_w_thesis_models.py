import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 14})


def decay(half_life, curr_nuclii, max_time, timestep, current_time, original_nuclii, rad_list):
    decay_const = log(0.5) / half_life
    if current_time < max_time:
        remaining_nuclii = curr_nuclii * exp(decay_const * timestep)
        rad_list.append(remaining_nuclii)
        return decay(half_life=half_life, curr_nuclii=remaining_nuclii, max_time=max_time, timestep=timestep,
                     current_time=current_time + timestep, original_nuclii=original_nuclii, rad_list=rad_list)
    else:
        return rad_list

def avg_vals_time(list_of_lists):
    avgs = []
    z_list = list(zip(*list_of_lists))
    for i in z_list:
        avg = sum(i) / len(i)
        avgs.append(avg)
    return avgs

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

hf_half_life = 8.9 * 10**6
al_half_life = 7.17 * 10**5
fe_half_life = 3.0 * 10**5
max_time = 100 * (10**6)
timestep = 1 * (10**6)
time_list = [i / (1 * 10**6) for i in np.arange(0, max_time + timestep, timestep)]
my_5_index = time_list.index(5)

eucrite_df = pd.read_excel("eucrites_kleine_2002.xlsx")

samples = []
decays = []
w_184_sample = []
epsilon_w = []
terrestrial_standards = []

for row in eucrite_df.index:
    sample = eucrite_df['Sample'][row]
    w_ppb = eucrite_df['W (ppb)'][row]
    hf_ppb = eucrite_df['Hf (ppb)'][row]
    hf_180_w_184 = eucrite_df['Hf (ppb)'][row]
    w_182_w_184 = eucrite_df['182W/184W'][row]
    w_183_w_184 = eucrite_df['183W/184W'][row]
    epsilon_w_final = eucrite_df['epsilon_W'][row]
    w_184 = w_ppb / (w_182_w_184 + w_183_w_184 + 1)
    w_182 = w_182_w_184 * w_184
    w_183 = w_183_w_184 * w_184
    terrestrial_standard = w_182_w_184 / ((epsilon_w_final / (10**4)) + 1)

    # assume that all 182W was initially 182Hf
    decay_hf_182 = decay(half_life=hf_half_life, curr_nuclii=w_182, timestep=timestep, max_time=max_time,
                        original_nuclii=w_182, current_time=0, rad_list=[w_182])

    samples.append(sample)
    decays.append(decay_hf_182)
    terrestrial_standards.append(terrestrial_standard)
    w_184_sample.append(w_184)
    epsilon_w.append(epsilon_w_final)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for index, i in enumerate(samples):
    if i == "Average":
        ax1.plot(time_list, decays[index], linewidth=3.0, color="black", linestyle="--", label=i)
    else:
        ax1.plot(time_list, decays[index], linewidth=2.0, label=i)
ax1.axvspan(0, 5, alpha=0.2, color='red', label='Vesta Core Formation Period')
ax1.set_xlabel("Time (Ma)")
ax1.set_ylabel("$^{182}$Hf (ppb)")
ax1.set_title("$^{182}$Hf Over Time in Eucrites")
ax1.grid()
ax1.legend(loc='upper right')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for index, i in enumerate(samples):
    if i == "Average":
        ax2.plot(time_list, [decays[index][0] - j for j in decays[index]], linewidth=3.0, color="black", linestyle="--", label=i)
    else:
        ax2.plot(time_list, [decays[index][0] - j for j in decays[index]], linewidth=2.0, label=i)
ax2.axvspan(0, 5, alpha=0.2, color='red', label='Vesta Core Formation Period')
ax2.set_xlabel("Time (Ma)")
ax2.set_ylabel("$^{182}$W (ppb)")
ax2.set_title("$^{182}$W Over Time in Eucrites")
ax2.grid()
ax2.legend(loc='upper left')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
for index, i in enumerate(samples):
    w_184 = w_184_sample[index]
    if i == "Average":
        ax3.plot(time_list, [(decays[index][0] - j) / w_184 for j in decays[index]], linewidth=3.0, color="black", linestyle="--", label=i)
    else:
        ax3.plot(time_list, [(decays[index][0] - j) / w_184 for j in decays[index]], linewidth=2.0, label=i)
ax3.axvspan(0, 5, alpha=0.2, color='red', label='Vesta Core Formation Period')
ax3.set_xlabel("Time (Ma)")
ax3.set_ylabel("$^{182}$W/$^{184}$W")
ax3.set_title("$$^{182}$W/$^{184}$W Over Time in Eucrites")
ax3.grid()
ax3.legend(loc='upper left')



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


fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for index, i in enumerate(samples):
    w_184 = w_184_sample[index]
    if i == "Average":
        ax4.plot(time_list, [k * bulk_D_vesta_1[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 1")
        ax4.plot(time_list, [k * bulk_D_vesta_2[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 2")
        ax4.plot(time_list, [k * bulk_D_vesta_3[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 3")
        ax4.plot(time_list, [k * bulk_D_vesta_4[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 4")
        ax4.plot(time_list, [k * bulk_D_vesta_5[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 5")
        ax4.plot(time_list, [k * bulk_D_vesta_6[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 6")
        ax4.plot(time_list, [k * bulk_D_vesta_7[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 7")
        ax4.plot(time_list, [k * bulk_D_vesta_8[-1] / terrestrial_standards[-1] for k in [(decays[index][0] - j) / w_184 for j in decays[index]]],
                 linewidth=2.0, label="Vesta Model 8")
ax4.axhline(epsilon_w[-1], 0, 1, linewidth=2.0, color='red', linestyle="--", label="Average $\epsilon$$^{182}$W")
ax4.axvspan(0, 5, alpha=0.2, color='red', label='Vesta Core Formation Period')
ax4.set_xlabel("Time (My)")
ax4.set_ylabel("$\epsilon$$^{182}$W")
ax4.set_title("$\epsilon$$^{182}$W in Vesta's Core")
ax4.grid()
ax4.legend(loc='upper right')


plt.show()


