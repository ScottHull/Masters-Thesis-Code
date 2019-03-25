import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 14})


def calcBulkIsotope(half_life, D, initial_time, timestep, hf_182_mantle_at_time, w_182_mantle_at_time,
                    w_182_mantle_present):
    decay_const = log(0.5) / half_life
    w_182_add_post_core_formation = w_182_mantle_present - w_182_mantle_at_time
    current_time = initial_time
    hf_182_mantle = [hf_182_mantle_at_time]
    w_182_mantle = [w_182_mantle_at_time]
    w_182_core = [D * w_182_mantle_at_time]
    while 0 < current_time <= initial_time:
        hf_182_mantle_now = hf_182_mantle[-1] / exp(timestep * decay_const)
        w_182_mantle_now = w_182_mantle[-1] * exp(timestep * decay_const)
        w_182_core_now = D * w_182_mantle_now
        hf_182_mantle.append(hf_182_mantle_now)
        w_182_mantle.append(w_182_mantle_now)
        w_182_core.append(w_182_core_now)
        current_time -= timestep
    bulk_hf_182 = hf_182_mantle_at_time + w_182_mantle_at_time + sum(w_182_core)
    return hf_182_mantle, w_182_mantle, w_182_core, bulk_hf_182


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
w_atomic_mass = 183.84

eucrite_df = pd.read_excel("eucrites_kleine_2002.xlsx")

samples = []
hf_182_decays = []
w_184_sample = []
w_183_sample = []
w_182_sample = []
w_initial_ppb_sample = []
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
    hf_182_decays.append(decay_hf_182)
    terrestrial_standards.append(terrestrial_standard)
    w_184_sample.append(w_184)
    w_183_sample.append(w_183)
    w_182_sample.append(w_182)
    w_initial_ppb_sample.append(w_ppb)
    epsilon_w.append(epsilon_w_final)

w_182_decays = []
for index, i in enumerate(samples):
    w_182_decays_temp = [(hf_182_decays[index][0] - j) for j in hf_182_decays[index]]
    w_182_decays.append(w_182_decays_temp)


hf_182_mantles = []
w_182_mantles = []
w_182_cores = []
bulk_hf_182s = []

for index, i in enumerate(samples):
    w_182_at_5my = w_182_decays[index][my_5_index]
    hf_182_at_5my = hf_182_decays[index][my_5_index]
    print(i, w_182_at_5my, hf_182_at_5my, w_182_decays[index][-1] - w_182_at_5my)

    hf_182_mantle, w_182_mantle, w_182_core, bulk_hf_182 = calcBulkIsotope(half_life=hf_half_life, D=100,
                   initial_time=(5 * 10**6), timestep=timestep, hf_182_mantle_at_time=hf_182_decays[index][my_5_index],
                   w_182_mantle_at_time=w_182_decays[index][my_5_index],
                   w_182_mantle_present=w_initial_ppb_sample[index])
    hf_182_mantles.append(hf_182_mantle)
    w_182_mantles.append(w_182_mantle)
    w_182_cores.append(w_182_core)
    bulk_hf_182s.append(bulk_hf_182)

epsilon_182W_mantles = []
epsilon_182W_cores = []

for index, i in enumerate(samples):
    terrestrial_standard = terrestrial_standards[index]
    w_184 = w_184_sample[index]
    hf_182_mantle = hf_182_mantles[index]
    w_182_mantle = w_182_mantles[index]
    w_182_core = w_182_cores[index]
    bulk_hf_182 = bulk_hf_182s[index]
    epsilon_182W_mantle_temp = [(((j / w_184) / terrestrial_standard) - 1) * 10**4 for j in w_182_mantle]
    epsilon_182W_core_temp = [(((j / w_184) / terrestrial_standard) - 1) * 10**4 for j in w_182_core]

fig1 = plt.figure()
fig2 = plt.figure()
# fig3 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
# ax3 = fig3.add_subplot(111)
for index, i in enumerate(samples):
    terrestrial_standard = terrestrial_standards[index]
    w_184 = w_184_sample[index]
    hf_182_mantle = hf_182_mantles[index]
    w_182_mantle = w_182_mantles[index]
    w_182_core = w_182_cores[index]
    bulk_hf_182 = bulk_hf_182s[index]
    epsilon_182W_mantle_temp = [(((j / w_184) / terrestrial_standard) - 1) * 10**4 for j in w_182_mantle]
    epsilon_182W_core_temp = [(((j / w_184) / terrestrial_standard) - 1) * 10**4 for j in w_182_core]
    ax1.plot(time_list[0:my_5_index + 1], epsilon_182W_mantle_temp, linewidth=2.0, label=i)
    ax2.plot(time_list[0:my_5_index + 1], epsilon_182W_core_temp, linewidth=2.0, label=i)
    # ax3.plot(time_list[0: my_5_index + 1], bulk_hf_182, linewidth=2.0, label=i)
ax1.set_xlabel("Time (Ma)")
ax1.set_ylabel("$\epsilon^{182}$W")
ax1.set_title("$\epsilon^{182}$W in Vestian Mantle with Time")
ax1.grid()
ax1.legend(loc='upper right')
ax2.set_xlabel("Time (Ma)")
ax2.set_ylabel("$\epsilon^{182}$W")
ax2.set_title("$\epsilon^{182}$W in Vestian Core with Time")
ax2.grid()
ax2.legend(loc='upper right')
# ax3.set_xlabel("Time (Ma)")
# ax3.set_ylabel("$\^{182}$Hf (ppb)")
# ax3.set_title("$^{182}$Hf$_{bulk}$ in Vesta with Time")
# ax3.grid()
# ax3.legend(loc='upper right')


fig4 = plt.figure()
ax4_0 = fig4.add_subplot(111)
ax4_1 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
for index, i in enumerate(samples):
    w_182_mantle = w_182_mantles[index]
    w_182_core = w_182_cores[index]
    ax4_1.plot(time_list[0:my_5_index + 1], [j * w_atomic_mass for j in w_182_mantle], linewidth=2.0, label=i)
    ax4_2.plot(time_list[0:my_5_index + 1], [j * w_atomic_mass for j in w_182_core], linewidth=2.0, label=i)
ax4_0.spines['top'].set_color('none')
ax4_0.spines['bottom'].set_color('none')
ax4_0.spines['left'].set_color('none')
ax4_0.spines['right'].set_color('none')
ax4_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax4_0.xaxis.labelpad = 20
ax4_0.yaxis.labelpad = 20
ax4_0.set_xlabel("Time (Ma)")
ax4_0.set_ylabel("Mass W (g)")
ax4_0.set_title("182W Mass for Vestian Core & Mantle")
ax4_1.grid()
ax4_2.grid()

plt.show()


