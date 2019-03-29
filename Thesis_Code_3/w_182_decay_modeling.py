import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def calcEpsilon182W(w182_at_time, w184_at_time, terretrial_standard):
    epsilon = (((w182_at_time / w184_at_time) / terretrial_standard) - 1) * (10**4)
    return epsilon


def decayW182(timestep, core_formation_max_time, inf_time, w182_at_wt, hf182_half_life, hf182_at_wt,
              initial_hf182_conc, mass_vesta, mass_vesta_core, partition_coeff):

    core_frac_per_timestep = timestep / core_formation_max_time
    decay_const = log(0.5) / hf182_half_life

    fraction_core_accumulated = [0]
    core_mass_added = [0]
    mantle_mass_depleted = [0]
    mass_core = [0]
    mantle_mass = [mass_vesta]

    initial_bulk_moles_w182 = (((initial_hf182_conc * (10**-9)) * mass_vesta) * 1000) / hf182_at_wt

    moles_182hf = [initial_bulk_moles_w182]
    bulk_mass_182w = [0]
    bulk_conc_182w = [0]
    bulk_mass_182w_added = [0]
    bulk_conc_182w_added = [0]
    mantle_conc_182w_added = [0]
    core_conc_182w_added = [0]
    core_mass_182w_added = [0]
    mantle_mass_182w_added = [0]
    bulk_core_mass_182w = [0]
    bulk_mantle_mass_182w = [0]
    bulk_mass_182w_check = [0]
    core_bulk_conc_182w = [0]
    mantle_bulk_conc_182w = [0]

    max_modeling_time_range = list(np.arange(0, inf_time + timestep, timestep))

    for time_index, time in enumerate(max_modeling_time_range):
        if not time_index == 0:
            if time <= core_formation_max_time:
                fraction_core_bulk = fraction_core_accumulated[-1] + core_frac_per_timestep
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

                fraction_core_accumulated.append(fraction_core_bulk)
                core_mass_added.append(mass_core_added_at_time)
                mantle_mass.append(mass_mantle_at_time)
                mantle_mass_depleted.append(mass_mantle_depleted_at_time)
                mass_core.append(mass_core[-1] + mass_core_added_at_time)

                moles_182hf_remaining = moles_182hf[-1] * exp(
                    (time - max_modeling_time_range[time_index - 1]) * decay_const)
                bulk_mass_w182_at_time = (moles_182hf[0] - moles_182hf_remaining) * (w182_at_wt / 1000)  # kg
                bulk_conc_w182_at_time = (bulk_mass_w182_at_time / mass_vesta) * (10 ** 9)
                bulk_mass_w182_at_time_added = bulk_mass_w182_at_time - bulk_mass_182w[-1]
                bulk_conc_w182_at_time_added = bulk_conc_w182_at_time - bulk_conc_182w[-1]

                mantle_conc_182w_at_time_added = bulk_conc_w182_at_time_added / (partition_coeff + 1)

                bulk_mass_182w_added.append(bulk_mass_w182_at_time_added)
                bulk_conc_182w_added.append(bulk_conc_w182_at_time_added)
                mantle_conc_182w_added.append(mantle_conc_182w_at_time_added)

                core_conc_w182_at_time_added = bulk_conc_w182_at_time_added - mantle_conc_182w_at_time_added
                core_mass_w182_at_time_added = (core_conc_w182_at_time_added * (10 ** -9)) * mass_core_added_at_time
                mass_mantle_w182_at_time_added = bulk_mass_w182_at_time_added - core_mass_w182_at_time_added

                core_conc_182w_added.append(core_conc_w182_at_time_added)
                core_mass_182w_added.append(core_mass_w182_at_time_added)
                mantle_mass_182w_added.append(mass_mantle_w182_at_time_added)

                mantle_bulk_mass_182w_at_time = sum(mantle_mass_182w_added)
                core_bulk_mass_182w_at_time = sum(core_mass_182w_added)
                mantle_bulk_conc_182w_at_time = (mantle_bulk_mass_182w_at_time / mass_mantle_at_time) * (10 ** 9)
                core_bulk_conc_182w_at_time = (core_bulk_mass_182w_at_time / mass_core_at_time) * (10 ** 9)
                bulk_conc_182w_at_time = ((mantle_bulk_mass_182w_at_time + core_bulk_mass_182w_at_time) / mass_vesta) * (
                        10 ** 9)

                moles_182hf.append(moles_182hf_remaining)
                bulk_core_mass_182w.append(core_bulk_mass_182w_at_time)
                bulk_mantle_mass_182w.append(mantle_bulk_mass_182w_at_time)
                bulk_mass_182w_check.append(mantle_bulk_mass_182w_at_time + core_bulk_mass_182w_at_time)
                mantle_bulk_conc_182w.append(mantle_bulk_conc_182w_at_time)
                core_bulk_conc_182w.append(core_bulk_conc_182w_at_time)
                bulk_conc_182w.append(bulk_conc_182w_at_time)
                bulk_mass_182w.append(bulk_mass_w182_at_time)



            else:
                fraction_core_bulk = 1
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

                fraction_core_accumulated.append(fraction_core_bulk)
                core_mass_added.append(mass_core_added_at_time)
                mantle_mass.append(mass_mantle_at_time)
                mantle_mass_depleted.append(mass_mantle_depleted_at_time)
                mass_core.append(mass_core[-1] + mass_core_added_at_time)

                moles_182hf_remaining = moles_182hf[-1] * exp(
                    (time - max_modeling_time_range[time_index - 1]) * decay_const)
                bulk_mass_w182_at_time = (moles_182hf[0] - moles_182hf_remaining) * (w182_at_wt / 1000)  # kg
                bulk_conc_w182_at_time = (bulk_mass_w182_at_time / mass_vesta) * (10 ** 9)
                bulk_mass_w182_at_time_added = bulk_mass_w182_at_time - bulk_mass_182w[-1]
                bulk_conc_w182_at_time_added = bulk_conc_w182_at_time - bulk_conc_182w[-1]

                mantle_conc_182w_at_time_added = bulk_conc_w182_at_time_added / (partition_coeff + 1)

                bulk_mass_182w_added.append(bulk_mass_w182_at_time_added)
                bulk_conc_182w_added.append(bulk_conc_w182_at_time_added)
                mantle_conc_182w_added.append(mantle_conc_182w_at_time_added)

                core_conc_w182_at_time_added = bulk_conc_w182_at_time_added - mantle_conc_182w_at_time_added
                core_mass_w182_at_time_added = (core_conc_w182_at_time_added * (10 ** -9)) * mass_core_added_at_time
                mass_mantle_w182_at_time_added = bulk_mass_w182_at_time_added - core_mass_w182_at_time_added

                core_conc_182w_added.append(core_conc_w182_at_time_added)
                core_mass_182w_added.append(core_mass_w182_at_time_added)
                mantle_mass_182w_added.append(mass_mantle_w182_at_time_added)

                mantle_bulk_mass_182w_at_time = sum(mantle_mass_182w_added)
                core_bulk_mass_182w_at_time = sum(core_mass_182w_added)
                mantle_bulk_conc_182w_at_time = (mantle_bulk_mass_182w_at_time / mass_mantle_at_time) * (10 ** 9)
                core_bulk_conc_182w_at_time = (core_bulk_mass_182w_at_time / mass_core_at_time) * (10 ** 9)
                bulk_conc_182w_at_time = ((
                                                      mantle_bulk_mass_182w_at_time + core_bulk_mass_182w_at_time) / mass_vesta) * (
                                                 10 ** 9)

                moles_182hf.append(moles_182hf_remaining)
                bulk_core_mass_182w.append(core_bulk_mass_182w_at_time)
                bulk_mantle_mass_182w.append(mantle_bulk_mass_182w_at_time)
                bulk_mass_182w_check.append(mantle_bulk_mass_182w_at_time + core_bulk_mass_182w_at_time)
                mantle_bulk_conc_182w.append(mantle_bulk_conc_182w_at_time)
                core_bulk_conc_182w.append(core_bulk_conc_182w_at_time)
                bulk_conc_182w.append(bulk_conc_182w_at_time)
                bulk_mass_182w.append(bulk_mass_w182_at_time)


    return fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core, mantle_mass, moles_182hf, \
           bulk_mass_182w, bulk_conc_182w, bulk_mass_182w_added, bulk_conc_182w_added, mantle_conc_182w_added, \
           core_conc_182w_added, core_mass_182w_added, mantle_mass_182w_added, bulk_core_mass_182w, \
           bulk_mantle_mass_182w, bulk_mass_182w_check, core_bulk_conc_182w, mantle_bulk_conc_182w





def decayW184(initial_conc_w184, mass_vesta, mass_vesta_core, core_formation_max_time, inf_time, timestep,
              partition_coeff):

    core_frac_per_timestep = timestep / core_formation_max_time
    bulk_mass_w184 = (initial_conc_w184 * (10**-9)) * mass_vesta
    bulk_conc_w184 = (bulk_mass_w184 / mass_vesta) * (10**9)

    fraction_core_accumulated = [0]
    core_mass_added = [0]
    mantle_mass_depleted = [0]
    mass_core = [0]
    mantle_mass = [mass_vesta]

    core_mass_w184_added = [0]
    current_mantle_mass_w184 = [bulk_mass_w184]
    core_mass_w184_at_time = [0]
    mantle_mass_w184_at_time = [bulk_mass_w184]
    bulk_mass_w184_check = [bulk_mass_w184]
    core_bulk_conc_184w = [0]
    mantle_bulk_conc_184w = [initial_conc_w184]
    bulk_conc_184w_at_time = [bulk_conc_w184]
    bulk_mass_w184_at_time = [bulk_mass_w184]

    max_modeling_time_range = list(np.arange(0, inf_time + timestep, timestep))

    for time_index, time in enumerate(max_modeling_time_range):
        if not time_index == 0:
            if time < core_formation_max_time:
                fraction_core_bulk = fraction_core_accumulated[-1] + core_frac_per_timestep
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

                fraction_core_accumulated.append(fraction_core_bulk)
                core_mass_added.append(mass_core_added_at_time)
                mantle_mass.append(mass_mantle_at_time)
                mantle_mass_depleted.append(mass_mantle_depleted_at_time)
                mass_core.append(mass_core_at_time)

                mantle_conc_w184_at_time = bulk_conc_w184 / (partition_coeff + 1)
                core_conc_w184_at_time = bulk_conc_w184 - mantle_conc_w184_at_time
                core_mass_w184_added_at_time = (core_conc_w184_at_time * (10**-9)) * mass_core_added_at_time
                mantle_mass_w184_remaining_at_time = current_mantle_mass_w184[-1] - core_mass_w184_added_at_time
                bulk_core_mass_w184_at_time = sum(core_mass_w184_added) + core_mass_w184_added_at_time
                bulk_mass_w184_check_at_time = mantle_mass_w184_remaining_at_time + bulk_core_mass_w184_at_time
                core_bulk_conc_w184_at_time = (bulk_core_mass_w184_at_time / mass_core_at_time) * (10**9)
                mantle_bulk_conc_184w_at_time = (mantle_mass_w184_remaining_at_time / mass_mantle_at_time) * (10**9)

                current_mantle_mass_w184.append(mantle_mass_w184_remaining_at_time)
                core_mass_w184_added.append(core_mass_w184_added_at_time)
                core_mass_w184_at_time.append(bulk_core_mass_w184_at_time)
                mantle_mass_w184_at_time.append(mantle_mass_w184_remaining_at_time)
                bulk_mass_w184_check.append(bulk_core_mass_w184_at_time + mantle_mass_w184_remaining_at_time)
                core_bulk_conc_184w.append(core_bulk_conc_w184_at_time)
                mantle_bulk_conc_184w.append(mantle_bulk_conc_184w_at_time)
                bulk_conc_184w_at_time.append(bulk_conc_w184)
                bulk_mass_w184_at_time.append(bulk_mass_w184)

            else:
                fraction_core_bulk = 1
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

                fraction_core_accumulated.append(fraction_core_bulk)
                core_mass_added.append(mass_core_added_at_time)
                mantle_mass.append(mass_mantle_at_time)
                mantle_mass_depleted.append(mass_mantle_depleted_at_time)
                mass_core.append(mass_core_at_time)

                mantle_conc_w184_at_time = bulk_conc_w184 / (partition_coeff + 1)
                core_conc_w184_at_time = bulk_conc_w184 - mantle_conc_w184_at_time
                core_mass_w184_added_at_time = (core_conc_w184_at_time * (10 ** -9)) * mass_core_added_at_time
                mantle_mass_w184_remaining_at_time = current_mantle_mass_w184[-1] - core_mass_w184_added_at_time
                bulk_core_mass_w184_at_time = sum(core_mass_w184_added) + core_mass_w184_added_at_time
                bulk_mass_w184_check_at_time = mantle_mass_w184_remaining_at_time + bulk_core_mass_w184_at_time
                core_bulk_conc_w184_at_time = (bulk_core_mass_w184_at_time / mass_core_at_time) * (10 ** 9)
                mantle_bulk_conc_184w_at_time = (mantle_mass_w184_remaining_at_time / mass_mantle_at_time) * (10 ** 9)

                current_mantle_mass_w184.append(mantle_mass_w184_remaining_at_time)
                core_mass_w184_added.append(core_mass_w184_added_at_time)
                core_mass_w184_at_time.append(bulk_core_mass_w184_at_time)
                mantle_mass_w184_at_time.append(mantle_mass_w184_remaining_at_time)
                bulk_mass_w184_check.append(bulk_core_mass_w184_at_time + mantle_mass_w184_remaining_at_time)
                core_bulk_conc_184w.append(core_bulk_conc_w184_at_time)
                mantle_bulk_conc_184w.append(mantle_bulk_conc_184w_at_time)
                bulk_conc_184w_at_time.append(bulk_conc_w184)
                bulk_mass_w184_at_time.append(bulk_mass_w184)




    return fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core,  mantle_mass, \
           core_mass_w184_added, current_mantle_mass_w184, core_mass_w184_at_time, mantle_mass_w184_at_time, \
           bulk_mass_w184_check, core_bulk_conc_184w, mantle_bulk_conc_184w, bulk_conc_184w_at_time, \
           bulk_mass_w184_at_time



timestep = 250000
core_formation_max_time = 5 * (10**6)
inf_time = 100 * (10**6)
w182_at_wt = 183.84
hf182_half_life = 8.9 * (10**6)
hf182_at_wt = 178.49
initial_hf182_conc = 20.26
mass_vesta = 2.59 * (10**20)
mass_vesta_core = 4.662 * (10**19)
partition_coeff = 100
initial_conc_w184 = 24.068
time_list = list(np.arange(0, inf_time + timestep, timestep))
time_list_ma = [i / (10**6) for i in list(np.arange(0, inf_time + timestep, timestep))]
core_formation_max_time_index = time_list.index(core_formation_max_time)


fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core, mantle_mass, moles_182hf, \
           bulk_mass_182w, bulk_conc_182w, bulk_mass_182w_added, bulk_conc_182w_added, mantle_conc_182w_added, \
           core_conc_182w_added, core_mass_182w_added, mantle_mass_182w_added, bulk_core_mass_182w, \
           bulk_mantle_mass_182w, bulk_mass_182w_check, core_bulk_conc_182w, mantle_bulk_conc_182w = \
    decayW182(timestep=timestep, core_formation_max_time=core_formation_max_time, inf_time=inf_time,
              w182_at_wt=w182_at_wt, hf182_half_life=hf182_half_life, hf182_at_wt=hf182_at_wt,
              initial_hf182_conc=initial_hf182_conc, mass_vesta=mass_vesta, mass_vesta_core=mass_vesta_core,
              partition_coeff=partition_coeff)

# print(
#     "Fraction Core Accumulated: {}\n"
#     "Mass Core Added: {}\n"
#     "Mass Core: {}\n"
#     "Mantle Mass: {}\n"
#     "Moles 182Hf: {}\n"
#     "Bulk Mass 182W: {}\n"
#     "Bulk Conc 182W: {}\n"
#     "Bulk Mass 182W added: {}\n"
#     "Bulk Conc 182W added: {}\n"
#     "Mantle Conc 182W added: {}\n"
#     "Core Conc 182W added: {}\n"
#     "Core Mass 182W added: {}\n"
#     "Mantle mass 182W added: {}\n"
#     "Bulk Core Mass 182W: {}\n"
#     "Bulk Mantle Mass 182W: {}\n"
#     "Bulk Mass 182W check: {}\n"
#     "Core Bulk Conc 182W: {}\n"
#     "Mantle Bulk Conc 182W: {}\n".format(
#     fraction_core_accumulated[core_formation_max_time_index],
#     core_mass_added[core_formation_max_time_index],
#     mass_core[core_formation_max_time_index],
#     mantle_mass[core_formation_max_time_index],
#     moles_182hf[core_formation_max_time_index],
#     bulk_mass_182w[core_formation_max_time_index],
#     bulk_conc_182w[core_formation_max_time_index],
#     bulk_mass_182w_added[core_formation_max_time_index],
#     bulk_conc_182w_added[core_formation_max_time_index],
#     mantle_conc_182w_added[core_formation_max_time_index],
#     core_conc_182w_added[core_formation_max_time_index],
#     core_mass_182w_added[core_formation_max_time_index],
#     mantle_mass_182w_added[core_formation_max_time_index],
#     bulk_core_mass_182w[core_formation_max_time_index],
#     bulk_mantle_mass_182w[core_formation_max_time_index],
#     bulk_mass_182w_check[core_formation_max_time_index],
#     core_bulk_conc_182w[core_formation_max_time_index],
#     mantle_bulk_conc_182w[core_formation_max_time_index]
#     )
# )

fraction_core_accumulated2, core_mass_added2, mantle_mass_depleted2, mass_core2,  mantle_mass2, \
           core_mass_w184_added, current_mantle_mass_w184, core_mass_w184_at_time, mantle_mass_w184_at_time, \
           bulk_mass_w184_check, core_bulk_conc_184w, mantle_bulk_conc_184w, bulk_conc_184w_at_time, bulk_mass_w184_at_time = \
    decayW184(initial_conc_w184=initial_conc_w184, mass_vesta=mass_vesta, mass_vesta_core=mass_vesta_core,
              core_formation_max_time=core_formation_max_time, inf_time=inf_time, timestep=timestep,
              partition_coeff=partition_coeff)

# print(
#     "Fraction Core Accumulated: {}\n"
#     "Core Mass Added: {}\n"
#     "Mantle Mass Depleted: {}\n"
#     "Mass Core: {}\n"
#     "Mass Mantle: {}\n"
#     "Core Mass W184 added: {}\n"
#     "Current Mantle Mass W184: {}\n"
#     "Core Mass 184W: {}\n"
#     "Mantle Mass 184W: {}\n"
#     "Bulk Mass W184 Check: {}\n"
#     "Core Bulk Conc 184W: {}\n"
#     "Mantle Bulk Conc: {}\n".format(
#         fraction_core_accumulated2[core_formation_max_time_index],
#         core_mass_added2[core_formation_max_time_index],
#         mantle_mass_depleted2[core_formation_max_time_index],
#         mass_core2[core_formation_max_time_index],
#         mantle_mass2[core_formation_max_time_index],
#         core_mass_w184_added[core_formation_max_time_index],
#         current_mantle_mass_w184[core_formation_max_time_index],
#         core_mass_w184_at_time[core_formation_max_time_index],
#         mantle_mass_w184_at_time[core_formation_max_time_index],
#         bulk_mass_w184_check[core_formation_max_time_index],
#         core_bulk_conc_184w[core_formation_max_time_index],
#         mantle_bulk_conc_184w[core_formation_max_time_index]
#     )
# )


terrestrial_standard = 0.864900
bulk_epsilon_w182_vesta = [calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
                          for i, j in zip(bulk_conc_182w, bulk_conc_184w_at_time)]
mantle_epsilon_w182_vesta = [calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
                          for i, j in zip(mantle_bulk_conc_182w, mantle_bulk_conc_184w)]
core_epsilon_w182_vesta = [calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
                          for i, j in zip(core_bulk_conc_182w[1:], core_bulk_conc_184w[1:])]

bulk_ratio_wt_w182_w184 = [i / j for i, j in zip(bulk_mass_182w, bulk_mass_w184_at_time)]
mantle_ratio_wt_w182_w184 = [i / j for i, j in zip(bulk_mantle_mass_182w, mantle_mass_w184_at_time)]
core_ratio_wt_w182_w184 = [i / j for i, j in zip(bulk_core_mass_182w[1:], core_mass_w184_at_time[1:])]
bulk_ratio_conc_w182_w184 = [i / j for i, j in zip(bulk_conc_182w, bulk_conc_184w_at_time)]
mantle_ratio_conc_w182_w184 = [i / j for i, j in zip(mantle_bulk_conc_182w, mantle_bulk_conc_184w)]
core_ratio_conc_w182_w184 = [i / j for i, j in zip(core_bulk_conc_182w[1:], core_bulk_conc_184w[1:])]
pct_mass_in_core_w182 = [i / bulk_mass_182w[-1] for i in bulk_core_mass_182w]
pct_mass_in_mantle_w182 = [i / bulk_mass_182w[-1] for i in bulk_mantle_mass_182w]
pct_mass_in_core_w184 = [i / bulk_mass_w184_at_time[-1] for i in core_mass_w184_at_time]
pct_mass_in_mantle_w184 = [i / bulk_mass_w184_at_time[-1] for i in mantle_mass_w184_at_time]



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(time_list_ma, bulk_epsilon_w182_vesta, linewidth=2.0, label='Bulk')
ax1.plot(time_list_ma, mantle_epsilon_w182_vesta, linewidth=2.0, label='Mantle')
ax1.plot(time_list_ma[1:], core_epsilon_w182_vesta, linewidth=2.0, label='Core')
ax1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax1.grid()
ax1.set_title("$\epsilon^{182}$W on Vesta Through Time")
ax1.set_xlabel("Time (Ma)")
ax1.set_ylabel("$\epsilon^{182}$W")
ax1.legend(loc='upper left')


fig2 = plt.figure()
ax2_0 = fig2.add_subplot(111)
ax2_1 = fig2.add_subplot(211)
ax2_2 = fig2.add_subplot(212)
ax2_1.plot(time_list_ma, mantle_bulk_conc_182w, linewidth=2.0, color='black', label='Mantle')
ax2_2.plot(time_list_ma, core_bulk_conc_182w, linewidth=2.0, color='black', label='Core')
ax2_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax2_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax2_1.grid()
ax2_2.grid()
ax2_1.legend(loc='lower right')
ax2_2.legend(loc='lower right')
ax2_0.spines['top'].set_color('none')
ax2_0.spines['bottom'].set_color('none')
ax2_0.spines['left'].set_color('none')
ax2_0.spines['right'].set_color('none')
ax2_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax2_0.xaxis.labelpad = 20
ax2_0.yaxis.labelpad = 20
ax2_0.set_xlabel("Time (Ma)")
ax2_0.set_ylabel("Concentration (ppb)")
ax2_0.set_title("$^{182}$W Concentration on Vesta over Time")


fig3 = plt.figure()
ax3_0 = fig3.add_subplot(111)
ax3_0.plot(time_list_ma, mass_core, linewidth=2.0, label='Core')
ax3_0.plot(time_list_ma, mantle_mass, linewidth=2.0, label='Mantle')
ax3_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax3_0.grid()
ax3_0.set_title("Core and Mantle Mass on Vesta Over Time")
ax3_0.set_xlabel("Time (Ma)")
ax3_0.set_ylabel("Mass (kg)")
ax3_0.legend(loc='lower right')

fig4 = plt.figure()
ax4_0 = fig4.add_subplot(111)
ax4_1 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
ax4_1.plot(time_list_ma, mantle_bulk_conc_184w, linewidth=2.0, color='black', label='Mantle')
ax4_2.plot(time_list_ma, core_bulk_conc_184w, color='black', linewidth=2.0, label='Core')
ax4_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax4_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax4_1.grid()
ax4_2.grid()
ax4_1.legend(loc='lower right')
ax4_2.legend(loc='lower right')
ax4_0.spines['top'].set_color('none')
ax4_0.spines['bottom'].set_color('none')
ax4_0.spines['left'].set_color('none')
ax4_0.spines['right'].set_color('none')
ax4_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax4_0.xaxis.labelpad = 20
ax4_0.yaxis.labelpad = 20
ax4_0.set_xlabel("Time (Ma)")
ax4_0.set_ylabel("Concentration (ppb)")
ax4_0.set_title("$^{184}$W Concentration on Vesta over Time")

fig5 = plt.figure()
ax5_0 = fig5.add_subplot(111)
ax5_0.plot(time_list_ma, bulk_ratio_conc_w182_w184, linewidth=2.0, label='Bulk')
ax5_0.plot(time_list_ma, mantle_ratio_conc_w182_w184, linewidth=2.0, label='Mantle')
ax5_0.plot(time_list_ma[1:], core_ratio_conc_w182_w184, linewidth=2.0, label='Core')
ax5_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax5_0.set_title("$^{182}$W/$^{184}$W Concentrations Over Time")
ax5_0.set_xlabel("Time (Ma)")
ax5_0.set_ylabel("$^{182}$W/$^{184}$W")
ax5_0.grid()
ax5_0.legend(loc='lower right')

fig6 = plt.figure()
ax6_0 = fig6.add_subplot(111)
ax6_0.plot(time_list_ma, bulk_ratio_wt_w182_w184, linewidth=2.0, label='Bulk')
ax6_0.plot(time_list_ma, mantle_ratio_wt_w182_w184, linewidth=2.0, label='Mantle')
ax6_0.plot(time_list_ma[1:], core_ratio_wt_w182_w184, linewidth=2.0, label='Core')
ax6_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax6_0.set_title("$^{182}$W/$^{184}$W Mass Over Time")
ax6_0.set_xlabel("Time (Ma)")
ax6_0.set_ylabel("$^{182}$W/$^{184}$W")
ax6_0.grid()
ax6_0.legend(loc='lower right')

fig7 = plt.figure()
ax7_0 = fig7.add_subplot(111)
ax7_1 = fig7.add_subplot(211)
ax7_2 = fig7.add_subplot(212)
ax7_1.plot(time_list_ma, pct_mass_in_mantle_w182, linewidth=2.0, label='$^{182}$ in Mantle')
ax7_1.plot(time_list_ma, pct_mass_in_core_w182, linewidth=2.0, label='$^{182}$ in Core')
ax7_2.plot(time_list_ma, pct_mass_in_mantle_w182, linewidth=2.0, label='$^{182}$ in Mantle')
ax7_2.plot(time_list_ma, pct_mass_in_core_w182, linewidth=2.0, label='$^{182}$ in Core')



plt.show()
