import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




def decayW182(timestep, core_formation_max_time, inf_time, w182_at_wt, hf182_half_life, initial_hf182_conc, mass_vesta,
              core_frac_per_timestep, mass_vesta_core, partition_coeff):

    decay_const = log(0.5) / hf182_half_life

    fraction_core_accumulated = [0]
    core_mass_added = [0]
    mantle_mass_depleted = [0]
    mass_core = [0]
    mantle_mass = [mass_vesta]


    moles_182hf = [0]
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

    for time_index, time in max_modeling_time_range:
        if not time_index == 0:

            fraction_core_bulk = fraction_core_accumulated[-1] + core_frac_per_timestep
            mass_core_at_time = fraction_core_bulk * mass_vesta_core
            mass_core_added_at_time = mass_core_at_time - mass_core[-1]
            mass_mantle_at_time = mass_vesta - mass_core_at_time
            mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

            fraction_core_accumulated.append(fraction_core_bulk)
            bulk_mass_182w.append(mass_core_at_time)
            core_mass_added.append(mass_core_added_at_time)
            mantle_mass.append(mass_mantle_at_time)
            mantle_mass_depleted.append(mass_mantle_depleted_at_time)



            moles_182hf_remaining = moles_182hf[-1] * exp((time -  max_modeling_time_range[time_index - 1]) * decay_const)
            bulk_mass_w182_at_time = (moles_182hf[0] - moles_182hf_remaining) * (w182_at_wt / 1000)  # kg
            bulk_conc_w182_at_time = (bulk_mass_w182_at_time / mass_vesta) * (10**9)
            bulk_mass_w182_at_time_added = bulk_mass_w182_at_time - bulk_mass_182w[-1]
            bulk_conc_w182_at_time_added = bulk_conc_w182_at_time - bulk_conc_182w[-1]

            mantle_conc_182w_at_time_added = bulk_conc_w182_at_time_added / (partition_coeff + 1)

            bulk_mass_182w_added.append(bulk_mass_w182_at_time_added)
            bulk_conc_182w_added.append(bulk_conc_w182_at_time_added)
            mantle_conc_182w_added.append(mantle_conc_182w_at_time_added)

            if time <= core_formation_max_time:
                core_conc_w182_at_time_added = bulk_conc_w182_at_time_added - mantle_conc_182w_at_time_added
                core_mass_w182_at_time_added = (core_conc_w182_at_time_added * (10**-9)) * mass_core_added_at_time
                mass_mantle_w182_added = bulk_mass_w182_at_time_added - core_mass_w182_at_time_added



            else:
                core_conc_w182_at_time_added = 0
                core_mass_w182_at_time_added = 0
                mass_mantle_w182_added = bulk_mass_w182_at_time_added - core_mass_w182_at_time_added






