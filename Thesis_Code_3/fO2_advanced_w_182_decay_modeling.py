import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp, sqrt, pi, log10
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

plt.rcParams.update({'font.size': 16})

def calcRadMantle(mass_core, density_metal, radius_body):
    volume_core = mass_core / density_metal
    radius_core = ((3 * volume_core) / (4 * pi))**(1/3)
    radius_mantle = radius_body - radius_core
    return radius_mantle, radius_core

def calcHydrostaticPressure(depth, gravity, density_melt, surface_pressure):
    p = surface_pressure + ((density_melt * gravity * depth) * 10**(-9))
    return p

def calcAdiabaticTemperature(thermal_expansivity, heat_capacity, surface_temperature, gravity, depth):
    t = surface_temperature + ((thermal_expansivity * gravity * surface_temperature * depth) / heat_capacity)
    return t

def collectCoeffsSimple(pressure, temperature):

    # pressure in GPa
    # temperature in degK
    # Cottrell et al 2009

    coeffs = {
        'alpha': 0,
          'beta': 0,
          'chi': 0,
          'delta': 0,
          'epsilon': 0
    }

    if 0 <= pressure <= 2:
        coeffs['alpha'] = 1.11
        coeffs['beta'] = -1.18
        coeffs['chi'] = -0.85
        coeffs['delta'] = 1680
        coeffs['epsilon'] = 487

    elif 2 < pressure:
        coeffs['alpha'] = 1.05
        coeffs['beta'] = -1.10
        coeffs['chi'] = -0.84
        coeffs['delta'] = 3588
        coeffs['epsilon'] = -102

    return coeffs

def partition(pressure, temperature, deltaIW):
    nbo_t = 2.6
    coeffs = collectCoeffsSimple(pressure=pressure, temperature=temperature)
    alpha = coeffs['alpha']
    beta = coeffs['beta']
    chi = coeffs['chi']
    delta = coeffs['delta']
    epsilon = coeffs['epsilon']
    logD = alpha + (beta * deltaIW) + (chi * nbo_t) + (delta * (1/temperature)) + (epsilon * (pressure/temperature))
    D = 10**logD
    return D

def calcEpsilon182W(w182_at_time, w184_at_time, terretrial_standard):
    epsilon = (((w182_at_time / w184_at_time) / terretrial_standard) - 1) * (10**4)
    return epsilon


def decayW182(timestep, core_formation_max_time, inf_time, w182_at_wt, hf182_half_life, hf182_at_wt,
              initial_hf182_conc, mass_vesta, mass_vesta_core, density_metal, density_melt, fO2,
              temperature_surf, pressure_surf, radius_body, gravity, thermal_expansivity, heat_capacity):

    core_frac_per_timestep = timestep / core_formation_max_time
    decay_const = log(0.5) / hf182_half_life

    fraction_core_accumulated = [0]
    core_mass_added = [0]
    mantle_mass_depleted = [0]
    mass_core = [0]
    mantle_mass = [mass_vesta]

    initial_bulk_moles_w182 = (((initial_hf182_conc * (10**-9)) * mass_vesta) * 1000) / hf182_at_wt
    initial_radius_mantle, initial_radius_core = calcRadMantle(mass_core=mass_core[0], density_metal=density_metal,
                                                               radius_body=radius_body)
    initial_temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                           depth=radius_body, gravity=gravity,
                                           heat_capacity=heat_capacity,
                                           surface_temperature=temperature_surf)
    initial_pressure = calcHydrostaticPressure(depth=radius_body, density_melt=density_melt, gravity=gravity,
                                       surface_pressure=pressure_surf)
    initial_dist_coeff = partition(pressure=initial_pressure, temperature=initial_temperature, deltaIW=fO2)

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
    radius_mantle = [initial_radius_mantle]
    radius_core = [initial_radius_core]
    cmb_d = [initial_dist_coeff]
    temperatures = [initial_temperature]
    pressures = [initial_pressure]

    max_modeling_time_range = list(np.arange(0, inf_time + timestep, timestep))

    for time_index, time in enumerate(max_modeling_time_range):
        if not time_index == 0:
            if time <= core_formation_max_time:
                fraction_core_bulk = fraction_core_accumulated[-1] + core_frac_per_timestep
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time
                radius_mantle_at_time, radius_core_at_time = calcRadMantle(mass_core=mass_core_at_time,
                                                                           density_metal=density_metal,
                                                                           radius_body=radius_body)
                temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                                       depth=radius_mantle_at_time, gravity=gravity,
                                                       heat_capacity=heat_capacity,
                                                       surface_temperature=temperature_surf)
                pressure = calcHydrostaticPressure(depth=radius_mantle_at_time, density_melt=density_melt, gravity=gravity,
                                                   surface_pressure=pressure_surf)
                partition_coeff = partition(pressure=pressure, temperature=temperature, deltaIW=fO2)


                radius_mantle.append(radius_mantle_at_time)
                radius_core.append(radius_core_at_time)
                temperatures.append(temperature)
                pressures.append(pressure)
                cmb_d.append(partition_coeff)
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

                radius_mantle_at_time, radius_core_at_time = calcRadMantle(mass_core=mass_core_at_time,
                                                                           density_metal=density_metal,
                                                                           radius_body=radius_body)
                temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                                       depth=radius_mantle_at_time, gravity=gravity,
                                                       heat_capacity=heat_capacity,
                                                       surface_temperature=temperature_surf)
                pressure = calcHydrostaticPressure(depth=radius_mantle_at_time, density_melt=density_melt,
                                                   gravity=gravity,
                                                   surface_pressure=pressure_surf)
                partition_coeff = partition(pressure=pressure, temperature=temperature, deltaIW=fO2)

                radius_mantle.append(radius_mantle_at_time)
                radius_core.append(radius_core_at_time)
                cmb_d.append(partition_coeff)
                temperatures.append(temperature)
                pressures.append(pressure)

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

    bulk_moles_182w_added_at_time = [i / w182_at_wt for i in bulk_mass_182w_added]
    mantle_moles_182w_added_at_time = [i / w182_at_wt for i in mantle_mass_182w_added]
    core_moles_182w_added_at_time = [i / w182_at_wt for i in core_mass_182w_added]
    bulk_moles_182w = [i / w182_at_wt for i in bulk_mass_182w]
    core_bulk_moles_182w = [i / w182_at_wt for i in bulk_core_mass_182w]
    mantle_moles_182w = [i / w182_at_wt for i in bulk_mantle_mass_182w]



    return fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core, mantle_mass, moles_182hf, \
           bulk_mass_182w, bulk_conc_182w, bulk_mass_182w_added, bulk_conc_182w_added, mantle_conc_182w_added, \
           core_conc_182w_added, core_mass_182w_added, mantle_mass_182w_added, bulk_core_mass_182w, \
           bulk_mantle_mass_182w, bulk_mass_182w_check, core_bulk_conc_182w, mantle_bulk_conc_182w, radius_mantle, \
           radius_core, cmb_d, temperatures, pressures, bulk_moles_182w_added_at_time, mantle_moles_182w_added_at_time, \
           core_moles_182w_added_at_time, bulk_moles_182w, core_bulk_moles_182w, mantle_moles_182w



def decayW184(initial_conc_w184, mass_vesta, mass_vesta_core, core_formation_max_time, inf_time, timestep,
              density_metal, density_melt, fO2, temperature_surf, pressure_surf, radius_body, gravity,
              thermal_expansivity, heat_capacity, w184_at_wt):

    core_frac_per_timestep = timestep / core_formation_max_time
    bulk_mass_w184 = (initial_conc_w184 * (10**-9)) * mass_vesta
    bulk_conc_w184 = (bulk_mass_w184 / mass_vesta) * (10**9)

    fraction_core_accumulated = [0]
    core_mass_added = [0]
    mantle_mass_depleted = [0]
    mass_core = [0]
    mantle_mass = [mass_vesta]

    initial_radius_mantle, initial_radius_core = calcRadMantle(mass_core=mass_core[0], density_metal=density_metal,
                                                               radius_body=radius_body)
    initial_temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                                   depth=radius_body, gravity=gravity,
                                                   heat_capacity=heat_capacity,
                                                   surface_temperature=temperature_surf)
    initial_pressure = calcHydrostaticPressure(depth=radius_body, density_melt=density_melt, gravity=gravity,
                                               surface_pressure=pressure_surf)
    initial_dist_coeff = partition(pressure=initial_pressure, temperature=initial_temperature, deltaIW=fO2)

    core_mass_w184_added = [0]
    current_mantle_mass_w184 = [bulk_mass_w184]
    core_mass_w184_at_time = [0]
    mantle_mass_w184_at_time = [bulk_mass_w184]
    bulk_mass_w184_check = [bulk_mass_w184]
    core_bulk_conc_184w = [0]
    mantle_bulk_conc_184w = [initial_conc_w184]
    bulk_conc_184w_at_time = [bulk_conc_w184]
    bulk_mass_w184_at_time = [bulk_mass_w184]
    radius_mantle = [initial_radius_mantle]
    radius_core = [initial_radius_core]
    cmb_d = [initial_dist_coeff]
    temperatures = [initial_temperature]
    pressures = [initial_pressure]

    max_modeling_time_range = list(np.arange(0, inf_time + timestep, timestep))

    for time_index, time in enumerate(max_modeling_time_range):
        if not time_index == 0:
            if time < core_formation_max_time:
                fraction_core_bulk = fraction_core_accumulated[-1] + core_frac_per_timestep
                mass_core_at_time = fraction_core_bulk * mass_vesta_core
                mass_core_added_at_time = mass_core_at_time - mass_core[-1]
                mass_mantle_at_time = mass_vesta - mass_core_at_time
                mass_mantle_depleted_at_time = mantle_mass[-1] - mass_mantle_at_time

                radius_mantle_at_time, radius_core_at_time = calcRadMantle(mass_core=mass_core_at_time,
                                                                           density_metal=density_metal,
                                                                           radius_body=radius_body)
                temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                                       depth=radius_mantle_at_time, gravity=gravity,
                                                       heat_capacity=heat_capacity,
                                                       surface_temperature=temperature_surf)
                pressure = calcHydrostaticPressure(depth=radius_mantle_at_time, density_melt=density_melt,
                                                   gravity=gravity,
                                                   surface_pressure=pressure_surf)
                partition_coeff = partition(pressure=pressure, temperature=temperature, deltaIW=fO2)

                radius_mantle.append(radius_mantle_at_time)
                radius_core.append(radius_core_at_time)
                cmb_d.append(partition_coeff)
                temperatures.append(temperature)
                pressures.append(pressure)


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

                radius_mantle_at_time, radius_core_at_time = calcRadMantle(mass_core=mass_core_at_time,
                                                                           density_metal=density_metal,
                                                                           radius_body=radius_body)
                temperature = calcAdiabaticTemperature(thermal_expansivity=thermal_expansivity,
                                                       depth=radius_mantle_at_time, gravity=gravity,
                                                       heat_capacity=heat_capacity,
                                                       surface_temperature=temperature_surf)
                pressure = calcHydrostaticPressure(depth=radius_mantle_at_time, density_melt=density_melt,
                                                   gravity=gravity,
                                                   surface_pressure=pressure_surf)
                partition_coeff = partition(pressure=pressure, temperature=temperature, deltaIW=fO2)

                radius_mantle.append(radius_mantle_at_time)
                radius_core.append(radius_core_at_time)
                cmb_d.append(partition_coeff)
                temperatures.append(temperature)
                pressures.append(pressure)

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

    bulk_moles_182w_added_at_time = [i / w182_at_wt for i in bulk_mass_182w_added]
    mantle_moles_184w_remaining_at_time = [i / w184_at_wt for i in mantle_mass_w184_at_time]
    core_moles_184w_added_at_time = [i / w184_at_wt for i in core_mass_w184_added]
    bulk_moles_184w = [i / w184_at_wt for i in bulk_mass_w184_at_time]
    core_bulk_moles_184w = [i / w184_at_wt for i in core_mass_w184_at_time]
    mantle_moles_184w = [i / w184_at_wt for i in mantle_mass_w184_at_time]



    return fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core,  mantle_mass, \
           core_mass_w184_added, current_mantle_mass_w184, core_mass_w184_at_time, mantle_mass_w184_at_time, \
           bulk_mass_w184_check, core_bulk_conc_184w, mantle_bulk_conc_184w, bulk_conc_184w_at_time, \
           bulk_mass_w184_at_time, radius_mantle, radius_core, cmb_d, temperatures, pressures, \
           bulk_moles_182w_added_at_time, mantle_moles_184w_remaining_at_time, core_moles_184w_added_at_time, \
           bulk_moles_184w, core_bulk_moles_184w, mantle_moles_184w


density_metal = 7800
density_melt = 3750
timestep = 250000
core_formation_max_time = 5 * (10**6)
inf_time = 100 * (10**6)
w182_at_wt = 183.84
hf182_half_life = 8.9 * (10**6)
hf182_at_wt = 178.49
mass_vesta = 2.59076 * (10**20)
radius_vesta = (262.7) * 1000
volume_vesta = (4/3) * pi * (radius_vesta**3)
radius_vesta_core = 113 * 1000
volume_vesta_core = (4/3) * pi * (radius_vesta_core**3)
mass_vesta_core = density_metal * volume_vesta_core
mass_vesta_mantle = mass_vesta - mass_vesta_core
# initial_hf182_conc = [16.5605, 16.57299, 16.58399, 16.584425]
# initial_conc_w184 = [23.29356, 23.47345, 23.85845, 23.866664]
initial_hf182_conc = [16.4539, 16.5605, 16.57299, 16.58399, 16.584425, 16.584556]
initial_conc_w184 = [20.0194, 23.29356, 23.47345, 23.85845, 23.866664, 23.87737]
# terrestrial_standard = 0.864900  # Kleine et al 2017
terrestrial_standard = 0.864680  # Kleine et al 2004
time_list = list(np.arange(0, inf_time + timestep, timestep))
time_list_ma = [i / (10**6) for i in list(np.arange(0, inf_time + timestep, timestep))]
core_formation_max_time_index = time_list.index(core_formation_max_time)
gravity = 0.25
thermal_expansivity = 6 * (10**(-5))
heat_capacity = (10**3)
# fO2 = [-0.8, -1.10, -2.25, -2.45]
fO2 = [0.5, -0.8, -1.10, -2.25, -2.45, -3.5]
temperature_surf = 2000
pressure_surf = 0
radius_body = (262.7 * 1000)

avg_eucrite_w182_w184_ratio = 0.866555125
# avg_eucrite_epsilon_w182 = 19.13660539


print(
    "Mass Vesta Core: {}\n"
    "Mass Vesta Mantle: {}\n".format(mass_vesta_core, mass_vesta_mantle)
)

bulk_d_w182_list = []
bulk_d_w184_list = []
bulk_d_w_list = []
cmbd_fO2 = []
bulk_core_mass_182w_list = []
bulk_mantle_mass_182w_list = []
bulk_mass_182w_list = []
bulk_core_mass_184w_list = []
bulk_mantle_mass_184w_list = []
mantle_w182_w184_ratios = []
core_w182_w184_ratios = []
bulk_w182_w184_ratios = []
epsilon_core = []
epsilon_mantle = []
epsilon_bulk = []
bulk_mantle_conc_w182_list = []
bulk_core_conc_w182_list = []
bulk_mantle_conc_w184_list = []
bulk_core_conc_w184_list = []



for index, i in enumerate(fO2):
    hf182_conc = initial_hf182_conc[index]
    w_184_conc = initial_conc_w184[index]

    fraction_core_accumulated, core_mass_added, mantle_mass_depleted, mass_core, mantle_mass, moles_182hf, \
               bulk_mass_182w, bulk_conc_182w, bulk_mass_182w_added, bulk_conc_182w_added, mantle_conc_182w_added, \
               core_conc_182w_added, core_mass_182w_added, mantle_mass_182w_added, bulk_core_mass_182w, \
               bulk_mantle_mass_182w, bulk_mass_182w_check, core_bulk_conc_182w, mantle_bulk_conc_182w, radius_mantle, \
                radius_core, cmb_d, temperatures, pressures, bulk_moles_182w_added_at_time, \
                mantle_moles_182w_added_at_time, core_moles_182w_added_at_time, bulk_moles_182w, core_bulk_moles_182w, \
                mantle_moles_182w = \
                    decayW182(timestep=timestep, core_formation_max_time=core_formation_max_time, inf_time=inf_time,
                  w182_at_wt=w182_at_wt, hf182_half_life=hf182_half_life, hf182_at_wt=hf182_at_wt,
                  initial_hf182_conc=hf182_conc, mass_vesta=mass_vesta, mass_vesta_core=mass_vesta_core,
                  density_metal=density_metal, density_melt=density_melt, fO2=i, temperature_surf=temperature_surf,
                  pressure_surf=pressure_surf, radius_body=radius_body, gravity=gravity,
                  thermal_expansivity=thermal_expansivity, heat_capacity=heat_capacity)

    fraction_core_accumulated2, core_mass_added2, mantle_mass_depleted2, mass_core2,  mantle_mass2, \
               core_mass_w184_added, current_mantle_mass_w184, core_mass_w184_at_time, mantle_mass_w184_at_time, \
               bulk_mass_w184_check, core_bulk_conc_184w, mantle_bulk_conc_184w, bulk_conc_184w_at_time, bulk_mass_w184_at_time, radius_mantle2, \
                radius_core2, cmb_d2, temperatures2, pressures2, bulk_moles_184w_added_at_time, \
                mantle_moles_184w_remaining_at_time, core_moles_184w_added_at_time, bulk_moles_184w, core_bulk_moles_184w, \
                mantle_moles_184w = \
                decayW184(initial_conc_w184=w_184_conc, mass_vesta=mass_vesta, mass_vesta_core=mass_vesta_core,
                  core_formation_max_time=core_formation_max_time, inf_time=inf_time, timestep=timestep,
                  density_metal=density_metal, density_melt=density_melt, fO2=i, temperature_surf=temperature_surf,
                  pressure_surf=pressure_surf, radius_body=radius_body, gravity=gravity,
                  thermal_expansivity=thermal_expansivity, heat_capacity=heat_capacity, w184_at_wt=w182_at_wt)

    bulk_epsilon_w182_vesta = [calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
                               for i, j in zip(bulk_conc_182w, bulk_conc_184w_at_time)]
    mantle_epsilon_w182_vesta = [
        calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
        for i, j in zip(mantle_bulk_conc_182w, mantle_bulk_conc_184w)]
    core_epsilon_w182_vesta = [calcEpsilon182W(w182_at_time=i, w184_at_time=j, terretrial_standard=terrestrial_standard)
                               for i, j in zip(core_bulk_conc_182w[1:], core_bulk_conc_184w[1:])]

    bulk_d_182w = [i / j for i, j in zip(bulk_core_mass_182w[1:], bulk_mantle_mass_182w[1:])]
    bulk_d_184w = [i / j for i, j in zip(core_mass_w184_at_time[1:], mantle_mass_w184_at_time[1:])]
    bulk_d_w = [(i + j) / (l + k)  for i, j, l, k in zip(bulk_core_mass_182w[1:], bulk_mantle_mass_182w[1:],
                                                         core_mass_w184_at_time, mantle_mass_w184_at_time)]

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

    bulk_d_w182_list.append(bulk_d_182w)
    bulk_d_w184_list.append(bulk_d_184w)
    bulk_d_w_list.append(bulk_d_w)
    bulk_core_mass_182w_list.append(bulk_core_mass_182w)
    bulk_mantle_mass_182w_list.append(bulk_mantle_mass_182w)
    bulk_core_mass_184w_list.append(core_mass_w184_at_time)
    bulk_mantle_mass_184w_list.append(mantle_mass_w184_at_time)
    bulk_mass_182w_list.append(bulk_mass_182w)
    cmbd_fO2.append(cmb_d)
    core_w182_w184_ratios.append([i / j for i, j in zip(bulk_core_mass_182w[1:], core_mass_w184_at_time[1:])])
    mantle_w182_w184_ratios.append([i / j for i, j in zip(bulk_mantle_mass_182w, mantle_mass_w184_at_time)])
    bulk_w182_w184_ratios.append([i / j for i, j in zip(bulk_mass_182w, bulk_mass_w184_at_time)])
    epsilon_bulk.append(bulk_epsilon_w182_vesta)
    epsilon_core.append(core_epsilon_w182_vesta)
    epsilon_mantle.append(mantle_epsilon_w182_vesta)
    bulk_mantle_conc_w182_list.append(mantle_bulk_conc_182w)
    bulk_core_conc_w182_list.append(core_bulk_conc_182w)
    bulk_mantle_conc_w184_list.append(mantle_bulk_conc_184w)
    bulk_core_conc_w184_list.append(core_bulk_conc_184w)

    print("fO2: {}, 182W/184W: {}, Eucrite 182W/184W: {}\n".format(i, mantle_w182_w184_ratios[-1][-1], avg_eucrite_w182_w184_ratio))


w182_mantle_dict = {'time': time_list[0:core_formation_max_time_index + 1]}
w184_mantle_dict = {'time': time_list[0:core_formation_max_time_index + 1]}
for index, i in enumerate(fO2):
    temp_mantle_w182 = {'fO2_{}'.format(i): bulk_mantle_conc_w182_list[index][0:core_formation_max_time_index + 1]}
    temp_mantle_w184 = {'fO2_{}'.format(i): bulk_mantle_conc_w184_list[index][0:core_formation_max_time_index + 1]}
    w182_mantle_dict.update(temp_mantle_w182)
    w184_mantle_dict.update(temp_mantle_w184)
w182_df = pd.DataFrame(w182_mantle_dict)
w184_df = pd.DataFrame(w184_mantle_dict)
w182_df.to_csv("182w_mantle.csv")
w184_df.to_csv("184w_mantle.csv")



# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.plot(time_list_ma, bulk_epsilon_w182_vesta, linewidth=2.0, label='Bulk')
# ax1.plot(time_list_ma, mantle_epsilon_w182_vesta, linewidth=2.0, label='Mantle')
# ax1.plot(time_list_ma[1:], core_epsilon_w182_vesta, linewidth=2.0, label='Core')
# ax1.axhline(avg_eucrite_epsilon_w182, linestyle='--', color='red', label="Avg. Eucrite")
# ax1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax1.grid()
# ax1.set_title("$\epsilon^{182}$W on Vesta Through Time")
# ax1.set_xlabel("Time (Ma)")
# ax1.set_ylabel("$\epsilon^{182}$W")
# ax1.legend(loc='upper left')


fig2 = plt.figure()
ax2_0 = fig2.add_subplot(111)
ax2_1 = fig2.add_subplot(211)
ax2_2 = fig2.add_subplot(212)
for index, i in enumerate(fO2):
    ax2_1.plot(time_list_ma, bulk_mantle_mass_182w_list[index], linewidth=2.0, label='Mantle, fO$_2$ = IW{}'.format(i))
    ax2_2.plot(time_list_ma, bulk_core_mass_184w_list[index], linewidth=2.0, label='Core, fO$_2$ = IW{}'.format(i))
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
ax2_0.set_ylabel("Mass (kg)")
ax2_0.set_title("$^{182}$W Mass on Vesta over Time")


# fig3 = plt.figure()
# ax3_0 = fig3.add_subplot(111)
# ax3_0.plot(time_list_ma, mass_core, linewidth=2.0, label='Core')
# ax3_0.plot(time_list_ma, mantle_mass, linewidth=2.0, label='Mantle')
# ax3_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax3_0.grid()
# ax3_0.set_title("Core and Mantle Mass on Vesta Over Time")
# ax3_0.set_xlabel("Time (Ma)")
# ax3_0.set_ylabel("Mass (kg)")
# ax3_0.legend(loc='lower right')

fig4 = plt.figure()
ax4_0 = fig4.add_subplot(111)
ax4_1 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
for index, i in enumerate(fO2):
    ax4_1.plot(time_list_ma, bulk_mantle_mass_184w_list[index], linewidth=2.0, label='Mantle, fO$_2$ = IW{}'.format(i))
    ax4_2.plot(time_list_ma, bulk_core_mass_184w_list[index], linewidth=2.0, label='Core, fO$_2$ = IW{}'.format(i))
ax4_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax4_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax4_1.grid()
ax4_2.grid()
ax4_0.spines['top'].set_color('none')
ax4_0.spines['bottom'].set_color('none')
ax4_0.spines['left'].set_color('none')
ax4_0.spines['right'].set_color('none')
ax4_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax4_0.xaxis.labelpad = 20
ax4_0.yaxis.labelpad = 20
ax4_0.set_xlabel("Time (Ma)")
ax4_0.set_ylabel("Mass (kg)")
ax4_0.set_title("$^{184}$W Mass on Vesta over Time")
ax4_1.legend(loc='upper right')
ax4_2.legend(loc='lower right')

# fig5 = plt.figure()
# ax5_0 = fig5.add_subplot(111)
# ax5_0.plot(time_list_ma, bulk_ratio_conc_w182_w184, linewidth=2.0, label='Bulk')
# ax5_0.plot(time_list_ma, mantle_ratio_conc_w182_w184, linewidth=2.0, label='Mantle')
# ax5_0.plot(time_list_ma[1:], core_ratio_conc_w182_w184, linewidth=2.0, label='Core')
# ax5_0.axhline(avg_eucrite_w182_w184_ratio, linestyle='--', color='red', label="Avg. Eucrite")
# ax5_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax5_0.set_title("$^{182}$W/$^{184}$W Concentrations On Vesta Over Time")
# ax5_0.set_xlabel("Time (Ma)")
# ax5_0.set_ylabel("$^{182}$W/$^{184}$W")
# ax5_0.grid()
# ax5_0.legend(loc='lower right')
#
# fig6 = plt.figure()
# ax6_0 = fig6.add_subplot(111)
# ax6_0.plot(time_list_ma, bulk_ratio_wt_w182_w184, linewidth=2.0, label='Bulk')
# ax6_0.plot(time_list_ma, mantle_ratio_wt_w182_w184, linewidth=2.0, label='Mantle')
# ax6_0.plot(time_list_ma[1:], core_ratio_wt_w182_w184, linewidth=2.0, label='Core')
# ax6_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax6_0.set_title("$^{182}$W/$^{184}$W Mass on Vesta Over Time")
# ax6_0.set_xlabel("Time (Ma)")
# ax6_0.set_ylabel("$^{182}$W/$^{184}$W")
# ax6_0.grid()
# ax6_0.legend(loc='lower right')

# fig7 = plt.figure()
# ax7_0 = fig7.add_subplot(111)
# ax7_1 = fig7.add_subplot(211)
# ax7_2 = fig7.add_subplot(212)
# ax7_1.plot(time_list_ma, pct_mass_in_mantle_w182, linewidth=2.0, label='$^{182}$W in Mantle')
# ax7_1.plot(time_list_ma, pct_mass_in_core_w182, linewidth=2.0, label='$^{182}$W in Core')
# ax7_2.plot(time_list_ma, pct_mass_in_mantle_w184, linewidth=2.0, label='$^{184}$W in Mantle')
# ax7_2.plot(time_list_ma, pct_mass_in_core_w184, linewidth=2.0, label='$^{184}$W in Core')
# ax7_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax7_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax7_1.grid()
# ax7_2.grid()
# ax7_0.spines['top'].set_color('none')
# ax7_0.spines['bottom'].set_color('none')
# ax7_0.spines['left'].set_color('none')
# ax7_0.spines['right'].set_color('none')
# ax7_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
# ax7_0.xaxis.labelpad = 20
# ax7_0.yaxis.labelpad = 20
# ax7_0.set_ylabel("Percent (%)")
# ax7_0.set_xlabel("Time (Ma)")
# ax7_0.set_title("$^{182}$W Mass Percent in Vesta Over Time (Relative to Infinite Time $^{182}$W)")
# ax7_1.legend(loc='lower right')
# ax7_2.legend(loc='lower right')
#
# fig8 = plt.figure()
# ax8_0 = fig8.add_subplot(111)
# ax8_1 = fig8.add_subplot(211)
# ax8_2 = fig8.add_subplot(212)
# ax8_1.plot(time_list_ma, temperatures, linewidth=2.0, color='black', label='Temperature at CMB')
# ax8_2.plot(time_list_ma, pressures, linewidth=2.0, color='black', label='Pressure at CMB')
# ax8_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax8_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax8_1.set_ylabel("Temeperature (K)")
# ax8_2.set_ylabel("Pressure (GPa)")
# ax8_1.grid()
# ax8_2.grid()
# ax8_0.spines['top'].set_color('none')
# ax8_0.spines['bottom'].set_color('none')
# ax8_0.spines['left'].set_color('none')
# ax8_0.spines['right'].set_color('none')
# ax8_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
# ax8_0.xaxis.labelpad = 20
# ax8_0.yaxis.labelpad = 20
# ax8_0.set_xlabel("Time (Ma)")
# ax8_0.set_title("Adiabatic Temperature and Hydrostatic Pressure at Vesta CMB Over Time")
# ax8_1.legend(loc='upper right')
# ax8_2.legend(loc='upper right')

fig9 = plt.figure()
ax9_1 = fig9.add_subplot(111)
for index, i in enumerate(fO2):
    ax9_1.plot(time_list_ma, cmbd_fO2[index], linewidth=2.0, label='fO$_2$ = IW{}'.format(i))
ax9_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax9_1.grid()
ax9_1.set_title("Metal-Silicate Partitioning Coefficient (log(D)) at Vesta CMB Over Time")
ax9_1.set_xlabel("Time (Ma)")
ax9_1.set_ylabel("log(D)")
ax9_1.legend(loc='center right')
# ax9_1.set_yscale('log')

# fig9 = plt.figure()
# ax9_0 = fig9.add_subplot(111)
# ax9_0.plot(time_list_ma, bulk_ratio_conc_w182_w184, linewidth=2.0, label='Bulk')
# ax9_0.plot(time_list_ma, mantle_ratio_conc_w182_w184, linewidth=2.0, label='Mantle')
# ax9_0.plot(time_list_ma[1:], core_ratio_conc_w182_w184, linewidth=2.0, label='Core')
# ax9_0.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax9_0.set_title("$^{182}$W/$^{184}$W Concentrations on Vesta Over Time")
# ax9_0.set_xlabel("Time (Ma)")
# ax9_0.set_ylabel("$^{182}$W/$^{184}$W")
# ax9_0.grid()
# ax9_0.legend(loc='lower right')

fig10 = plt.figure()
ax10_0 = fig10.add_subplot(111)
ax10_1 = fig10.add_subplot(211)
ax10_2 = fig10.add_subplot(212)
for index, i in enumerate(fO2):
    ax10_1.plot(time_list_ma[1:], bulk_d_w182_list[index], linewidth=2.0, label='fO$_2$ = IW{}'.format(i))
    ax10_2.plot(time_list_ma[1:], bulk_d_w184_list[index], linewidth=2.0, label='fO$_2$ = IW{}'.format(i))
ax10_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax10_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax10_0.spines['top'].set_color('none')
ax10_0.spines['bottom'].set_color('none')
ax10_0.spines['left'].set_color('none')
ax10_0.spines['right'].set_color('none')
ax10_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax10_0.xaxis.labelpad = 20
ax10_0.yaxis.labelpad = 20
ax10_0.set_title("Bulk Partition Coefficient (D) on Vesta Over Time")
ax10_0.set_xlabel("Time (Ma)")
ax10_0.set_ylabel("D")
ax10_1.grid()
ax10_2.grid()
ax10_1.legend(loc='upper right')
ax10_2.legend(loc='lower right')
# ax10_1.set_yscale('log')
# ax10_2.set_yscale('log')


# fig11 = plt.figure()
# ax11_0 = fig11.add_subplot(111)
# ax11_1 = fig11.add_subplot(211)
# ax11_2 = fig11.add_subplot(212)
# ax11_1.plot(time_list_ma, mass_core, linewidth=2.0, color='black', label='Core Mass')
# ax11_2.plot(time_list_ma, mantle_mass, linewidth=2.0, color='black', label='Mantle Mass')
# ax11_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax11_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax11_0.spines['top'].set_color('none')
# ax11_0.spines['bottom'].set_color('none')
# ax11_0.spines['left'].set_color('none')
# ax11_0.spines['right'].set_color('none')
# ax11_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
# ax11_0.xaxis.labelpad = 20
# ax11_0.yaxis.labelpad = 20
# ax11_0.set_title("Core and Mantle Mass Evolution on Vesta Through Time")
# ax11_0.set_xlabel("Time (Ma)")
# ax11_0.set_ylabel("Mass (kg)")
# ax11_1.grid()
# ax11_2.grid()
# ax11_1.legend(loc='center right')
# ax11_2.legend(loc='center right')



fig12 = plt.figure()
ax12_0 = fig12.add_subplot(111)
ax12_1 = fig12.add_subplot(211)
ax12_2 = fig12.add_subplot(212)
ax12_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax12_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax12_1.fill_between(time_list_ma, bulk_mantle_mass_182w_list[0], bulk_mantle_mass_182w_list[1], label='Oxidizing Model')
ax12_1.fill_between(time_list_ma, bulk_mantle_mass_182w_list[2], bulk_mantle_mass_182w_list[3], label='Reducing Model')
ax12_2.fill_between(time_list_ma, bulk_mantle_mass_184w_list[0], bulk_mantle_mass_184w_list[1], label='Oxidizing Model')
ax12_2.fill_between(time_list_ma, bulk_mantle_mass_184w_list[2], bulk_mantle_mass_184w_list[3], label='Reducing Model')
# ax12_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_182w_list[0][1:]], [log10(i) for i in bulk_mantle_mass_182w_list[1][1:]], label='Oxidizing Model')
# ax12_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_182w_list[2][1:]], [log10(i) for i in bulk_mantle_mass_182w_list[3][1:]], label='Reducing Model')
# ax12_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_184w_list[0][1:]], [log10(i) for i in bulk_mantle_mass_184w_list[1][1:]], label='Oxidizing Model')
# ax12_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_184w_list[2][1:]], [log10(i) for i in bulk_mantle_mass_184w_list[3][1:]], label='Reducing Model')
ax12_1.grid()
ax12_2.grid()
ax12_1.legend(loc='lower right')
ax12_2.legend(loc='lower right')
ax12_0.spines['top'].set_color('none')
ax12_0.spines['bottom'].set_color('none')
ax12_0.spines['left'].set_color('none')
ax12_0.spines['right'].set_color('none')
ax12_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax12_0.xaxis.labelpad = 20
ax12_0.yaxis.labelpad = 20
ax12_0.set_xlabel("Time (Ma)")
ax12_0.set_ylabel("Mass (kg)")
ax12_0.set_title("$^{182}$W Mass on Vesta over Time")



fig13 = plt.figure()
ax13_0 = fig13.add_subplot(111)
ax13_1 = fig13.add_subplot(211)
ax13_2 = fig13.add_subplot(212)
ax13_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax13_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax13_1.fill_between(time_list_ma, bulk_mantle_mass_184w_list[0], bulk_mantle_mass_184w_list[1], label='Oxidizing Model')
ax13_1.fill_between(time_list_ma, bulk_mantle_mass_184w_list[2], bulk_mantle_mass_184w_list[3], label='Reducing Model')
ax13_2.fill_between(time_list_ma, bulk_core_mass_184w_list[0], bulk_core_mass_184w_list[1], label='Oxidizing Model')
ax13_2.fill_between(time_list_ma, bulk_core_mass_184w_list[2], bulk_core_mass_184w_list[3], label='Reducing Model')
# ax13_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_184w_list[0][1:]], [log10(i) for i in bulk_mantle_mass_184w_list[1][1:]], label='Oxidizing Model')
# ax13_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_mantle_mass_184w_list[2][1:]], [log10(i) for i in bulk_mantle_mass_184w_list[3][1:]], label='Reducing Model')
# ax13_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_core_mass_184w_list[0][1:]], [log10(i) for i in bulk_core_mass_184w_list[1][1:]], label='Oxidizing Model')
# ax13_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_core_mass_184w_list[2][1:]], [log10(i) for i in bulk_core_mass_184w_list[3][1:]], label='Reducing Model')
ax13_1.grid()
ax13_2.grid()
ax13_0.spines['top'].set_color('none')
ax13_0.spines['bottom'].set_color('none')
ax13_0.spines['left'].set_color('none')
ax13_0.spines['right'].set_color('none')
ax13_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax13_0.xaxis.labelpad = 20
ax13_0.yaxis.labelpad = 20
ax13_0.set_xlabel("Time (Ma)")
ax13_0.set_ylabel("Mass (kg)")
ax13_0.set_title("$^{184}$W Mass on Vesta over Time")
ax13_1.legend(loc='upper right')
ax13_2.legend(loc='lower right')

fig14 = plt.figure()
ax14_1 = fig14.add_subplot(111)
# ax14_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
# ax14_1.fill_between(time_list_ma[0:core_formation_max_time_index], [log10(i) for i in cmbd_fO2[0][0:core_formation_max_time_index]], [log10(i) for i in cmbd_fO2[1][0:core_formation_max_time_index]], label='Oxidizing Model')
# ax14_1.fill_between(time_list_ma[0:core_formation_max_time_index], [log10(i) for i in cmbd_fO2[2][0:core_formation_max_time_index]], [log10(i) for i in cmbd_fO2[3][0:core_formation_max_time_index]], label='Reducing Model')
ax14_1.fill_between(time_list_ma[0:core_formation_max_time_index], [log10(i) for i in cmbd_fO2[0][0:core_formation_max_time_index]], [log10(i) for i in cmbd_fO2[1][0:core_formation_max_time_index]], label='Oxidizing Model')
ax14_1.fill_between(time_list_ma[0:core_formation_max_time_index], [log10(i) for i in cmbd_fO2[2][0:core_formation_max_time_index]], [log10(i) for i in cmbd_fO2[3][0:core_formation_max_time_index]], label='Reducing Model')
ax14_1.grid()
ax14_1.set_title("Metal-Silicate Partitioning Coefficient (log(D))) at Vesta CMB Over Time")
ax14_1.set_xlabel("Time (Ma)")
ax14_1.set_ylabel("log(D)")
ax14_1.legend(loc='center right')
# ax14_1.set_yscale('log')

fig15 = plt.figure()
ax15_0 = fig15.add_subplot(111)
ax15_1 = fig15.add_subplot(211)
ax15_2 = fig15.add_subplot(212)
ax15_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax15_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax15_1.fill_between(time_list_ma[1:], bulk_d_w182_list[0], bulk_d_w182_list[1], label='Oxidizing Model')
ax15_1.fill_between(time_list_ma[1:], bulk_d_w182_list[2], bulk_d_w182_list[3], label='Reducing Model')
ax15_2.fill_between(time_list_ma[1:], bulk_d_w184_list[0], bulk_d_w184_list[1], label='Oxidizing Model')
ax15_2.fill_between(time_list_ma[1:], bulk_d_w184_list[2], bulk_d_w184_list[3], label='Reducing Model')
# ax15_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_d_w182_list[0]], [log10(i) for i in bulk_d_w182_list[1]], label='Oxidizing Model')
# ax15_1.fill_between(time_list_ma[1:], [log10(i) for i in bulk_d_w182_list[2]], [log10(i) for i in bulk_d_w182_list[3]], label='Reducing Model')
# ax15_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_d_w184_list[0]], [log10(i) for i in bulk_d_w184_list[1]], label='Oxidizing Model')
# ax15_2.fill_between(time_list_ma[1:], [log10(i) for i in bulk_d_w184_list[2]], [log10(i) for i in bulk_d_w184_list[3]], label='Reducing Model')
ax15_0.spines['top'].set_color('none')
ax15_0.spines['bottom'].set_color('none')
ax15_0.spines['left'].set_color('none')
ax15_0.spines['right'].set_color('none')
ax15_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax15_0.xaxis.labelpad = 20
ax15_0.yaxis.labelpad = 20
ax15_0.set_title("Bulk Partition Coefficient (D) on Vesta Over Time")
ax15_0.set_xlabel("Time (Ma)")
ax15_0.set_ylabel("D")
ax15_1.grid()
ax15_2.grid()
ax15_1.legend(loc='upper right')
ax15_2.legend(loc='lower right')
# ax15_1.set_yscale('log')
# ax15_2.set_yscale('log')

fig16 = plt.figure()
ax16_0 = fig16.add_subplot(111)
ax16_1 = fig16.add_subplot(311)
ax16_2 = fig16.add_subplot(312)
ax16_3 = fig16.add_subplot(313)
for index, i in enumerate(fO2):
    ax16_1.plot(time_list_ma[1:], core_w182_w184_ratios[index], label='fO$_2$ = IW{} (Core)'.format(i))
    ax16_2.plot(time_list_ma, mantle_w182_w184_ratios[index], label='fO$_2$ = IW{} (Mantle)'.format(i))
    ax16_3.plot(time_list_ma, bulk_w182_w184_ratios[index], label='fO$_2$ = IW{} (Bulk)'.format(i))
ax16_2.axhline(avg_eucrite_w182_w184_ratio, color='black', linestyle="--", label='Avg. Eucrite')
ax16_3.axhline(avg_eucrite_w182_w184_ratio, color='black', linestyle="--", label='Avg. Eucrite')
ax16_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax16_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax16_3.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax16_0.spines['top'].set_color('none')
ax16_0.spines['bottom'].set_color('none')
ax16_0.spines['left'].set_color('none')
ax16_0.spines['right'].set_color('none')
ax16_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax16_0.xaxis.labelpad = 20
ax16_0.yaxis.labelpad = 20
ax16_1.grid()
ax16_2.grid()
ax16_3.grid()
ax16_0.set_title("$^{182}$W/$^{184}$W On Vesta Over Time")
ax16_0.set_xlabel("Time (Ma)")
ax16_0.set_ylabel("$^{182}$W/$^{184}$W")
ax16_1.legend(loc='upper right')
ax16_2.legend(loc='lower right')
ax16_3.legend(loc='lower right')

fig17 = plt.figure()
ax17_0 = fig17.add_subplot(111)
ax17_1 = fig17.add_subplot(311)
ax17_2 = fig17.add_subplot(312)
ax17_3 = fig17.add_subplot(313)
for index, i in enumerate(fO2):
    ax17_1.plot(time_list_ma[1:], epsilon_core[index], label='fO$_2$ = IW{} (Core)'.format(i))
    ax17_2.plot(time_list_ma, epsilon_mantle[index], label='fO$_2$ = IW{} (Mantle)'.format(i))
    ax17_3.plot(time_list_ma, epsilon_bulk[index], label='fO$_2$ = IW{} (Bulk)'.format(i))
ax17_1.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax17_2.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax17_3.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax17_0.spines['top'].set_color('none')
ax17_0.spines['bottom'].set_color('none')
ax17_0.spines['left'].set_color('none')
ax17_0.spines['right'].set_color('none')
ax17_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax17_0.xaxis.labelpad = 20
ax17_0.yaxis.labelpad = 20
ax17_1.grid()
ax17_2.grid()
ax17_3.grid()
ax17_0.set_title("$\epsilon^{182}$W On Vesta Over Time")
ax17_0.set_xlabel("Time (Ma)")
ax17_0.set_ylabel("$\epsilon^{182}$W")
ax17_1.legend(loc='upper right')
ax17_2.legend(loc='lower right')
ax17_3.legend(loc='lower right')

fig18 = plt.figure()
fig19 = plt.figure()
fig20 = plt.figure()
ax18 = fig18.add_subplot(111)
ax19 = fig19.add_subplot(111)
ax20 = fig20.add_subplot(111)
for index, i in enumerate(fO2):
    ax18.plot(time_list_ma[1:], core_w182_w184_ratios[index], label='fO$_2$ = IW{} (Core)'.format(i))
    ax19.plot(time_list_ma, mantle_w182_w184_ratios[index], label='fO$_2$ = IW{} (Mantle)'.format(i))
    ax20.plot(time_list_ma, bulk_w182_w184_ratios[index], label='fO$_2$ = IW{} (Bulk)'.format(i))
ax19.axhline(avg_eucrite_w182_w184_ratio, color='black', linestyle="--", label='Avg. Eucrite')
ax20.axhline(avg_eucrite_w182_w184_ratio, color='black', linestyle="--", label='Avg. Eucrite')
ax18.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax19.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax20.axvspan(0, 5, color='red', alpha=0.2, label='Core Formation Time')
ax18.grid()
ax19.grid()
ax20.grid()
ax18.set_title("$^{182}$W/$^{184}$W In Vesta Core Over Time")
ax19.set_title("$^{182}$W/$^{184}$W In Vesta Mantle Over Time")
ax20.set_title("$^{182}$W/$^{184}$W In Bulk Vesta Over Time")
ax18.set_xlabel("Time (Ma)")
ax18.set_ylabel("$^{182}$W/$^{184}$W")
ax19.set_xlabel("Time (Ma)")
ax19.set_ylabel("$^{182}$W/$^{184}$W")
ax20.set_xlabel("Time (Ma)")
ax20.set_ylabel("$^{182}$W/$^{184}$W")
ax18.legend(loc='upper right')
ax19.legend(loc='lower right')
ax20.legend(loc='lower right')


fig21 = plt.figure()
ax21_1 = fig21.add_subplot(111)
ax21_1.plot(list(reversed(fO2)), list(reversed(initial_hf182_conc)), linewidth=2.0, color='black', linestyle="-", label="$^{182}$Hf")
ax21_1.plot(list(reversed(fO2)), list(reversed(initial_conc_w184)), linewidth=2.0, color='black', linestyle="--", label="$^{184}$W")
ax21_1.set_title("Initial Bulk Concentration on Vesta")
ax21_1.set_xlabel("$fO_2$")
ax21_1.set_ylabel("Concentration (ppb)")
ax21_1.grid()
ax21_1.legend(loc='center left')


plt.show()
