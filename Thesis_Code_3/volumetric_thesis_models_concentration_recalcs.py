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
        d = d * (density_metal / density_silicate)
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

density_metal = 7800
density_silicate = 3580
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
w_at_wt = 183.84



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

vesta_mantle_ppb = 6.67 * (10**-9)
earth_mantle_ppb = vesta_mantle_ppb
vesta_magma_ocean_moles_1 = (vesta_vol_mesh_1_thru_4 * density_silicate * vesta_mantle_ppb) / w_at_wt
vesta_magma_ocean_moles_2 = (vesta_vol_mesh_5_thru_8 * density_silicate * vesta_mantle_ppb) / w_at_wt
earth_magma_ocean_moles_1 = (earth_vol_mesh_1_thru_4 * density_silicate * earth_mantle_ppb) / w_at_wt
earth_magma_ocean_moles_2 = (earth_vol_mesh_5_thru_8 * density_silicate * earth_mantle_ppb) / w_at_wt

concs_mesh_vesta_1, concs_objs_vesta_1, moles_mesh_vesta_1, moles_objs_vesta_1, verify_D_vesta_1 = recalcConcentration(predicted_d=vesta_1['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_1, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_2, concs_objs_vesta_2, moles_mesh_vesta_2, moles_objs_vesta_2, verify_D_vesta_2 = recalcConcentration(predicted_d=vesta_2['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_1, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_3, concs_objs_vesta_3, moles_mesh_vesta_3, moles_objs_vesta_3, verify_D_vesta_3 = recalcConcentration(predicted_d=vesta_3['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_1, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_4, concs_objs_vesta_4, moles_mesh_vesta_4, moles_objs_vesta_4, verify_D_vesta_4 = recalcConcentration(predicted_d=vesta_4['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_1, original_moles_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=droplet_radius)
concs_mesh_vesta_5, concs_objs_vesta_5, moles_mesh_vesta_5, moles_objs_vesta_5, verify_D_vesta_5 = recalcConcentration(predicted_d=vesta_5['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_2, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_6, concs_objs_vesta_6, moles_mesh_vesta_6, moles_objs_vesta_6, verify_D_vesta_6 = recalcConcentration(predicted_d=vesta_6['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_2, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_7, concs_objs_vesta_7, moles_mesh_vesta_7, moles_objs_vesta_7, verify_D_vesta_7 = recalcConcentration(predicted_d=vesta_7['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_2, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)
concs_mesh_vesta_8, concs_objs_vesta_8, moles_mesh_vesta_8, moles_objs_vesta_8, verify_D_vesta_8 = recalcConcentration(predicted_d=vesta_8['D'],
                              original_moles_silicate=vesta_magma_ocean_moles_2, original_moles_metal=0, volume_mesh=vesta_vol_mesh_5_thru_8, radius_object=droplet_radius)

# reverse_D_vesta_1 = forIterReverseD(obj_concs=concs_objs_vesta_1, cell_concs=concs_mesh_vesta_1)
# reverse_D_vesta_2 = forIterReverseD(obj_concs=concs_objs_vesta_2, cell_concs=concs_mesh_vesta_2)
# reverse_D_vesta_3 = forIterReverseD(obj_concs=concs_objs_vesta_3, cell_concs=concs_mesh_vesta_3)
# reverse_D_vesta_4 = forIterReverseD(obj_concs=concs_objs_vesta_4, cell_concs=concs_mesh_vesta_4)
# reverse_D_vesta_5 = forIterReverseD(obj_concs=concs_objs_vesta_5, cell_concs=concs_mesh_vesta_5)
# reverse_D_vesta_6 = forIterReverseD(obj_concs=concs_objs_vesta_6, cell_concs=concs_mesh_vesta_6)
# reverse_D_vesta_7 = forIterReverseD(obj_concs=concs_objs_vesta_7, cell_concs=concs_mesh_vesta_7)
# reverse_D_vesta_8 = forIterReverseD(obj_concs=concs_objs_vesta_8, cell_concs=concs_mesh_vesta_8)


num_droplets_vesta = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=droplet_radius)
num_droplets_earth = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=earth_droplet_radius)


# earth_1 = pd.read_csv("thesis_model_outputs/Earth_1.csv")
# earth_2 = pd.read_csv("thesis_model_outputs/Earth_2.csv")
# earth_3 = pd.read_csv("thesis_model_outputs/Earth_3.csv")
# earth_4 = pd.read_csv("thesis_model_outputs/Earth_4.csv")
# earth_5 = pd.read_csv("thesis_model_outputs/Earth_5.csv")
# earth_6 = pd.read_csv("thesis_model_outputs/Earth_6.csv")
# earth_7 = pd.read_csv("thesis_model_outputs/Earth_7.csv")
# earth_8 = pd.read_csv("thesis_model_outputs/Earth_8.csv")
#
# depth_earth_1 = [i / 1000 for i in [0] + list(earth_1['z-depth'])]
# depth_earth_2 = [i / 1000 for i in [0] + list(earth_2['z-depth'])]
# depth_earth_3 = [i / 1000 for i in [0] + list(earth_3['z-depth'])]
# depth_earth_4 = [i / 1000 for i in [0] + list(earth_4['z-depth'])]
# depth_earth_5 = [i / 1000 for i in [0] + list(earth_5['z-depth'])]
# depth_earth_6 = [i / 1000 for i in [0] + list(earth_6['z-depth'])]
# depth_earth_7 = [i / 1000 for i in [0] + list(earth_7['z-depth'])]
# depth_earth_8 = [i / 1000 for i in [0] + list(earth_8['z-depth'])]
#
# concs_mesh_earth_1, concs_objs_earth_1, moles_mesh_earth_1, moles_objs_earth_1, verify_D_earth_1 = recalcConcentration(predicted_d=earth_1['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_1, original_moles_metal=moles_per_droplet_earth_1, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
# concs_mesh_earth_2, concs_objs_earth_2, moles_mesh_earth_2, moles_objs_earth_2, verify_D_earth_2 = recalcConcentration(predicted_d=earth_2['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_1, original_moles_metal=moles_per_droplet_earth_2, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
# concs_mesh_earth_3, concs_objs_earth_3, moles_mesh_earth_3, moles_objs_earth_3, verify_D_earth_3 = recalcConcentration(predicted_d=earth_3['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_1, original_moles_metal=moles_per_droplet_earth_3, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
# concs_mesh_earth_4, concs_objs_earth_4, moles_mesh_earth_4, moles_objs_earth_4, verify_D_earth_4 = recalcConcentration(predicted_d=earth_4['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_1, original_moles_metal=moles_per_droplet_earth_4, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=earth_droplet_radius)
# concs_mesh_earth_5, concs_objs_earth_5, moles_mesh_earth_5, moles_objs_earth_5, verify_D_earth_5 = recalcConcentration(predicted_d=earth_5['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_2, original_moles_metal=moles_per_droplet_earth_5, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
# concs_mesh_earth_6, concs_objs_earth_6, moles_mesh_earth_6, moles_objs_earth_6, verify_D_earth_6 = recalcConcentration(predicted_d=earth_6['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_2, original_moles_metal=moles_per_droplet_earth_6, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
# concs_mesh_earth_7, concs_objs_earth_7, moles_mesh_earth_7, moles_objs_earth_7, verify_D_earth_7 = recalcConcentration(predicted_d=earth_7['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_2, original_moles_metal=moles_per_droplet_earth_7, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)
# concs_mesh_earth_8, concs_objs_earth_8, moles_mesh_earth_8, moles_objs_earth_8, verify_D_earth_8 = recalcConcentration(predicted_d=earth_8['D'],
#                               original_moles_silicate=earth_magma_ocean_moles_2, original_moles_metal=moles_per_droplet_earth_8, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=earth_droplet_radius)



fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(depth_vesta_1, concs_objs_vesta_1, linewidth=2.0)
ax1.grid()





plt.show()