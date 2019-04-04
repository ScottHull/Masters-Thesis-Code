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

def recalcConcentration(predicted_d, original_mass_conc_silicate, original_mass_conc_metal, volume_mesh, radius_object,
                        density_metal, density_melt, mass_droplet, mass_silicate, w_at_wt):
    volume_obj = (4 / 3) * pi * (radius_object ** 3)
    original_moles_silicate = (original_mass_conc_silicate * (mass_silicate)) / (w_at_wt / 1000)
    original_moles_metal = (original_mass_conc_metal * (mass_droplet)) / (w_at_wt / 1000)

    original_conc_silicate = original_moles_silicate / volume_mesh
    original_conc_metal = original_moles_metal / volume_obj
    concs_mesh = [original_conc_silicate]
    concs_objs = [original_conc_metal]
    moles_mesh = [original_moles_silicate]
    moles_objs = [original_moles_metal]
    verify_D = [concs_objs[0] / concs_mesh[0]]
    mass_d_list = []
    volume_d_list = []
    for index, d_mole in enumerate(list(predicted_d)):
        d_volume = (density_metal / density_melt) * d_mole
        old_moles_obj = moles_objs[index - 1]
        old_moles_cell = moles_mesh[index - 1]

        adj_matrix = (old_moles_cell /
                      (1 + (3 * volume_mesh * (
                                  (4 * pi * (radius_object ** 3) * d_volume) ** (-1)))))
        adj_object = (old_moles_obj /
                      (1 + (4 * pi * (radius_object ** 3) * d_volume) * (
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
        mass_d_list.append(d_mole)
        volume_d_list.append(d_volume)

    return concs_mesh, concs_objs, moles_mesh, moles_objs, mass_d_list, volume_d_list



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



density_metal = 7800
density_silicate = 3580
vesta_droplet_radius = 0.0185
earth_droplet_radius = 0.002957762
vesta_droplet_volume = (4/3) * pi * (vesta_droplet_radius**3)
earth_droplet_volume = (4/3) * pi * (earth_droplet_radius**3)
vesta_droplet_mass = density_metal * vesta_droplet_volume
earth_droplet_mass = density_metal * earth_droplet_volume
diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.2580697580112788, droplet_radius=vesta_droplet_radius)
earth_diff_length = calcDiffusionLength(chem_diffusivity=10**-8, settling_velocity=0.646064514, droplet_radius=earth_droplet_radius)
vesta_z_eq_1_thru_4 = 240
vesta_z_eq_5_thru_8 = 600
vesta_vol_mesh_1_thru_4 = ((2 * (vesta_droplet_radius + diff_length))**2) * vesta_z_eq_1_thru_4
vesta_vol_mesh_5_thru_8 = ((2 * (vesta_droplet_radius + diff_length))**2) * vesta_z_eq_5_thru_8
vesta_mass_mesh_1_thru_4 = density_silicate * vesta_vol_mesh_1_thru_4
vesta_mass_mesh_5_thru_8 = density_silicate * vesta_vol_mesh_5_thru_8
earth_z_eq_1_thru_4 = 25
earth_z_eq_5_thru_8 = 65
earth_vol_mesh_1_thru_4 = ((2 * (earth_droplet_radius + earth_diff_length))**2) * earth_z_eq_1_thru_4
earth_vol_mesh_5_thru_8 = ((2 * (earth_droplet_radius + earth_diff_length))**2) * earth_z_eq_5_thru_8
w_at_wt = 183.84

modeled_mass_w182_in_vesta_core_08 = 10901550702.125
modeled_mass_w182_in_vesta_core_110 = 11973303629.4755
modeled_mass_w182_in_vesta_core_225 = 12938743437.6168
modeled_mass_w182_in_vesta_core_245 = 12959044131.1024



magma_ocean_concs_182w = pd.read_csv("182w_mantle.csv")
magma_ocean_concs_184w = pd.read_csv("184w_mantle.csv")

w182_fO2_08_time1 = list(magma_ocean_concs_182w['fO2_-0.8'])[1] * (10**-9)
w182_fO2_08_time5 = list(magma_ocean_concs_182w['fO2_-0.8'])[-1] * (10**-9)
w182_fO2_110_time1 = list(magma_ocean_concs_182w['fO2_-1.1'])[1] * (10**-9)
w182_fO2_110_time5 = list(magma_ocean_concs_182w['fO2_-1.1'])[-1] * (10**-9)
w182_fO2_225_time1 = list(magma_ocean_concs_182w['fO2_-2.25'])[1] * (10**-9)
w182_fO2_225_time5 = list(magma_ocean_concs_182w['fO2_-2.25'])[-1] * (10**-9)
w182_fO2_245_time1 = list(magma_ocean_concs_182w['fO2_-0.8'])[1] * (10**-9)
w182_fO2_245_time5 = list(magma_ocean_concs_182w['fO2_-2.45'])[-1] * (10**-9)


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



concs_mesh_vesta_1_time1, concs_objs_vesta_1_time1, moles_mesh_vesta_1_time1, moles_objs_vesta_1_time1, mass_d_list_vesta_1_time1, volume_d_list_vesta_1_time1 = recalcConcentration(predicted_d=vesta_1['D'],
                              original_mass_conc_silicate=w182_fO2_08_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)
concs_mesh_vesta_1_time5, concs_objs_vesta_1_time5, moles_mesh_vesta_1_time5, moles_objs_vesta_1_time5, mass_d_list_vesta_1_time5, volume_d_list_vesta_1_time5 = recalcConcentration(predicted_d=vesta_1['D'],
                              original_mass_conc_silicate=w182_fO2_08_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)

concs_mesh_vesta_2_time1, concs_objs_vesta_2_time1, moles_mesh_vesta_2_time1, moles_objs_vesta_2_time1, mass_d_list_vesta_2_time1, volume_d_list_vesta_2_time1 = recalcConcentration(predicted_d=vesta_2['D'],
                              original_mass_conc_silicate=w182_fO2_110_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)
concs_mesh_vesta_2_time5, concs_objs_vesta_2_time5, moles_mesh_vesta_2_time5, moles_objs_vesta_2_time5, mass_d_list_vesta_2_time5, volume_d_list_vesta_2_time5 = recalcConcentration(predicted_d=vesta_2['D'],
                              original_mass_conc_silicate=w182_fO2_110_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)

concs_mesh_vesta_3_time1, concs_objs_vesta_3_time1, moles_mesh_vesta_3_time1, moles_objs_vesta_3_time1, mass_d_list_vesta_3_time1, volume_d_list_vesta_3_time1 = recalcConcentration(predicted_d=vesta_3['D'],
                              original_mass_conc_silicate=w182_fO2_225_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)
concs_mesh_vesta_3_time5, concs_objs_vesta_3_time5, moles_mesh_vesta_3_time5, moles_objs_vesta_3_time5, mass_d_list_vesta_3_time5, volume_d_list_vesta_3_time5 = recalcConcentration(predicted_d=vesta_3['D'],
                              original_mass_conc_silicate=w182_fO2_225_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)

concs_mesh_vesta_4_time1, concs_objs_vesta_4_time1, moles_mesh_vesta_4_time1, moles_objs_vesta_4_time1, mass_d_list_vesta_4_time1, volume_d_list_vesta_4_time1 = recalcConcentration(predicted_d=vesta_4['D'],
                              original_mass_conc_silicate=w182_fO2_245_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)
concs_mesh_vesta_4_time5, concs_objs_vesta_4_time5, moles_mesh_vesta_4_time5, moles_objs_vesta_4_time5, mass_d_list_vesta_4_time5, volume_d_list_vesta_4_time5 = recalcConcentration(predicted_d=vesta_4['D'],
                              original_mass_conc_silicate=w182_fO2_245_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_1_thru_4, w_at_wt=w_at_wt)

concs_mesh_vesta_5_time1, concs_objs_vesta_5_time1, moles_mesh_vesta_5_time1, moles_objs_vesta_5_time1, mass_d_list_vesta_5_time1, volume_d_list_vesta_5_time1 = recalcConcentration(predicted_d=vesta_5['D'],
                              original_mass_conc_silicate=w182_fO2_08_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)
concs_mesh_vesta_5_time5, concs_objs_vesta_5_time5, moles_mesh_vesta_5_time5, moles_objs_vesta_5_time5, mass_d_list_vesta_5_time5, volume_d_list_vesta_5_time5 = recalcConcentration(predicted_d=vesta_5['D'],
                              original_mass_conc_silicate=w182_fO2_08_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)

concs_mesh_vesta_6_time1, concs_objs_vesta_6_time1, moles_mesh_vesta_6_time1, moles_objs_vesta_6_time1, mass_d_list_vesta_6_time1, volume_d_list_vesta_6_time1 = recalcConcentration(predicted_d=vesta_6['D'],
                              original_mass_conc_silicate=w182_fO2_110_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)
concs_mesh_vesta_6_time5, concs_objs_vesta_6_time5, moles_mesh_vesta_6_time5, moles_objs_vesta_6_time5, mass_d_list_vesta_6_time5, volume_d_list_vesta_6_time5 = recalcConcentration(predicted_d=vesta_6['D'],
                              original_mass_conc_silicate=w182_fO2_110_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)

concs_mesh_vesta_7_time1, concs_objs_vesta_7_time1, moles_mesh_vesta_7_time1, moles_objs_vesta_7_time1, mass_d_list_vesta_7_time1, volume_d_list_vesta_7_time1 = recalcConcentration(predicted_d=vesta_7['D'],
                              original_mass_conc_silicate=w182_fO2_225_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)
concs_mesh_vesta_7_time5, concs_objs_vesta_7_time5, moles_mesh_vesta_7_time5, moles_objs_vesta_7_time5, mass_d_list_vesta_7_time5, volume_d_list_vesta_7_time5 = recalcConcentration(predicted_d=vesta_7['D'],
                              original_mass_conc_silicate=w182_fO2_225_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)
concs_mesh_vesta_8_time1, concs_objs_vesta_8_time1, moles_mesh_vesta_8_time1, moles_objs_vesta_8_time1, mass_d_list_vesta_8_time1, volume_d_list_vesta_8_time1 = recalcConcentration(predicted_d=vesta_8['D'],
                              original_mass_conc_silicate=w182_fO2_245_time1, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)
concs_mesh_vesta_8_time5, concs_objs_vesta_8_time5, moles_mesh_vesta_8_time5, moles_objs_vesta_8_time5, mass_d_list_vesta_8_time5, volume_d_list_vesta_8_time5 = recalcConcentration(predicted_d=vesta_8['D'],
                              original_mass_conc_silicate=w182_fO2_245_time5, original_mass_conc_metal=0, volume_mesh=vesta_vol_mesh_1_thru_4, radius_object=vesta_droplet_radius, density_melt=density_silicate, density_metal=density_metal,
                                                mass_droplet=vesta_droplet_mass, mass_silicate=vesta_mass_mesh_5_thru_8, w_at_wt=w_at_wt)

mass_conc_droplet_vesta_1_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_1_time1]
mass_conc_droplet_vesta_1_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_1_time5]
mass_conc_droplet_vesta_2_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_2_time1]
mass_conc_droplet_vesta_2_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_2_time5]
mass_conc_droplet_vesta_3_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_3_time1]
mass_conc_droplet_vesta_3_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_3_time5]
mass_conc_droplet_vesta_4_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_4_time1]
mass_conc_droplet_vesta_4_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_4_time5]
mass_conc_droplet_vesta_5_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_5_time1]
mass_conc_droplet_vesta_5_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_5_time5]
mass_conc_droplet_vesta_6_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_6_time1]
mass_conc_droplet_vesta_6_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_6_time5]
mass_conc_droplet_vesta_7_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_7_time1]
mass_conc_droplet_vesta_7_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_7_time5]
mass_conc_droplet_vesta_8_time1 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_8_time1]
mass_conc_droplet_vesta_8_time5 = [(((i * w_at_wt) * vesta_droplet_volume / 1000) / vesta_droplet_mass) * (10**9) for i in concs_objs_vesta_8_time5]



num_droplets_vesta = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=vesta_droplet_radius)
num_droplets_earth = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=earth_droplet_radius)

# print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
#     moles_objs_vesta_1[-1], moles_objs_vesta_2[-1], moles_objs_vesta_3[-1],  moles_objs_vesta_4[-1],
#     moles_objs_vesta_5[-1], moles_objs_vesta_6[-1], moles_objs_vesta_7[-1], moles_objs_vesta_8[-1]
# ))


moles_per_droplet_earth_1_time1 = (moles_objs_vesta_1_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_1_time5 = (moles_objs_vesta_1_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_2_time1 = (moles_objs_vesta_2_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_2_time5 = (moles_objs_vesta_2_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_3_time1 = (moles_objs_vesta_3_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_3_time5 = (moles_objs_vesta_3_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_4_time1 = (moles_objs_vesta_4_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_4_time5 = (moles_objs_vesta_4_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_5_time1 = (moles_objs_vesta_5_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_5_time5 = (moles_objs_vesta_5_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_6_time1 = (moles_objs_vesta_6_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_6_time5 = (moles_objs_vesta_6_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_7_time1 = (moles_objs_vesta_7_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_7_time5 = (moles_objs_vesta_7_time5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_8_time1 = (moles_objs_vesta_8_time1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_8_time5 = (moles_objs_vesta_8_time5[-1] * num_droplets_vesta) / num_droplets_earth


# print('1: {}\n2: {}\n3: {}\n4: {}\n5: {}\n6: {}\n7: {}\n8: {}\n'.format(
#     moles_per_droplet_earth_1, moles_per_droplet_earth_2, moles_per_droplet_earth_3,  moles_per_droplet_earth_4,
#     moles_per_droplet_earth_5, moles_per_droplet_earth_6, moles_per_droplet_earth_7, moles_per_droplet_earth_8
# ))
#
# print("Num Droplets Vesta: {}\nNum Droplets Earth: {}".format(num_droplets_vesta, num_droplets_earth))


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

# concs_mesh_earth_1, concs_objs_earth_1, moles_mesh_earth_1, moles_objs_earth_1, mass_d_list_earth_1, volume_d_list_earth_1  = recalcConcentration(predicted_d=earth_1['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_1, original_mass_conc_metal=moles_per_droplet_earth_1, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_2, concs_objs_earth_2, moles_mesh_earth_2, moles_objs_earth_2,  mass_d_list_earth_2, volume_d_list_earth_2 = recalcConcentration(predicted_d=earth_2['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_1, original_mass_conc_metal=moles_per_droplet_earth_2, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_3, concs_objs_earth_3, moles_mesh_earth_3, moles_objs_earth_3, mass_d_list_earth_3, volume_d_list_earth_3 = recalcConcentration(predicted_d=earth_3['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_1, original_mass_conc_metal=moles_per_droplet_earth_3, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_4, concs_objs_earth_4, moles_mesh_earth_4, moles_objs_earth_4, mass_d_list_earth_4, volume_d_list_earth_4 = recalcConcentration(predicted_d=earth_4['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_1, original_mass_conc_metal=moles_per_droplet_earth_4, volume_mesh=earth_vol_mesh_1_thru_4, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_5, concs_objs_earth_5, moles_mesh_earth_5, moles_objs_earth_5, mass_d_list_earth_5, volume_d_list_earth_5 = recalcConcentration(predicted_d=earth_5['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_2, original_mass_conc_metal=moles_per_droplet_earth_5, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_6, concs_objs_earth_6, moles_mesh_earth_6, moles_objs_earth_6,  mass_d_list_earth_6, volume_d_list_earth_6 = recalcConcentration(predicted_d=earth_6['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_2, original_mass_conc_metal=moles_per_droplet_earth_6, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_7, concs_objs_earth_7, moles_mesh_earth_7, moles_objs_earth_7, mass_d_list_earth_7, volume_d_list_earth_7 = recalcConcentration(predicted_d=earth_7['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_2, original_mass_conc_metal=moles_per_droplet_earth_7, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)
# concs_mesh_earth_8, concs_objs_earth_8, moles_mesh_earth_8, moles_objs_earth_8, mass_d_list_earth_8, volume_d_list_earth_8 = recalcConcentration(predicted_d=earth_8['D'],
#                               original_mass_conc_silicate=earth_magma_ocean_moles_2, original_mass_conc_metal=moles_per_droplet_earth_8, volume_mesh=earth_vol_mesh_5_thru_8, radius_object=droplet_radius,
#                               density_melt=density_silicate, density_metal=density_metal)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(depth_vesta_1, mass_conc_droplet_vesta_1_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax1.plot(depth_vesta_1, mass_conc_droplet_vesta_1_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax1.fill_between(depth_vesta_1, mass_conc_droplet_vesta_1_time1, mass_conc_droplet_vesta_1_time5, color='red', alpha=0.2)
ax1.grid()
ax1.set_xlabel("Depth (km)")
ax1.set_ylabel("Concentration (ppb)")
ax1.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 1)")
ax1.legend(loc='upper left')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(depth_vesta_2, mass_conc_droplet_vesta_2_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax2.plot(depth_vesta_2, mass_conc_droplet_vesta_2_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax2.fill_between(depth_vesta_2, mass_conc_droplet_vesta_2_time1, mass_conc_droplet_vesta_2_time5, color='red', alpha=0.2)
ax2.grid()
ax2.set_xlabel("Depth (km)")
ax2.set_ylabel("Concentration (ppb)")
ax2.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 2)")
ax2.legend(loc='upper left')

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(depth_vesta_3, mass_conc_droplet_vesta_3_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax3.plot(depth_vesta_3, mass_conc_droplet_vesta_3_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax3.fill_between(depth_vesta_3, mass_conc_droplet_vesta_3_time1, mass_conc_droplet_vesta_3_time5, color='red', alpha=0.2)
ax3.grid()
ax3.set_xlabel("Depth (km)")
ax3.set_ylabel("Concentration (ppb)")
ax3.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 3)")
ax3.legend(loc='upper left')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(depth_vesta_4, mass_conc_droplet_vesta_4_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax4.plot(depth_vesta_4, mass_conc_droplet_vesta_4_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax4.fill_between(depth_vesta_4, mass_conc_droplet_vesta_4_time1, mass_conc_droplet_vesta_4_time5, color='red', alpha=0.2)
ax4.grid()
ax4.set_xlabel("Depth (km)")
ax4.set_ylabel("Concentration (ppb)")
ax4.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 4)")
ax4.legend(loc='upper left')

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
ax5.plot(depth_vesta_5, mass_conc_droplet_vesta_5_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax5.plot(depth_vesta_5, mass_conc_droplet_vesta_5_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax5.fill_between(depth_vesta_5, mass_conc_droplet_vesta_5_time1, mass_conc_droplet_vesta_5_time5, color='red', alpha=0.2)
ax5.grid()
ax5.set_xlabel("Depth (km)")
ax5.set_ylabel("Concentration (ppb)")
ax5.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 5)")
ax5.legend(loc='upper left')

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot(depth_vesta_6, mass_conc_droplet_vesta_6_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax6.plot(depth_vesta_6, mass_conc_droplet_vesta_6_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax6.fill_between(depth_vesta_6, mass_conc_droplet_vesta_6_time1, mass_conc_droplet_vesta_6_time5, color='red', alpha=0.2)
ax6.grid()
ax6.set_xlabel("Depth (km)")
ax6.set_ylabel("Concentration (ppb)")
ax6.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 6)")
ax6.legend(loc='upper left')

fig7 = plt.figure()
ax7 = fig7.add_subplot(111)
ax7.plot(depth_vesta_7, mass_conc_droplet_vesta_7_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax7.plot(depth_vesta_7, mass_conc_droplet_vesta_7_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax7.fill_between(depth_vesta_7, mass_conc_droplet_vesta_7_time1, mass_conc_droplet_vesta_7_time5, color='red', alpha=0.2)
ax7.grid()
ax7.set_xlabel("Depth (km)")
ax7.set_ylabel("Concentration (ppb)")
ax7.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 7)")
ax7.legend(loc='upper left')

fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.plot(depth_vesta_8, mass_conc_droplet_vesta_8_time1, linewidth=2, linestyle="--", color='black', label='0 Ma')
ax8.plot(depth_vesta_8, mass_conc_droplet_vesta_8_time5, linewidth=2, linestyle="-", color='black', label='5 Ma')
ax8.fill_between(depth_vesta_8, mass_conc_droplet_vesta_8_time1, mass_conc_droplet_vesta_8_time5, color='red', alpha=0.2)
ax8.grid()
ax8.set_xlabel("Depth (km)")
ax8.set_ylabel("Concentration (ppb)")
ax8.set_title("Mass Concentration of $^{182}$W in an Iron Droplet (Vesta Model 8)")
ax8.legend(loc='upper left')








fig9 = plt.figure()
ax9 = fig9.add_subplot(111)
ax9.plot(depth_vesta_1, [i * (10**-9) * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_1_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax9.plot(depth_vesta_1, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_1_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax9.fill_between(depth_vesta_1, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_1_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_1_time5], color='red', alpha=0.2)
ax9.grid()
ax9.set_xlabel("Depth (km)")
ax9.set_ylabel("Mass (kg)")
ax9.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 1)")
ax9.legend(loc='upper left')

fig10 = plt.figure()
ax10 = fig10.add_subplot(111)
ax10.plot(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_2_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax10.plot(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_2_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax10.fill_between(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_2_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_2_time5], color='red', alpha=0.2)
ax10.grid()
ax10.set_xlabel("Depth (km)")
ax10.set_ylabel("Mass (kg)")
ax10.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 2)")
ax10.legend(loc='upper left')

fig11 = plt.figure()
ax11 = fig11.add_subplot(111)
ax11.plot(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_3_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax11.plot(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_3_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax11.fill_between(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_3_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_3_time5], color='red', alpha=0.2)
ax11.grid()
ax11.set_xlabel("Depth (km)")
ax11.set_ylabel("Mass (kg)")
ax11.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 3)")
ax11.legend(loc='upper left')

fig12 = plt.figure()
ax12 = fig12.add_subplot(111)
ax12.plot(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_4_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax12.plot(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_4_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax12.fill_between(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_4_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_4_time5], color='red', alpha=0.2)
ax12.grid()
ax12.set_xlabel("Depth (km)")
ax12.set_ylabel("Mass (kg)")
ax12.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 4)")
ax12.legend(loc='upper left')

fig13 = plt.figure()
ax13 = fig13.add_subplot(111)
ax13.plot(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_5_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax13.plot(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_5_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax13.fill_between(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_5_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_5_time5], color='red', alpha=0.2)
ax13.grid()
ax13.set_xlabel("Depth (km)")
ax13.set_ylabel("Mass (kg)")
ax13.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 5)")
ax13.legend(loc='upper left')

fig14 = plt.figure()
ax14 = fig14.add_subplot(111)
ax14.plot(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_6_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax14.plot(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_6_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax14.fill_between(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_6_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_6_time5], color='red', alpha=0.2)
ax14.grid()
ax14.set_xlabel("Depth (km)")
ax14.set_ylabel("Mass (kg)")
ax14.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 6)")
ax14.legend(loc='upper left')

fig15 = plt.figure()
ax15 = fig15.add_subplot(111)
ax15.plot(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_7_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax15.plot(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_7_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax15.fill_between(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_7_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_7_time5], color='red', alpha=0.2)
ax15.grid()
ax15.set_xlabel("Depth (km)")
ax15.set_ylabel("Mass (kg)")
ax15.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 7)")
ax15.legend(loc='upper left')

fig16 = plt.figure()
ax16 = fig16.add_subplot(111)
ax16.plot(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_8_time1], linewidth=2, linestyle="--", color='black', label='0 Ma')
ax16.plot(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_8_time5], linewidth=2, linestyle="-", color='black', label='5 Ma')
ax16.fill_between(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_8_time1], [i * (10**-9) * vesta_droplet_mass for i in mass_conc_droplet_vesta_8_time5], color='red', alpha=0.2)
ax16.grid()
ax16.set_xlabel("Depth (km)")
ax16.set_ylabel("Mass (kg)")
ax16.set_title("Mass of $^{182}$W in an Iron Droplet (Vesta Model 8)")
ax16.legend(loc='upper left')






fig17 = plt.figure()
ax17 = fig17.add_subplot(111)
ax17.plot(depth_vesta_1, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_1_time1], label="Time = 0.25 Ma")
ax17.plot(depth_vesta_1, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_1_time5], label="Time = 5.0 Ma")
ax17.fill_between(depth_vesta_1, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_1_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_1_time5], color='red', alpha=0.2)
ax17.axhline(modeled_mass_w182_in_vesta_core_08, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax17.grid()
ax17.set_xlabel("Depth (km)")
ax17.set_ylabel("Mass (kg)")
ax17.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 1)".format(Decimal(str(num_droplets_vesta)))))
ax17.legend(loc='center right')

fig18 = plt.figure()
ax18 = fig18.add_subplot(111)
ax18.plot(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_2_time1], label="Time = 0.25 Ma")
ax18.plot(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_2_time5], label="Time = 5.0 Ma")
ax18.fill_between(depth_vesta_2, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_2_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_2_time5], color='red', alpha=0.2)
ax18.axhline(modeled_mass_w182_in_vesta_core_110, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax18.grid()
ax18.set_xlabel("Depth (km)")
ax18.set_ylabel("Mass (kg)")
ax18.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 2)".format(Decimal(str(num_droplets_vesta)))))
ax18.legend(loc='center right')

fig19 = plt.figure()
ax19 = fig19.add_subplot(111)
ax19.plot(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_3_time1], label="Time = 0.25 Ma")
ax19.plot(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_3_time5], label="Time = 5.0 Ma")
ax19.fill_between(depth_vesta_3, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_3_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_3_time5], color='red', alpha=0.2)
ax19.axhline(modeled_mass_w182_in_vesta_core_225, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax19.grid()
ax19.set_xlabel("Depth (km)")
ax19.set_ylabel("Mass (kg)")
ax19.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 3)".format(Decimal(str(num_droplets_vesta)))))
ax19.legend(loc='center right')

fig20 = plt.figure()
ax20 = fig20.add_subplot(111)
ax20.plot(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_4_time1], label="Time = 0.25 Ma")
ax20.plot(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_4_time5], label="Time = 5.0 Ma")
ax20.fill_between(depth_vesta_4, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_4_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_4_time5], color='red', alpha=0.2)
ax20.axhline(modeled_mass_w182_in_vesta_core_245, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax20.grid()
ax20.set_xlabel("Depth (km)")
ax20.set_ylabel("Mass (kg)")
ax20.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 4)".format(Decimal(str(num_droplets_vesta)))))
ax20.legend(loc='center right')

fig21 = plt.figure()
ax21 = fig21.add_subplot(111)
ax21.plot(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_5_time1], label="Time = 0.25 Ma")
ax21.plot(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_5_time5], label="Time = 5.0 Ma")
ax21.fill_between(depth_vesta_5, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_5_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_5_time5], color='red', alpha=0.2)
ax21.axhline(modeled_mass_w182_in_vesta_core_08, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax21.grid()
ax21.set_xlabel("Depth (km)")
ax21.set_ylabel("Mass (kg)")
ax21.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 5)".format(Decimal(str(num_droplets_vesta)))))
ax21.legend(loc='center right')

fig22 = plt.figure()
ax22 = fig22.add_subplot(111)
ax22.plot(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_6_time1], label="Time = 0.25 Ma")
ax22.plot(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_6_time5], label="Time = 5.0 Ma")
ax22.fill_between(depth_vesta_6, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_6_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_6_time5], color='red', alpha=0.2)
ax22.axhline(modeled_mass_w182_in_vesta_core_110, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax22.grid()
ax22.set_xlabel("Depth (km)")
ax22.set_ylabel("Mass (kg)")
ax22.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 6)".format(Decimal(str(num_droplets_vesta)))))
ax22.legend(loc='center right')

fig23 = plt.figure()
ax23 = fig23.add_subplot(111)
ax23.plot(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_7_time1], label="Time = 0.25 Ma")
ax23.plot(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_7_time5], label="Time = 5.0 Ma")
ax23.fill_between(depth_vesta_7, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_7_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_7_time5], color='red', alpha=0.2)
ax23.axhline(modeled_mass_w182_in_vesta_core_225, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax23.grid()
ax23.set_xlabel("Depth (km)")
ax23.set_ylabel("Mass (kg)")
ax23.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 7)".format(Decimal(str(num_droplets_vesta)))))
ax23.legend(loc='center right')

fig24 = plt.figure()
ax24 = fig24.add_subplot(111)
ax24.plot(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_8_time1], label="Time = 0.25 Ma")
ax24.plot(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_8_time5], label="Time = 5.0 Ma")
ax24.fill_between(depth_vesta_8, [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_8_time1], [i * (10**-9) * vesta_droplet_mass * num_droplets_vesta for i in mass_conc_droplet_vesta_8_time5], color='red', alpha=0.2)
ax24.axhline(modeled_mass_w182_in_vesta_core_245, color='black', linestyle="--", label="Modeled Mass $^{182}$W in Vesta Core")
ax24.grid()
ax24.set_xlabel("Depth (km)")
ax24.set_ylabel("Mass (kg)")
ax24.set_title(("Mass of $^{182}$W " + "Transported to Vesta's Core by Iron Rain (N={:.2E}) (Vesta Model 8)".format(Decimal(str(num_droplets_vesta)))))
ax24.legend(loc='center right')
plt.show()