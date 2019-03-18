from math import pi, log, sqrt
import pandas as pd
import numpy as np
from math import pi, log, exp
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 16})


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



def adiabat(current_depth, current_temperature, depth_increment, depth, gravity, thermal_expansion, heat_capacity,
            depth_gradient, adiabat_gradient):
    current_depth += depth_increment
    # dT/dh = (alpha * g * T) / c_p
    current_temperature += ((thermal_expansion * gravity * current_temperature) / heat_capacity) * depth_increment
    adiabat_gradient.append(current_temperature)
    depth_gradient.append(current_depth / 1000)
    if current_depth < depth:  # recursion!
        return adiabat(current_depth, current_temperature, depth_increment, depth, gravity, thermal_expansion,
                       heat_capacity, depth_gradient, adiabat_gradient)
    else:
        return depth_gradient, adiabat_gradient

def hydrostatic(current_depth, depth_increment, depth, gravity, density_melt, current_pressure, depth_gradient,
                pressure_gradient):
    current_depth += depth_increment
    # dP/dh = rho * g
    current_pressure += (density_melt * gravity) * depth_increment
    pressure_gradient.append(current_pressure * (10**(-9)))  # gPa
    depth_gradient.append(current_depth / 1000)
    if current_depth < depth:
        return hydrostatic(current_depth, depth_increment, depth, gravity, density_melt, current_pressure, depth_gradient,
                pressure_gradient)
    else:
        return depth_gradient, pressure_gradient


def calcNumTotalDroplets(core_radius, droplet_radius):
    core_volume = (4 / 3) * pi * (core_radius ** 3)
    droplet_volume = (4 / 3) * pi * (droplet_radius ** 3)
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
        return iterReverseD(obj_concs=obj_concs, cell_concs=cell_concs, index=(index + 1),
                            iterReverseDList=iterReverseDList)
    else:
        return iterReverseDList


def forIterReverseD(obj_concs, cell_concs):
    iterReverseDList = []
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


earth_magma_ocean_depth = 520 * 1000
earth_magma_ocean_depth_increment = 260 * 10
earth_surface_gravity = 9.8

vesta_magma_ocean_depth = 262.8 * 1000
vesta_magma_ocean_depth_increment = 400
vesta_surface_gravity = 0.25

silicate_density = 3750
metal_density = 7800
w_diffusivity = 10**-8
thermal_expansivity = 6 * 10**(-5)
heat_capacity = 10**3


earth_adiabatic_depths, earth_adiabatic = adiabat(
    current_depth=0 / 1000,  # begin at surface
    current_temperature=2000,  # degrees K
    depth_increment=earth_magma_ocean_depth_increment,  # 10 km interval
    depth=earth_magma_ocean_depth,  # to depth of 250 km
    gravity=earth_surface_gravity,  # m/s^2, Vesta
    thermal_expansion=thermal_expansivity,
    heat_capacity=heat_capacity,
    depth_gradient=[0],  # returns this list at index=0
    adiabat_gradient=[2000],  # returns this list at index=1
)

earth_hydrostatic_depths, earth_hydrostat = hydrostatic(
    current_depth=0 / 1000,
    depth_increment=earth_magma_ocean_depth_increment,
    depth=earth_magma_ocean_depth,
    gravity=earth_surface_gravity,
    current_pressure=0,
    density_melt=silicate_density,
    depth_gradient=[0],
    pressure_gradient=[0],
)

vesta_adiabatic_depths, vesta_adiabatic = adiabat(
    current_depth=0 / 1000,  # begin at surface
    current_temperature=2000,  # degrees K
    depth_increment=vesta_magma_ocean_depth_increment,  # 10 km interval
    depth=vesta_magma_ocean_depth,  # to depth of 250 km
    gravity=vesta_surface_gravity,  # m/s^2, Vesta
    thermal_expansion=thermal_expansivity,
    heat_capacity=heat_capacity,
    depth_gradient=[0],  # returns this list at index=0
    adiabat_gradient=[2000],  # returns this list at index=1
)

vesta_hydrostatic_depths, vesta_hydrostat = hydrostatic(
    current_depth=0 / 1000,
    depth_increment=vesta_magma_ocean_depth_increment,
    depth=vesta_magma_ocean_depth,
    gravity=vesta_surface_gravity,
    current_pressure=0,
    density_melt=silicate_density,
    depth_gradient=[0],
    pressure_gradient=[0],
)

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

reverse_D_vesta_1 = forIterReverseD(obj_concs=concs_objs_vesta_1, cell_concs=concs_mesh_vesta_1)
reverse_D_vesta_2 = forIterReverseD(obj_concs=concs_objs_vesta_2, cell_concs=concs_mesh_vesta_2)
reverse_D_vesta_3 = forIterReverseD(obj_concs=concs_objs_vesta_3, cell_concs=concs_mesh_vesta_3)
reverse_D_vesta_4 = forIterReverseD(obj_concs=concs_objs_vesta_4, cell_concs=concs_mesh_vesta_4)
reverse_D_vesta_5 = forIterReverseD(obj_concs=concs_objs_vesta_5, cell_concs=concs_mesh_vesta_5)
reverse_D_vesta_6 = forIterReverseD(obj_concs=concs_objs_vesta_6, cell_concs=concs_mesh_vesta_6)
reverse_D_vesta_7 = forIterReverseD(obj_concs=concs_objs_vesta_7, cell_concs=concs_mesh_vesta_7)
reverse_D_vesta_8 = forIterReverseD(obj_concs=concs_objs_vesta_8, cell_concs=concs_mesh_vesta_8)

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

num_droplets_vesta = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=droplet_radius)
num_droplets_earth = calcNumTotalDroplets(core_radius=113*1000, droplet_radius=earth_droplet_radius)

moles_per_droplet_earth_1 = (moles_objs_vesta_1[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_2 = (moles_objs_vesta_2[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_3 = (moles_objs_vesta_3[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_4 = (moles_objs_vesta_4[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_5 = (moles_objs_vesta_5[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_6 = (moles_objs_vesta_6[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_7 = (moles_objs_vesta_7[-1] * num_droplets_vesta) / num_droplets_earth
moles_per_droplet_earth_8 = (moles_objs_vesta_8[-1] * num_droplets_vesta) / num_droplets_earth

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

reverse_D_earth_1 = forIterReverseD(obj_concs=concs_objs_earth_1, cell_concs=concs_mesh_earth_1)
reverse_D_earth_2 = forIterReverseD(obj_concs=concs_objs_earth_2, cell_concs=concs_mesh_earth_2)
reverse_D_earth_3 = forIterReverseD(obj_concs=concs_objs_earth_3, cell_concs=concs_mesh_earth_3)
reverse_D_earth_4 = forIterReverseD(obj_concs=concs_objs_earth_4, cell_concs=concs_mesh_earth_4)
reverse_D_earth_5 = forIterReverseD(obj_concs=concs_objs_earth_5, cell_concs=concs_mesh_earth_5)
reverse_D_earth_6 = forIterReverseD(obj_concs=concs_objs_earth_6, cell_concs=concs_mesh_earth_6)
reverse_D_earth_7 = forIterReverseD(obj_concs=concs_objs_earth_7, cell_concs=concs_mesh_earth_7)
reverse_D_earth_8 = forIterReverseD(obj_concs=concs_objs_earth_8, cell_concs=concs_mesh_earth_8)

cottrell_model_vesta = []
cottrell_model_earth = []
fO2_list = [-0.8, -1.10, -2.25, -2.45]
for i in fO2_list:
    temp_d_vesta = []
    temp_d_earth = []
    for index, t in enumerate(vesta_hydrostatic_depths):
        pressure = vesta_hydrostat[index]
        temperature = vesta_adiabatic[index]
        d = partition(pressure=pressure, temperature=temperature, deltaIW=i)
        temp_d_vesta.append(d)
    for index, t in enumerate(earth_hydrostatic_depths):
        pressure = earth_hydrostat[index]
        temperature = earth_adiabatic[index]
        d = partition(pressure=pressure, temperature=temperature, deltaIW=i)
        temp_d_earth.append(d)
    cottrell_model_vesta.append(temp_d_vesta)
    cottrell_model_earth.append(temp_d_earth)


# adiabatic/hydrostatic gradient figure
fig1 = plt.figure()
ax1_0 = fig1.add_subplot(111)
ax1_1 = fig1.add_subplot(211)
ax1_2 = fig1.add_subplot(212)
ax1_1.plot(vesta_adiabatic_depths, vesta_adiabatic, linewidth=2.0, color='black', label='Vesta')
ax1_2.plot(earth_adiabatic_depths, earth_adiabatic, linewidth=2.0, color='black', label='Earth')
# Turn off axis lines and ticks of the big subplot
ax1_0.spines['top'].set_color('none')
ax1_0.spines['bottom'].set_color('none')
ax1_0.spines['left'].set_color('none')
ax1_0.spines['right'].set_color('none')
ax1_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax1_0.xaxis.labelpad = 20
ax1_0.yaxis.labelpad = 20
ax1_0.set_xlabel("Depth (km)")
ax1_0.set_ylabel("Temperature ($\degree$K)")
ax1_0.set_title("Adiabatic Temperature Profile for Vestian and Earth Magma Ocean")
ax1_1.grid()
ax1_2.grid()
ax1_1.legend(loc='upper left')
ax1_2.legend(loc='upper left')
ax1_1.set_xlim(0, max(vesta_adiabatic_depths))
ax1_2.set_xlim(0, max(earth_adiabatic_depths))

# cottrell partition models for Vesta
fig2 = plt.figure()
ax2_0 = fig2.add_subplot(111)
ax2_1 = fig2.add_subplot(211)
ax2_2 = fig2.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax2_0.spines['top'].set_color('none')
ax2_0.spines['bottom'].set_color('none')
ax2_0.spines['left'].set_color('none')
ax2_0.spines['right'].set_color('none')
ax2_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax2_0.xaxis.labelpad = 20
ax2_0.yaxis.labelpad = 20
ax2_0.set_xlabel("Depth (km)")
ax2_0.set_ylabel("D")
ax2_0.set_title("Predicted Partition Coefficients for a Vestian Magma Ocean")
ax2_1.plot(vesta_hydrostatic_depths, cottrell_model_vesta[0], color='black')
ax2_1.plot(vesta_hydrostatic_depths, cottrell_model_vesta[1], color='black')
ax2_1.fill_between(vesta_hydrostatic_depths, cottrell_model_vesta[0], cottrell_model_vesta[1], color='red', alpha=0.4,
                   label="Vesta Oxidizing Model")
ax2_1.grid()
ax2_1.legend(loc='upper right')
ax2_2.plot(vesta_hydrostatic_depths, cottrell_model_vesta[2], color='black')
ax2_2.plot(vesta_hydrostatic_depths, cottrell_model_vesta[3], color='black')
ax2_2.fill_between(vesta_hydrostatic_depths, cottrell_model_vesta[2], cottrell_model_vesta[3], color='red', alpha=0.4,
                   label="Vesta Reducing Model")
ax2_2.grid()
ax2_2.legend(loc='upper right')

# cottrell partition models for Earth
fig3 = plt.figure()
ax3_0 = fig3.add_subplot(111)
ax3_1 = fig3.add_subplot(211)
ax3_2 = fig3.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax3_0.spines['top'].set_color('none')
ax3_0.spines['bottom'].set_color('none')
ax3_0.spines['left'].set_color('none')
ax3_0.spines['right'].set_color('none')
ax3_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax3_0.xaxis.labelpad = 20
ax3_0.yaxis.labelpad = 20
ax3_0.set_xlabel("Depth (km)")
ax3_0.set_ylabel("D")
ax3_0.set_title("Predicted Partition Coefficients for an Earth Magma Ocean")
ax3_1.plot(earth_hydrostatic_depths, cottrell_model_earth[0], color='black')
ax3_1.plot(earth_hydrostatic_depths, cottrell_model_earth[1], color='black')
ax3_1.fill_between(earth_hydrostatic_depths, cottrell_model_earth[0], cottrell_model_earth[1], color='red', alpha=0.4,
                   label="Earth Oxidizing Model")
ax3_1.grid()
ax3_1.legend(loc='upper right')
ax3_2.plot(earth_hydrostatic_depths, cottrell_model_earth[2], color='black')
ax3_2.plot(earth_hydrostatic_depths, cottrell_model_earth[3], color='black')
ax3_2.fill_between(earth_hydrostatic_depths, cottrell_model_earth[2], cottrell_model_earth[3], color='red', alpha=0.4,
                   label="Earth Reducing Model")
ax3_2.grid()
ax3_2.legend(loc='upper right')

# eta = 10^-3.5 model mapped to earth
fig4 = plt.figure()
ax4_0 = fig4.add_subplot(111)
ax4_1 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax4_0.spines['top'].set_color('none')
ax4_0.spines['bottom'].set_color('none')
ax4_0.spines['left'].set_color('none')
ax4_0.spines['right'].set_color('none')
ax4_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax4_0.xaxis.labelpad = 20
ax4_0.yaxis.labelpad = 20
ax4_0.set_xlabel("Depth (km)")
ax4_0.set_ylabel("D")
ax4_0.set_title("Predicted Partition Coefficients for an Earth Magma Ocean ($\eta$=10$^{-3.5}$)")
ax4_1.plot(earth_hydrostatic_depths, cottrell_model_earth[0], color='black')
ax4_1.plot(earth_hydrostatic_depths, cottrell_model_earth[1], color='black')
ax4_1.fill_between(earth_hydrostatic_depths, cottrell_model_earth[0], cottrell_model_earth[1], color='red', alpha=0.4,
                   label="Earth Oxidizing Model")
ax4_1.fill_between(earth_hydrostatic_depths, [list(vesta_1['D'])[-1] for i in earth_hydrostatic_depths],
                   [list(vesta_2['D'])[-1] for i in earth_hydrostatic_depths], color='blue', label='Final D on Vesta')
ax4_1.fill_between(earth_hydrostatic_depths, [verify_D_earth_1[0] for i in earth_hydrostatic_depths], [verify_D_earth_2[0] for i in earth_hydrostatic_depths], alpha=0.2, color='green', label='Initial D on Earth')
ax4_1.grid()
ax4_1.legend(loc='upper right')
ax4_2.plot(earth_hydrostatic_depths, cottrell_model_earth[2], color='black')
ax4_2.plot(earth_hydrostatic_depths, cottrell_model_earth[3], color='black')
ax4_2.fill_between(earth_hydrostatic_depths, cottrell_model_earth[2], cottrell_model_earth[3], color='red', alpha=0.4,
                   label="Earth Reducing Model")
ax4_2.fill_between(earth_hydrostatic_depths, [list(vesta_3['D'])[-1] for i in earth_hydrostatic_depths],
                   [list(vesta_4['D'])[-1] for i in earth_hydrostatic_depths], color='blue', label='Final D on Vesta')
ax4_2.fill_between(earth_hydrostatic_depths, [verify_D_earth_3[0] for i in earth_hydrostatic_depths], [verify_D_earth_4[0] for i in earth_hydrostatic_depths], alpha=0.2, color='green', label='Initial D on Earth')
ax4_2.grid()
ax4_2.legend(loc='upper right')

fig5 = plt.figure()
ax5_0 = fig5.add_subplot(111)
ax5_1 = fig5.add_subplot(211)
ax5_2 = fig5.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax5_0.spines['top'].set_color('none')
ax5_0.spines['bottom'].set_color('none')
ax5_0.spines['left'].set_color('none')
ax5_0.spines['right'].set_color('none')
ax5_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax5_0.xaxis.labelpad = 20
ax5_0.yaxis.labelpad = 20
ax5_0.set_xlabel("Depth (km)")
ax5_0.set_ylabel("D")
ax5_0.set_title("Predicted Partition Coefficients for an Earth Magma Ocean ($\eta$=10$^{-1.0}$)")
ax5_1.plot(earth_hydrostatic_depths, cottrell_model_earth[0], color='black')
ax5_1.plot(earth_hydrostatic_depths, cottrell_model_earth[1], color='black')
ax5_1.fill_between(earth_hydrostatic_depths, cottrell_model_earth[0], cottrell_model_earth[1], color='red', alpha=0.4,
                   label="Earth Oxidizing Model")
ax5_1.fill_between(earth_hydrostatic_depths, [list(vesta_1['D'])[-1] for i in earth_hydrostatic_depths],
                   [list(vesta_2['D'])[-1] for i in earth_hydrostatic_depths], color='blue', label='Final D on Vesta')
ax5_1.fill_between(earth_hydrostatic_depths, [verify_D_earth_5[0] for i in earth_hydrostatic_depths], [verify_D_earth_6[0] for i in earth_hydrostatic_depths], alpha=0.2, color='green', label='Initial D on Earth')
ax5_1.grid()
ax5_1.legend(loc='upper right')
ax5_2.plot(earth_hydrostatic_depths, cottrell_model_earth[2], color='black')
ax5_2.plot(earth_hydrostatic_depths, cottrell_model_earth[3], color='black')
ax5_2.fill_between(earth_hydrostatic_depths, cottrell_model_earth[2], cottrell_model_earth[3], color='red', alpha=0.4,
                   label="Earth Reducing Model")
ax5_2.fill_between(earth_hydrostatic_depths, [list(vesta_3['D'])[-1] for i in earth_hydrostatic_depths],
                   [list(vesta_4['D'])[-1] for i in earth_hydrostatic_depths], color='blue', label='Final D on Vesta')
ax5_2.fill_between(earth_hydrostatic_depths, [verify_D_earth_7[0] for i in earth_hydrostatic_depths], [verify_D_earth_8[0] for i in earth_hydrostatic_depths], alpha=0.2, color='green', label='Initial D on Earth')
ax5_2.grid()
ax5_2.legend(loc='upper right')

plt.show()

