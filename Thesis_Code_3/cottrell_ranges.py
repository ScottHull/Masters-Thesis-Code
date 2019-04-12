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

    if 0 <= pressure <= 3:
        coeffs['alpha'] = 1.11
        coeffs['beta'] = -1.18
        coeffs['chi'] = -0.85
        coeffs['delta'] = 1680
        coeffs['epsilon'] = 487

    elif 3 < pressure:
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

earth_z_eq_eta_35_actual = 23.93 / 1000
earth_z_eq_eta_10_actual = 62.97 / 1000

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


print(vesta_adiabatic[-1], earth_adiabatic[-1], vesta_hydrostat[-1], earth_hydrostat[-1])


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

fig4 = plt.figure()
ax4_0 = fig4.add_subplot(111)
ax4_1 = fig4.add_subplot(211)
ax4_2 = fig4.add_subplot(212)
ax4_1.plot(vesta_hydrostatic_depths, vesta_hydrostat, linewidth=2.0, color='black', label='Vesta')
ax4_2.plot(earth_hydrostatic_depths, earth_hydrostat, linewidth=2.0, color='black', label='Earth')
# Turn off axis lines and ticks of the big subplot
ax4_0.spines['top'].set_color('none')
ax4_0.spines['bottom'].set_color('none')
ax4_0.spines['left'].set_color('none')
ax4_0.spines['right'].set_color('none')
ax4_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax4_0.xaxis.labelpad = 20
ax4_0.yaxis.labelpad = 20
ax4_0.set_xlabel("Depth (km)")
ax4_0.set_ylabel("Pressure (GPa)")
ax4_0.set_title("Hydrostatic Pressure Profile for Vestian and Earth Magma Ocean")
ax4_1.grid()
ax4_2.grid()
ax4_1.legend(loc='upper left')
ax4_2.legend(loc='upper left')
ax4_1.set_xlim(0, max(vesta_hydrostatic_depths))
ax4_2.set_xlim(0, max(earth_hydrostatic_depths))

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
ax3_1.fill_between(earth_hydrostatic_depths, cottrell_model_earth[0], cottrell_model_earth[1], color='red', alpha=0.5,
                   label="Earth Oxidizing Model")
ax3_1.axvline(earth_z_eq_eta_10_actual, linewidth=2.0, linestyle="--", color='black', label="z$_{eq}$ ($\eta$=10$^{-1.0}$)")
ax3_1.axvline(earth_z_eq_eta_35_actual, linewidth=2.0, linestyle="--", color='red', label="z$_{eq}$ ($\eta$=10$^{-3.5}$)")
ax3_2.axvline(earth_z_eq_eta_10_actual, linewidth=2.0, linestyle="--", color='black', label="z$_{eq}$ ($\eta$=10$^{-1.0}$)")
ax3_2.axvline(earth_z_eq_eta_35_actual, linewidth=2.0, linestyle="--", color='red', label="z$_{eq}$ ($\eta$=10$^{-3.5}$)")
ax3_1.grid()
ax3_1.legend(loc='upper right')
ax3_2.plot(earth_hydrostatic_depths, cottrell_model_earth[2], color='black')
ax3_2.plot(earth_hydrostatic_depths, cottrell_model_earth[3], color='black')
ax3_2.fill_between(earth_hydrostatic_depths, cottrell_model_earth[2], cottrell_model_earth[3], color='red', alpha=0.5,
                   label="Earth Reducing Model")
ax3_2.grid()
ax3_2.legend(loc='upper right')








plt.show()

