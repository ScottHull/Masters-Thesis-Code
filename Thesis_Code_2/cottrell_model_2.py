import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import log, exp
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

    elif 2 < pressure < 18:
        coeffs['alpha'] = 1.05
        coeffs['beta'] = -1.10
        coeffs['chi'] = -0.84
        coeffs['delta'] = 3588
        coeffs['epsilon'] = -102

    return coeffs

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


def partition(pressure, temperature, deltaIW):
    nbo_t = 2.6
    coeffs = collectCoeffsSimple(pressure=pressure, temperature=temperature)
    alpha = coeffs['alpha']
    beta = coeffs['beta']
    chi = coeffs['chi']
    delta = coeffs['delta']
    epsilon = coeffs['epsilon']
    logD = alpha + (beta * deltaIW) + (chi * nbo_t) + (delta * (1/temperature)) + (epsilon * (pressure/temperature))
    D = exp(logD)
    return D

def nearestValue(list, nearest_value):
    sort_list = sorted(list)
    for index, i in enumerate(sort_list):
        if i < nearest_value:
            pass
        else:
            return index, i




radius = 400 * 1000
melt_density = 3750
depth_increment = 100 * 1000 # km
surface_gravity = 9.8

adiabatic_depths, adiabatic = adiabat(
    current_depth=0 / 1000,  # begin at surface
    current_temperature=2000,  # degrees K
    depth_increment=depth_increment,  # 10 km interval
    depth=radius,  # to depth of 250 km
    gravity=surface_gravity,  # m/s^2, Vesta
    thermal_expansion=6 * 10**(-5),
    heat_capacity=10**3,
    depth_gradient=[0],  # returns this list at index=0
    adiabat_gradient=[2000],  # returns this list at index=1
)

hydrostatic_depths, hydrostat = hydrostatic(
    current_depth=0 / 1000,
    depth_increment=depth_increment,
    depth=radius,
    gravity=surface_gravity,
    current_pressure=0,
    density_melt=melt_density,
    depth_gradient=[0],
    pressure_gradient=[0],
)

cottrellD = []
fO2s = [-2.45, -2.25, -1.10, -0.8]

for f in fO2s:
    cottrellD_temp = []
    for index, i in enumerate(adiabatic_depths):
        d = partition(pressure=hydrostat[index], temperature=adiabatic[index], deltaIW=f)
        cottrellD_temp.append(d)
    cottrellD.append(cottrellD_temp)


fig = plt.figure()
ax = fig.add_subplot(111)
for index, i in enumerate(cottrellD):
    ax.plot(hydrostatic_depths, i, linewidth=2.0, label="fO2: IW-{}".format(fO2s[index]))
ax.fill_between(hydrostatic_depths, cottrellD[0], cottrellD[1], color='red', alpha=0.2, label='Reducing fO2 Model')
ax.fill_between(hydrostatic_depths, cottrellD[2], cottrellD[3], color='blue', alpha=0.2, label='Oxidizing fO2 Model')
ax.axvline(6.8, 0, 1, color='red', linestyle="--", label="Related Depth to Pressure at Base of Vesta")
ax2 = ax.twiny()
ax2.set_xticks(hydrostat)
ax2.xaxis.set_ticks(np.arange(hydrostat[0], hydrostat[-1], 2))
ax2.set_xlabel("Pressure (GPa)")
ax.grid()
ax.legend(loc='upper right')
ax.set_xlabel("Depth (km)")
ax.set_ylabel("Partition Coefficient (D)")
ax.set_title("W Partition Coefficients in a Shallow Earth Magma Ocean (Cottrell et al. 2009)")
# nearest_in_earth_index, nearest_in_earth_value = nearestValue(list=hydrostat, nearest_value=0.25)
# print(nearest_in_earth_index, nearest_in_earth_value)
# ax.axvspan(xmin=0, xmax=hydrostatic_depths[nearest_in_earth_index], color='orange', alpha=0.2)

plt.show()
