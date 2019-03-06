import os
from math import sqrt, pi
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

def surfaceGravity(radius, melt_density):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius**3 * melt_density) / radius**2
    return gSurface

def gravity(radius, melt_density, depth):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius ** 3 * melt_density) / radius**2
    gAccel = gSurface * (1 - (depth/radius))
    return gAccel

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


radius = 265 * 1000
melt_density = 3750
depth_increment = 5 * 1000 # km
surface_gravity = surfaceGravity(radius=radius, melt_density=melt_density)

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


print(hydrostat[-1], hydrostat[0])
print((hydrostat[-1] - hydrostat[0]) / hydrostat[0])
# plot it!
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(hydrostatic_depths, hydrostat)
ax.set_xlabel("Depth (km)")
ax.set_ylabel("Pressure (gPa)")
ax.set_title("Hydrostatic Pressure Gradient on Vesta")
ax.grid()
plt.show()
