from math import pi, log, sqrt, log10, exp
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt



def eqTime(A, B, C_m_t, C_m_0):
    t = (1 / A) * log((C_m_t + (B/A)) / (C_m_0 + (B/A)))
    return t

def eqDistance(eqTime, velocity):
    d = velocity * eqTime
    return d

def coeffA(velocity, radius, F_m, F_s, D_ms):
    A = (velocity / (2 * radius)) * ((F_m / (F_m + (F_s / D_ms))) - 1)
    return A

def coeffB(velocity, radius, F_m, F_s, D_ms, C_0_s):
    B = (velocity / (2 * radius)) * ((F_s * C_0_s) / (F_m + (F_s / D_ms)))
    return B

def calcFm(radius, density_droplet, density_melt, D_s, velocity):
    vol_metal = (radius**3) * density_droplet / 3
    vol_silicate = (radius**2) * density_melt * sqrt((2 * D_s * radius) / velocity)
    Fm = vol_metal / (vol_metal + vol_silicate)
    return Fm

def calcFs(Fm):
    Fs = 1 - Fm
    return Fs

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

def turbulentVelocity(gravAccel, radius_droplet, density_droplet, density_melt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * (2 * radius_droplet)) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def laminarVelocity(density_droplet, density_melt, gravity, droplet_radius, dynamic_viscosity):
    v = ((density_droplet - density_melt) * gravity * ((2 * droplet_radius)**2)) / (18 * dynamic_viscosity)
    return v


def rubieCmt(A, B, time, C_m_0):
    if C_m_0 > 0:
        c_m_t = ((C_m_0 + (B / A)) * exp(A * time)) - (B / A)
        return c_m_t
    else:
        return 0



if __name__ == "__main__":

    radius_vesta = 265 * 1000  # 265 km
    current_vesta_core_radius_min = 107 * 1000
    current_vesta_core_radius_max = 113 * 1000
    depth_increment = 5 * 1000  # 5 km
    droplet_radius = 0.01  # m
    densityDroplet = 7800  # kg m3
    densityMelt = 3750  # kg m3
    c_0_m = 32700  # ppm
    c_0_s = 42  # ppm
    chemical_diffusivity = 10 ** (-8)
    thermal_expansion = 6 * 10 ** (-5)
    heat_capacity = 10 ** 3
    dynamic_viscosity = 10**(-2)
    thermal_diffusivity = 10 ** (-6)
    distrib_coeff = 28
    surface_gravity = surfaceGravity(radius=radius_vesta, melt_density=densityMelt)

    adiabatic_depths, adiabatic = adiabat(
        current_depth=0 / 1000,  # begin at surface
        current_temperature=2000,  # degrees K
        depth_increment=depth_increment,  # 10 km interval
        depth=radius_vesta,  # to depth of 250 km
        gravity=surface_gravity,  # m/s^2, Vesta
        thermal_expansion=6 * 10 ** (-5),
        heat_capacity=10 ** 3,
        depth_gradient=[0],  # returns this list at index=0
        adiabat_gradient=[2000],  # returns this list at index=1
    )

    hydrostatic_depths, hydrostat = hydrostatic(
        current_depth=0 / 1000,
        depth_increment=depth_increment,
        depth=radius_vesta,
        gravity=surface_gravity,
        current_pressure=0,
        density_melt=densityMelt,
        depth_gradient=[0],
        pressure_gradient=[0],
    )

    for index, depth in enumerate(adiabatic_depths):
        pass

