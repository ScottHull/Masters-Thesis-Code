from math import pi, log, sqrt
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt


def dropletSize(density_melt, density_droplet, settling_velosity, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = (1/2) * ((weber_number * surface_tension) / ((density_droplet - density_melt) * (settling_velosity**2)))
    return r

def surfaceGravity(radius, melt_density):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius**3 * melt_density) / radius**2
    return gSurface

def gravity(radius, melt_density, depth):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius ** 3 * melt_density) / radius**2
    gAccel = gSurface * (1 - (depth/radius))
    return gAccel

def turbulentVelocity(gravAccel, droplet_radius, density_droplet, density_melt, Cd=0.2):
    diameter = 2 * droplet_radius
    v = sqrt(((4 * gravAccel * diameter) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def laminarVelocity(density_melt, density_droplet, gravity, droplet_radius, dynamic_viscosity):
    diameter = 2 * droplet_radius
    v = ((density_droplet - density_melt) * gravity * (diameter**2)) / (18 * dynamic_viscosity)
    return v

def rFromWeberTurbulent(density_melt, density_droplet, gravity, Cd=0.2, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = sqrt(weber_number * (((16 * gravity * (density_droplet - density_melt)) / (3 * Cd * surface_tension)) *
                        ((density_droplet - density_melt) / density_melt))**(-1))
    return r

def rFromWeberLaminar(dynamic_viscosity, density_droplet, density_melt, gravity, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = ((81 * weber_number * surface_tension * (dynamic_viscosity**2)) / (8 * (gravity**2) * ((density_droplet - density_melt)**3)))**(1/5)
    return r

def WeberNumber(density_melt, density_droplet, droplet_radius, velocity, surface_tension=1):
    we = ((density_droplet - density_melt) * (2 * droplet_radius) * (velocity**2)) / surface_tension
    return we

def critLaminarViscosity(density_melt, density_droplet, gravity, droplet_radius, f=10):
    eta = ((6 * f * density_melt) / (pi * gravity * (density_melt**2) * ((2 * droplet_radius)**3) * (density_droplet - density_melt)))**(-1/2)
    return eta

def frictionCoeff(density_melt, density_droplet, dynamic_viscosity, gravity, droplet_radius):
    f = (pi / 6) * ((density_droplet - density_melt) / density_melt) * ((density_melt / dynamic_viscosity)**2) * \
        gravity * ((2 * droplet_radius)**3)
    return f

def rayleigh(density_melt, thermal_expansivity, gravity, temp_surface, temp_depth, depth, thermal_diffusivity, dynamic_viscosity):
    temp_difference = temp_depth - temp_surface
    ra = (density_melt * thermal_expansivity * gravity * temp_difference * (depth**3)) / (thermal_diffusivity * dynamic_viscosity)
    return ra

def nusselt(rayleigh, beta=(1/3)):
    nu = (0.089) * (rayleigh**beta)
    return nu

def convectionVelocity(thermal_expansion, gravity, thermal_diffusivity, temp_surface, temp_depth, nusselt):
    temperature_difference = temp_depth - temp_surface
    v_c = ((thermal_expansion * gravity * thermal_diffusivity * temperature_difference)**(1/3)) * (nusselt**(1/3))
    return v_c

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



if __name__ == "__main__":

    radius_vesta = 263 * 1000  # 265 km
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
    dynamic_viscosity = 10 ** (-2)
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

    etas = [10 ** (-5), 10 ** (-4), 10 ** (-3), 10 ** (-2), 10 ** (-1), 10 ** (0)]

    ra_list = []
    nu_list = []
    v_c_list = []
    r_stable_list = []
    v_s_list = []
    v_s_plus_v_c_list = []
    v_s_minus_v_c_list = []

    for i in etas:
        ra = rayleigh(density_melt=densityMelt, thermal_expansivity=thermal_expansion, gravity=surface_gravity,
                      temp_surface=adiabatic[0], temp_depth=adiabatic[-1], depth=radius_vesta,
                      thermal_diffusivity=thermal_diffusivity, dynamic_viscosity=i)
        nu = nusselt(rayleigh=ra)
        v_c = convectionVelocity(thermal_expansion=thermal_expansion, gravity=surface_gravity,
                                 thermal_diffusivity=thermal_diffusivity, temp_surface=adiabatic[0],
                                 temp_depth=adiabatic[-1], nusselt=nu)

        r_stable = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity)
        v_s = turbulentVelocity(gravAccel=surface_gravity, droplet_radius=r_stable, density_droplet=densityDroplet,
                                density_melt=densityMelt)

        plus = v_s + v_c
        minus = v_s - v_c

        ra_list.append(ra)
        nu_list.append(nu)
        v_c_list.append(v_c)
        r_stable_list.append(r_stable)
        v_s_list.append(v_s)
        v_s_plus_v_c_list.append(plus)
        v_s_minus_v_c_list.append(minus)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(etas, ra_list, linewidth=2.0)
    ax1_2 = ax1.twinx()
    ax1_2.plot(etas, nu_list, linewidth=2.0, linestyle="--")
    ax1.set_xlabel("Dynamic Viscosity (Pa s)")
    ax1.set_ylabel("Rayleigh Number (solid)")
    ax1_2.set_ylabel("Nusselt Number (dashed)")
    ax1.set_title("Rayleigh Number and Nusselt Number vs. Dynamic Viscosity")
    ax1.set_xscale('log')
    ax1.axvspan(10**-3.5, 10**-1, alpha=0.2, color='red')
    ax1.grid()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(etas, v_c_list, linewidth=2.0)
    ax2.set_xlabel("Dynamic Viscosity (Pa s)")
    ax2.set_ylabel("Convective Velocity (m/s)")
    ax2.set_title("Convective Velocity vs. Dynamic Viscosity")
    ax2.set_xscale('log')
    ax2.axvspan(10**-3.5, 10**-1, alpha=0.2, color='red')
    ax2.grid()

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(etas, [i * 100 for i in r_stable_list], linewidth=2.0)
    ax3_2 = ax3.twinx()
    ax3_2.plot(etas, v_s_list, linewidth=2.0, linestyle="--")
    ax3.set_xlabel("Dynamic Viscosity (Pa s)")
    ax3.set_ylabel("Stable Droplet Radius (cm) (solid)")
    ax3_2.set_ylabel("Settling Velocity (m/s) (dashed)")
    ax3.set_xscale('log')
    ax3.set_title("Stable Droplet Radius and Settling Velocity vs. Dynamic Viscosity")
    ax3.axvspan(10**-3.5, 10**-1, alpha=0.2, color='red')
    ax3.grid()

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    ax4.plot(etas, v_s_plus_v_c_list, linewidth=2.0)
    ax4_2 = ax4.twinx()
    ax4_2.plot(etas, v_s_minus_v_c_list, linewidth=2.0, linestyle="--")
    ax4.set_xlabel("Dynamic Viscosity (Pa s)")
    ax4.set_ylabel("Downward Convective Regime Settling Velocity (m/s) (solid)")
    ax4_2.set_ylabel("Upward Convective Regime Settling Velocity (m/s) (dashed)")
    ax4.set_title("Settling Velocity of Stable Droplet In Upward and Downward Convective Regimes")
    ax4.axvspan(10 ** -3.5, 10 ** -1, alpha=0.2, color='red')
    ax4.set_xscale('log')
    ax4.grid()

    plt.show()

