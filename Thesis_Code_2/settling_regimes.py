from math import pi, log, sqrt
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

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

def velocity(gravAccel, radius_droplet, density_droplet, density_melt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * (2 * radius_droplet)) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def surfaceGravity(radius, melt_density):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius**3 * melt_density) / radius**2
    return gSurface

def gravity(radius, melt_density, depth):
    G = 6.67408 * 10**(-11) # m3 kg-1 s-2
    gSurface = (G * (4 / 3) * pi * radius ** 3 * melt_density) / radius**2
    gAccel = gSurface * (1 - (depth/radius))
    return gAccel

def frictionCoeff(density_melt, density_droplet, dynamic_viscosity, gravity, droplet_radius):
    f = (pi / 6) * ((density_droplet - density_melt) / density_melt) * ((density_melt / dynamic_viscosity)**2) * \
        gravity * ((2 * droplet_radius)**3)
    return f

def critLaminarRadius(density_melt, density_droplet, dynamic_viscosity, gravity, f):
    r = (1/2) * ((6 * density_melt * (dynamic_viscosity**2) * f) / (pi * (density_melt**2) * gravity * (density_droplet -
                                                                                               density_melt)))**(1/3)
    return r

def critLaminarViscosity(density_melt, density_droplet, gravity, droplet_radius, f=10):
    eta = ((6 * f * density_melt) / (pi * gravity * (density_melt**2) * ((2 * droplet_radius)**3) * (density_droplet - density_melt)))**(-1/2)
    return eta


if __name__ == "__main__":

    radius_vesta = 265 * 1000  # 265 km
    current_vesta_core_radius_min = 107 * 1000
    current_vesta_core_radius_max = 113 * 1000
    depth_increment = 5 * 1000  # 5 km
    droplet_radius = 0.01 # m
    densityDroplet = 7800  # kg m3
    densityMelt = 3750  # kg m3
    c_0_m = 32700  # ppm
    c_0_s = 42  # ppm
    diffusivity = 10**(-8)
    thermal_expansion = 6 * 10**(-5)
    heat_capacity = 10**3
    # dynamic_viscosity = 10**(-2)
    thermal_diffusivity = 10**(-6)
    surface_gravity = surfaceGravity(radius=radius_vesta, melt_density=densityMelt)

    crit_viscosities = []
    droplet_radii = np.arange(0.000, 0.011, 0.001)
    for i in droplet_radii:
        eta = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity, droplet_radius=i, f=10)
        crit_viscosities.append(eta)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot([i * 100 for i in droplet_radii], crit_viscosities, linewidth=2.0, label="f=10")
    ax.plot([i * 100 for i in droplet_radii],
            [critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity,
                                  droplet_radius=i, f=5) for i in droplet_radii], linewidth=2.0, linestyle="--",
            label="f=5")
    ax.plot([i * 100 for i in droplet_radii],
            [critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity,
                                  droplet_radius=i, f=15) for i in droplet_radii], linewidth=2.0, linestyle="--",
            label="f=15")
    # ax.axvline(0.01 * 100, 0, 1, linestyle="--")
    # ax.axhline(10**-2, 0, 1, linestyle="--")
    ax.fill_between([i * 100 for i in droplet_radii], 0, crit_viscosities, color='red', alpha=0.2)
    ax.fill_between([i * 100 for i in droplet_radii], crit_viscosities, 2, color='blue', alpha=0.2)
    ax.set_title("Critical Dynamic Viscosity for Turbulent (f>10) to Laminar (f<10) Flow Regime Change (f=10)")
    ax.set_xlabel("Droplet Radius (cm)")
    ax.set_ylabel("Critical Dynamic Viscosity (Pa s)")
    ax.grid()
    ax.text(0.8, 0.25, "Turbulent", fontsize=12)
    ax.text(0.4, 1.50, "Laminar", fontsize=12)
    ax.legend(loc='upper left')

    f_static_eta = []
    for i in droplet_radii:
        f = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=10**(-2), gravity=surface_gravity, droplet_radius=i)
        f_static_eta.append(f)

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot([i * 100 for i in droplet_radii], f_static_eta, linewidth=2.0)
    ax2.set_title("Friction Coefficient (f) vs. Droplet Radius for Dynamic Visocisty 10**-2 Pa s")
    ax2.set_xlabel("Droplet Radius (cm)")
    ax2.set_ylabel("Friction Coefficient f")
    ax2.grid()

    etas = [10 ** (-5), 10 ** (-4), 10 ** (-3), 10 ** (-2), 10 ** (-1), 10 ** (-0), 10 ** (1), 10 ** (2), 10 ** (3),
            10 ** (4), 10 ** (5)]

    radii = [0.001, 0.01, 0.1]

    f_varied_rad = {}
    all_f_varied_rad = []

    # for rad in radii:
    #     f_static_rad = []
    #     for i in etas:
    #         f = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=i, gravity=surface_gravity, droplet_radius=rad)
    #         f_static_rad.append(f)
    #         all_f_varied_rad.append(f)
    #     f_varied_rad.update({str(rad): f_static_rad})
    #
    # fig3 = plt.figure()
    # ax3 = fig3.add_subplot(111)
    # for key in f_varied_rad.keys():
    #     ax3.plot(etas, f_varied_rad[key], linewidth=2.0, label="Droplet Radius: {} cm".format(float(key) * 100))
    # ax3.set_title("Friction Coefficient (f) vs. Dynamic Viscosity for Varied Droplet Radii")
    # ax3.axhspan(10, max(all_f_varied_rad), alpha=0.2, color='red')
    # ax3.axhspan(min(all_f_varied_rad), 10, alpha=0.2, color='blue')
    # ax3.text(10**2, 10**7, "Turbulent Settling Regime", fontsize=10)
    # ax3.text(10 ** -2, 10 ** -5, "Laminar Settling Regime", fontsize=10)
    # ax3.set_xlabel("Dynamic Viscosity (Pa s)")
    # ax3.set_ylabel("Friction Coefficient f")
    # ax3.set_xscale('log')
    # ax3.set_yscale('log')
    # ax3.grid()
    # ax3.legend(loc='upper right')
    #
    # m = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity, droplet_radius=0.001, f=10)
    # l = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity, droplet_radius=0.1, f=10)
    #
    # print(m, l)



    plt.show()