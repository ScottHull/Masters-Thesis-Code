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

if __name__ == "__main__":

    radius_vesta = 263 * 1000  # 265 km
    current_vesta_core_radius_min = 107 * 1000
    current_vesta_core_radius_max = 113 * 1000
    depth_increment = 5 * 1000  # 5 km
    densityDroplet = 7800  # kg m3
    densityMelt = 3750  # kg m3
    c_0_m = 32700  # ppm
    c_0_s = 42  # ppm
    diffusivity = 10**(-8)
    thermal_expansion = 6 * 10**(-5)
    heat_capacity = 10**3
    thermal_diffusivity = 10**(-6)
    surface_gravity = 0.25

    # etas = [10 ** (-5), 10 ** (-4), 10 ** (-3), 10 ** (-2), 10 ** (-1), 10 ** (-0), 10 ** (1), 10 ** (2), 10 ** (3),
    #         10 ** (4), 10 ** (5)]

    etas = [10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0)]

    v_turb_list = []
    v_lam_list = []
    r_weber_turb_list = []
    r_weber_lam_list = []

    for i in etas:

        r_weber_turb = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet,
                                           gravity=surface_gravity)
        v_turb = turbulentVelocity(gravAccel=surface_gravity, droplet_radius=r_weber_turb, density_droplet=densityDroplet,
                          density_melt=densityMelt)
        v_lam = laminarVelocity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity,
                                droplet_radius=r_weber_lam, dynamic_viscosity=i)

        r_weber_turb_list.append(r_weber_turb)
        r_weber_lam_list.append(r_weber_lam)
        v_turb_list.append(v_turb)
        v_lam_list.append(v_lam)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(etas, [i * 100 for i in r_weber_turb_list], linewidth=2.0, label="Turbulent Settling Regime")
    ax.plot(etas, [i * 100 for i in r_weber_lam_list], linewidth=2.0, label="Laminar Settling Regime")
    ax.grid()
    ax.set_title("Dynamic Viscosity vs. Stable Droplet Radius (We=10)")
    ax.set_xlabel("Dynamic Viscosity (Pa s)")
    ax.set_ylabel("Stable Droplet Radius (cm)")
    ax.set_xscale('log')
    ax.legend(loc='upper left')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(etas, v_turb_list, linewidth=2.0, label="Turbulent Settling Regime", color='blue')
    ax2.plot(etas, v_lam_list, linewidth=2.0, label="Laminar Settling Regime", color='red')
    ax3 = ax2.twinx()
    ax3.plot(etas, [i * 100 for i in r_weber_turb_list], linewidth=2.0, label="Turbulent Settling Regime", color='blue', linestyle="--")
    ax3.plot(etas, [i * 100 for i in r_weber_lam_list], linewidth=2.0, label="Laminar Settling Regime", color='red', linestyle="--")
    ax2.set_title("Dynamic Viscosity vs. Velocity at Stable Droplet Radius (We=10)")
    ax2.set_xlabel("Dynamic Viscosity (Pa s)")
    ax2.set_ylabel("Stable Droplet Velocity (m/s)")
    ax3.set_ylabel("Stable Droplet Radius (cm) (dashed)")
    ax2.set_xscale('log')
    ax2.legend(loc='upper right')
    ax2.grid()


    m_lam = rFromWeberLaminar(dynamic_viscosity=10**(-5), density_droplet=densityDroplet, density_melt=densityMelt,
                                            gravity=surface_gravity)
    l_lam = rFromWeberLaminar(dynamic_viscosity=10 ** (0), density_droplet=densityDroplet, density_melt=densityMelt,
                              gravity=surface_gravity)
    m_turb = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet,
                                           gravity=surface_gravity)
    l_turb = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet,
                                           gravity=surface_gravity)

    print(m_lam * 100, l_lam * 100, m_turb * 100, l_turb * 100)


    plt.show()






