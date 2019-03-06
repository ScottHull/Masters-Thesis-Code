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

def turbulentVelocity(gravity, droplet_radius, density_droplet, density_melt, Cd=0.2):
    diameter = 2 * droplet_radius
    v = sqrt(((4 * gravity * diameter) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
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

def checkCritLaminarViscosity(density_melt, density_droplet, gravity, droplet_radius, f=10):
     eta = sqrt((4 * pi * density_melt * (density_droplet - density_melt) * gravity * (droplet_radius**3)) / (3 * f))
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


def collectCoeffs(pressure, temperature):

    # pressure in GPa
    # temperature in degK

    coeffs = {
        'alpha': 0,
          'beta': 0,
          'chi': 0,
          'delta': 0,
          'epsilon': 0
    }

    coeffs['alpha'] = 1.11
    coeffs['beta'] = -1.18
    coeffs['chi'] = -0.85
    coeffs['delta'] = 1680
    coeffs['epsilon'] = 487

    return coeffs


if __name__ == "__main__":

    radius_vesta = 263 * 1000  # 265 km
    current_vesta_core_radius_min = 107 * 1000
    current_vesta_core_radius_max = 113 * 1000
    depth_increment = 5 * 1000   # 5 km
    densityDroplet = 7800  # kg m3
    densityMelt = 3750  # kg m3
    diffusivity = 10**(-8)
    thermal_expansion = 6 * 10**(-5)
    heat_capacity = 10**3
    thermal_diffusivity = 10**(-6)
    surface_gravity = 0.25

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
        pressure_gradient=[0],  # returns pressures in GPa
    )

    etas = [10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0)]

    droplet_radius_turb_list = []
    droplet_radius_lam_list = []
    f_turb_list = []
    f_lam_list = []
    v_s_turb_list = []
    v_s_lam_list = []
    ra_list = []
    nu_list = []
    v_c_list = []
    eta_crit_turb_list = []
    eta_crit_lam_list = []
    v_s_plus_v_c_turb_list = []
    v_s_minus_v_c_turb_list = []
    v_s_plus_v_c_lam_list = []
    v_s_minus_v_c_lam_list = []

    logD_list = []

    for i in etas:

        droplet_radius_turb = rFromWeberTurbulent(density_droplet=densityDroplet, density_melt=densityMelt,
                                             gravity=surface_gravity)
        droplet_radius_lam = rFromWeberLaminar(dynamic_viscosity=i, density_droplet=densityDroplet,
                                               density_melt=densityMelt, gravity=surface_gravity)
        f_turb = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity,
                          dynamic_viscosity=i, droplet_radius=droplet_radius_turb)
        f_lam = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity,
                          dynamic_viscosity=i, droplet_radius=droplet_radius_lam)
        v_s_turb = turbulentVelocity(gravity=surface_gravity, density_droplet=densityDroplet, density_melt=densityMelt,
                                     droplet_radius=droplet_radius_turb)
        v_s_lam = laminarVelocity(gravity=surface_gravity, density_droplet=densityDroplet, density_melt=densityMelt,
                                     droplet_radius=droplet_radius_lam, dynamic_viscosity=i)
        ra = rayleigh(density_melt=densityMelt, gravity=surface_gravity, dynamic_viscosity=i, depth=radius_vesta,
                      temp_surface=adiabatic[0], temp_depth=adiabatic[-1], thermal_diffusivity=thermal_diffusivity,
                      thermal_expansivity=thermal_expansion)
        nu = nusselt(rayleigh=ra)
        v_c = convectionVelocity(thermal_expansion=thermal_expansion, thermal_diffusivity=thermal_diffusivity,
                                 gravity=surface_gravity, nusselt=nu, temp_surface=adiabatic[0],
                                 temp_depth=adiabatic[-1])
        eta_crit_turb = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet,
                                             droplet_radius=droplet_radius_turb, gravity=surface_gravity)
        eta_crit_lam = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet,
                                             droplet_radius=droplet_radius_lam, gravity=surface_gravity)

        v_s_plus_v_c_turb = v_s_turb + v_c
        v_s_plus_v_c_lam = v_s_lam + v_c
        v_s_minus_v_c_turb = v_s_turb - v_c
        v_s_minus_v_c_lam = v_s_lam - v_c


        droplet_radius_turb_list.append(droplet_radius_turb)
        droplet_radius_lam_list.append(droplet_radius_lam)
        f_turb_list.append(f_turb)
        f_lam_list.append(f_lam)
        v_s_turb_list.append(v_s_turb)
        v_s_lam_list.append(v_s_lam)
        ra_list.append(ra)
        nu_list.append(nu)
        v_c_list.append(v_c)
        eta_crit_turb_list.append(eta_crit_turb)
        eta_crit_lam_list.append(eta_crit_lam)
        v_s_plus_v_c_turb_list.append(v_s_plus_v_c_turb)
        v_s_minus_v_c_turb_list.append(v_s_minus_v_c_turb)
        v_s_plus_v_c_lam_list.append(v_s_plus_v_c_lam)
        v_s_minus_v_c_lam_list.append(v_s_minus_v_c_lam)


    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(etas, [i * 100 for i in droplet_radius_turb_list], linewidth=2.0, label="Turbulent Regime")
    ax1.plot(etas, [i * 100 for i in droplet_radius_lam_list], linewidth=2.0, label="Laminar Regime")
    ax1.set_xlabel("Dynamic Viscosity (Pa s)")
    ax1.set_ylabel("Radius (cm)")
    ax1.set_title("Stable Droplet Radius for Laminar and Turbulent Settling Regimes")
    ax1.grid()
    ax1.legend(loc='lower right')
    ax1.set_xscale('log')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(etas, f_turb_list, linewidth=2.0, label="Turbulent Regime")
    ax2.plot(etas, f_lam_list, linewidth=2.0, label="Laminar Regime")
    ax2.set_xlabel("Dynamic Viscosity (Pa s)")
    ax2.set_ylabel("Friction Coefficient")
    ax2.set_title("Friction for Laminar and Turbulent Settling Regimes")
    ax2.grid()
    ax2.legend(loc='upper right')
    ax2.set_xscale('log')

    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(etas, v_s_turb_list, linewidth=2.0, label="Turbulent Regime")
    ax3.plot(etas, v_s_lam_list, linewidth=2.0, label="Laminar Regime")
    ax3.set_xlabel("Dynamic Viscosity (Pa s)")
    ax3.set_ylabel("Velocity (m/s)")
    ax3.set_title("Settling Velocity for Laminar and Turbulent Settling Regimes")
    ax3.grid()
    ax3.legend(loc='upper right')
    ax3.set_xscale('log')

    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    ax4.plot(etas, v_c_list, linewidth=2.0)
    ax4.set_xlabel("Dynamic Viscosity (Pa s)")
    ax4.set_ylabel("Velocity (m/s)")
    ax4.set_title("Convective Velocities")
    ax4.grid()
    ax4.set_xscale('log')

    fig5 = plt.figure()
    ax5 = fig5.add_subplot(111)
    ax5.axhline(eta_crit_turb_list[0], 0, 1, linewidth=2.0, label="Turbulent Regime", color='blue')
    ax5.plot([i * 100 for i in droplet_radius_lam_list], eta_crit_lam_list, linewidth=2.0, label="Laminar Regime", color='orange')
    ax5.set_xlabel("Droplet Radius (cm)")
    ax5.set_ylabel("Critical Viscosity (Pa s)")
    ax5.set_title("Critical Viscosity for Laminar-Turbulent Transition")
    ax5.grid()
    ax5.legend(loc='upper right')
    ax5.set_yscale('log')


    print("T_surface: {} ... T_depth: {} ... Delta T: {}".format(adiabatic[0], adiabatic[-1],
                                                                 adiabatic[-1] - adiabatic[0]))
    print("P_surface: {} ... P_depth: {} ... Delta P: {}".format(hydrostat[0], hydrostat[-1],
                                                                 hydrostat[-1] - hydrostat[0]))
    print("Droplet Radius (turbulent): MIN: {}   MAX: {}".format(min(droplet_radius_turb_list),
                                                                 max(droplet_radius_turb_list)))
    print("Droplet Radius (laminar): MIN: {}   MAX: {}".format(min(droplet_radius_lam_list),
                                                                 max(droplet_radius_lam_list)))
    print("Friction Coefficient (turbulent): MIN: {}   MAX: {}".format(min(f_turb_list),
                                                                 max(f_turb_list)))
    print("Friction Coefficient (laminar): MIN: {}   MAX: {}".format(min(f_lam_list),
                                                                       max(f_lam_list)))
    print("Settling Velocity (turbulent): MIN: {}   MAX: {}".format(min(v_s_turb_list),
                                                                       max(v_s_turb_list)))
    print("Settling Velocity (laminar): MIN: {}   MAX: {}".format(min(v_s_lam_list),
                                                                    max(v_s_lam_list)))
    print("Critical Viscosity (turbulent): MIN: {}   MAX: {}".format(min(eta_crit_turb_list),
                                                                  max(eta_crit_turb_list)))
    print("Critical Viscosity (laminar): MIN: {}   MAX: {}".format(min(eta_crit_lam_list),
                                                                     max(eta_crit_lam_list)))
    print("Convective Velocity: MIN: {}   MAX: {}".format(min(v_c_list), max(v_c_list)))



    # plt.show()