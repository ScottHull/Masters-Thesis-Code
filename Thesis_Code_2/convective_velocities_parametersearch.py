from math import pi, log, sqrt
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

def critLaminarViscosity(density_melt, density_droplet, gravity, droplet_radius, f):
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

    # f = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=10**(-2), gravity=surface_gravity, droplet_radius=droplet_radius)
    # r = critLaminarRadius(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=10**(-2), gravity=surface_gravity, f=f)
    # eta = critLaminarViscosity(density_melt=densityMelt, density_droplet=densityDroplet, gravity=surface_gravity, droplet_radius=droplet_radius, f=f)
    # print(f, r, eta, droplet_radius)

    # adiabatic_depths, adiabatic = adiabat(
    #     current_depth=0 / 1000,  # begin at surface
    #     current_temperature=2000,  # degrees K
    #     depth_increment=depth_increment,  # 10 km interval
    #     depth=radius_vesta,  # to depth of 250 km
    #     gravity=surface_gravity,  # m/s^2, Vesta
    #     thermal_expansion=thermal_expansion,
    #     heat_capacity=heat_capacity,
    #     depth_gradient=[0],  # returns this list at index=0
    #     adiabat_gradient=[2000],  # returns this list at index=1
    # )
    #
    # v_dict = {}
    #
    # viscosities = [10**i for i in range(-5, 6)]
    #
    # for i in viscosities:
    #     ra = rayleigh(density_melt=densityMelt, thermal_expansivity=thermal_expansion, gravity=surface_gravity,
    #                   dynamic_viscosity=i, depth=radius_vesta, temp_depth=adiabatic[-1],
    #                   temp_surface=adiabatic[0], thermal_diffusivity=thermal_diffusivity)
    #     nu = nusselt(rayleigh=ra)
    #     v_convect = convectionVelocity(thermal_expansion=thermal_expansion, gravity=surface_gravity, nusselt=nu,
    #                                    temp_depth=adiabatic[-1], temp_surface=adiabatic[0],
    #                                    thermal_diffusivity=thermal_diffusivity)
    #
    #     v_dict.update({str(i): {
    #         'ra': ra,
    #         'nu': nu,
    #         'v_c': v_convect
    #     }})
    #
    # v2_dict = {}
    # droplet_radii = [0.001, 0.005, 0.01, 0.02]
    #
    # for rad in droplet_radii:
    #     v_list = []
    #     for index, i in enumerate(adiabatic_depths):
    #         g = gravity(radius=radius_vesta, melt_density=densityMelt, depth=i * 1000)
    #         v = velocity(gravAccel=g, density_droplet=densityDroplet, density_melt=densityMelt, radius_droplet=rad)
    #         v_list.append(v)
    #     v2_dict.update({str(rad): v_list})
    #
    #
    #
    # fig = plt.figure()
    # ax1 = fig.add_subplot(311)
    # ax2 = fig.add_subplot(312)
    # ax3 = fig.add_subplot(313)
    #
    # ax1.set_xscale('log')
    # ax2.set_xscale('log')
    # ax3.set_xscale('log')
    #
    # ax1.plot(viscosities, [v_dict[str(i)]['ra'] for i in viscosities])
    # ax2.plot(viscosities, [v_dict[str(i)]['nu'] for i in viscosities])
    # ax3.plot(viscosities, [v_dict[str(i)]['v_c'] for i in viscosities])
    #
    # ax3.set_xlabel("Dynamic Viscosity (Pa s)")
    # ax1.set_xticklabels([])
    # ax2.set_xticklabels([])
    #
    # ax1.set_ylabel("Rayleigh #")
    # ax2.set_ylabel("Nusselt #")
    # ax3.set_ylabel("Convective Velocity (m/s)")
    #
    # ax1.set_title("Rayleigh, Nusselt, Convective Velocities for Vesta")
    #
    # ax1.grid()
    # ax2.grid()
    # ax3.grid()
    #
    # fig2 = plt.figure()
    # ax4 = fig2.add_subplot(111)
    # ax4.plot(viscosities, [v_dict[str(i)]['v_c'] for i in viscosities], linewidth=2.0)
    # for key in v2_dict.keys():
    #     ax4.plot(viscosities, [v2_dict[str(key)][0] for i in viscosities], label="Terminal Sinking Velocity (r = {} cm)".format(float(key) * 100), linestyle="--", linewidth=2)
    # ax4.axvline(10**-2, 0, 1, label="Rubie et al 2003", linewidth=2, linestyle="--")
    # ax4.set_title("Convective Velocity on Vesta Compared to Droplet Sinking Velocities")
    # ax4.set_xlabel("Dynamic Viscosity (Pa s)")
    # ax4.set_ylabel("Velocity (m/s)")
    # ax4.legend(loc='upper right')
    # ax4.grid()
    # ax4.set_xscale('log')
    #
    # fig3 = plt.figure()
    # ax5 = fig3.add_subplot(111)
    #
    # for rad in droplet_radii:
    #     for i in viscosities:
    #         f = frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=i, droplet_radius=rad, gravity=surface_gravity)
    #
    #         ax5.plot(i, f, label="Droplet Radius: {} cm".format(rad * 100), linewidth=2.0)
    # ax5.grid()
    # ax5.legend(loc='upper right')
    # ax5.set_title("Friction Coefficient on Vesta (g=0.27 m/s^2) For Variable Droplet Radius and Dynamic Viscosity")
    # ax5.set_xlabel("Dynamic Viscosity (Pa s)")
    # ax5.set_ylabel("Friction Coefficient")
    # ax5.set_xscale('log')
    # plt.show()
