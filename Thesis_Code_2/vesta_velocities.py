from math import sqrt, pi, log, exp
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt

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
    dynamic_viscosity = 10**(-2)
    thermal_diffusivity = 10**(-6)
    surface_gravity = surfaceGravity(radius=radius_vesta, melt_density=densityMelt)

    adiabatic_depths, adiabatic = adiabat(
        current_depth=0 / 1000,  # begin at surface
        current_temperature=2000,  # degrees K
        depth_increment=depth_increment,  # 10 km interval
        depth=radius_vesta,  # to depth of 250 km
        gravity=surface_gravity,  # m/s^2, Vesta
        thermal_expansion=thermal_expansion,
        heat_capacity=heat_capacity,
        depth_gradient=[0],  # returns this list at index=0
        adiabat_gradient=[2000],  # returns this list at index=1
    )

    ra = rayleigh(density_melt=densityMelt, thermal_expansivity=thermal_expansion, gravity=surface_gravity,
                  dynamic_viscosity=dynamic_viscosity, depth=radius_vesta, temp_depth=adiabatic[-1],
                  temp_surface=adiabatic[0], thermal_diffusivity=thermal_diffusivity)
    nu = nusselt(rayleigh=ra)
    v_convect = convectionVelocity(thermal_expansion=thermal_expansion, gravity=surface_gravity, nusselt=nu,
                                   temp_depth=adiabatic[-1], temp_surface=adiabatic[0],
                                   thermal_diffusivity=thermal_diffusivity)

    v_dict = {}
    droplet_radii = [0.001, 0.005, 0.01, 0.05]

    for rad in droplet_radii:
        v_list = []
        for index, i in enumerate(adiabatic_depths):
            g = gravity(radius=radius_vesta, melt_density=densityMelt, depth=i * 1000)
            v = velocity(gravAccel=g, density_droplet=densityDroplet, density_melt=densityMelt, radius_droplet=rad)
            v_list.append(v)
        v_dict.update({str(rad): v_list})

    fig = plt.figure()
    ax = fig.add_subplot(111)
    print(v_dict.keys())
    for key in v_dict.keys():
        print(key)
        v_list = v_dict[key]
        ax.plot(adiabatic_depths, v_list, linewidth=2, label="Terminal Sinking Velocity (r = {} cm)".format(float(key) * 100))
    ax.plot(adiabatic_depths, [v_convect for i in adiabatic_depths], linestyle="--", linewidth=2,
            label="Convective Velocity (m/s)")
    ax.set_title("Sinking vs. Convective Velocity for Iron Droplets on Vesta")
    ax.set_ylabel("Velocity (m/s)")
    ax.set_xlabel("Depth Below Surface (km)")
    ax.axvspan((radius_vesta - current_vesta_core_radius_max) / 1000, radius_vesta / 1000, alpha=0.2)
    ax.legend(loc='upper right')
    ax.grid()
    plt.show()

