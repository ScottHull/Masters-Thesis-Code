from math import sqrt, pi, log, exp
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


def velocity(gravAccel, radius_droplet, density_droplet, density_melt):
    Cd = 0.2
    v = sqrt(((4 * gravAccel * (2 * radius_droplet)) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def coeff_a(droplet_radius, velocity, F_m, F_s, partition_coeff):
    a = - (velocity / (2 * droplet_radius)) * ((F_m / (F_m + (F_s / partition_coeff))) - 1)
    return a

def coeff_b(droplet_radius, velocity, F_m, F_s, partition_coeff, C_0_s):
    b = (velocity / (2 * droplet_radius)) * ((F_s * C_0_s) / (F_m + (F_s / partition_coeff)))
    return b

def t_eq(a, b, C_t_m, C_0_m):
    # t = (1 / a) * log((0.01 * C_0_m + (b/a)) / (C_0_m + (b/a)))
    t = (1 / a) * log((C_t_m + (b / a)) / (C_0_m + (b / a)))
    return t

def z_eq(t_eq, velocity):
    z = velocity * t_eq
    return z

def F_m(radius_droplet, density_melt, density_droplet, velocity, diffusivity):
    numerator = (droplet_radius**3) * (density_droplet / 3)
    denominator = numerator + ((droplet_radius**2) * density_melt * sqrt((2 * diffusivity * radius_droplet) / velocity))
    f = numerator / denominator
    return f

def F_s(F_m):
    f = 1 - F_m
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

def cottrellPartitioning(temperature, pressure, fO2, nbo_t=2.6):
    coefficients = {
        'alpha': 1.11,
        'beta': -1.18,
        'chi': -0.85,
        'delta': 1680,
        'epsilon': 487,
    }
    alpha = coefficients['alpha']
    beta = coefficients['beta']
    chi = coefficients['chi']
    delta = coefficients['delta']
    epsilon = coefficients['epsilon']

    logD = alpha + (beta * (fO2)) + (chi * (nbo_t)) + (delta * (1 / temperature)) + (epsilon * (pressure / temperature))
    D = exp(logD)

    return D



if __name__ == "__main__":

    radius_vesta = 265 * 1000  # 265 km
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

    ra = rayleigh(density_melt=densityMelt, thermal_expansivity=thermal_expansion, gravity=surface_gravity,
                  dynamic_viscosity=dynamic_viscosity, depth=radius_vesta, temp_depth=adiabatic[-1],
                  temp_surface=adiabatic[0], thermal_diffusivity=thermal_diffusivity)
    nu = nusselt(rayleigh=ra)
    v_convect = convectionVelocity(thermal_expansion=thermal_expansion, gravity=surface_gravity, nusselt=nu,
                                   temp_depth=adiabatic[-1], temp_surface=adiabatic[0],
                                   thermal_diffusivity=thermal_diffusivity)

    print(ra, nu, v_convect)


    # t_list = []
    # z_list = []
    #
    # for index, i in enumerate(adiabatic_depths):
    #     temperature = adiabatic[index]
    #     pressure = hydrostat[index]
    #     D = cottrellPartitioning(temperature=temperature, pressure=pressure, fO2=-1.5)
    #     v = velocity(gravAccel=surface_gravity, density_melt=densityMelt, density_droplet=densityDroplet, radius_droplet=droplet_radius)
    #     f_m = F_m(radius_droplet=droplet_radius, density_droplet=densityDroplet, density_melt=densityMelt, diffusivity=diffusivity, velocity=v)
    #     f_s = F_s(F_m=f_m)
    #     a = coeff_a(droplet_radius=droplet_radius, F_m=f_m, F_s=f_s, partition_coeff=D, velocity=v)
    #     b = coeff_b(droplet_radius=droplet_radius, C_0_s=c_0_s, F_m =f_m, F_s=f_s, partition_coeff=D, velocity=v)
    #     t = t_eq(a=a, b=b, C_0_m=c_0_m, C_t_m=(c_0_m * D))
    #     z = z_eq(t_eq=t, velocity=v)
    #
    #     t_list.append(t)
    #     z_list.append(z)
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(adiabatic_depths, z_list)
    # ax.set_xlabel("Depth (km)")
    # ax.set_ylabel("Equilibrium Time (s)")
    # ax.grid()
    # plt.show()




