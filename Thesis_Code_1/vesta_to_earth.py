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


def z_eq(radius_droplet, dynamic_viscosity, settling_velocity, diffusion_coeff, density_melt):
    kinematic_viscosity = dynamic_viscosity / density_melt
    num = ((2 * radius_droplet)**(3/2)) * (kinematic_viscosity**(1/6)) * (settling_velocity**(1/2))
    den = diffusion_coeff**(2/3)
    z = 4.624 * (num/den)
    return z




radius_earth = 1000 * 1000  # 265 km
current_vesta_core_radius_min = 107 * 1000
current_vesta_core_radius_max = 113 * 1000
depth_increment = 5 * 1000   # 5 km
densityDroplet = 7800  # kg m3
densityMelt = 3750  # kg m3
diffusivity = 10**(-8)
thermal_expansion = 6 * 10**(-5)
heat_capacity = 10**3
thermal_diffusivity = 10**(-6)
surface_gravity = 9.8
etas = [10 ** (-5.0), 10 ** (-4.5), 10 ** (-4.0), 10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0), 10 ** (-0.5), 10 ** (0)]

# calculate stable droplet radii
radius_laminar_list = []
radius_turbulent_list = []
for i in etas:
    radius_laminar = rFromWeberLaminar(dynamic_viscosity=i, density_droplet=densityDroplet, density_melt=densityMelt,
                                       gravity=surface_gravity)
    radius_turbulent = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet,
                                           gravity=surface_gravity)
    radius_laminar_list.append(radius_laminar)
    radius_turbulent_list.append(radius_turbulent)

# calculate velocities from stable droplet radii
velocity_laminar_list = []
velocity_turbulent_list = []
for index, i in enumerate(etas):
    radius_laminar = radius_laminar_list[index]
    radius_turbulent = radius_turbulent_list[index]
    velocity_laminar = laminarVelocity(density_melt=densityMelt, density_droplet=densityDroplet,
                                       gravity=surface_gravity, droplet_radius=radius_laminar, dynamic_viscosity=i)
    velocity_turbulent = turbulentVelocity(gravity=surface_gravity, droplet_radius=radius_turbulent,
                                           density_droplet=densityDroplet, density_melt=densityMelt)
    velocity_laminar_list.append(velocity_laminar)
    velocity_turbulent_list.append(velocity_turbulent)

# calculate equilibrium time and distance
z_eq_laminar_list = []
z_eq_turbulent_list = []
t_eq_laminar_list = []
t_eq_turbulent_list = []
for index, i in enumerate(etas):
    radius_laminar = radius_laminar_list[index]
    radius_turbulent = radius_turbulent_list[index]
    velocity_laminar = velocity_laminar_list[index]
    velocity_turbulent = velocity_turbulent_list[index]
    z_eq_laminar = z_eq(radius_droplet=radius_laminar, dynamic_viscosity=i, settling_velocity=velocity_laminar,
                     diffusion_coeff=diffusivity, density_melt=densityMelt)
    z_eq_turbulent = z_eq(radius_droplet=radius_turbulent, dynamic_viscosity=i, settling_velocity=velocity_turbulent,
                     diffusion_coeff=diffusivity, density_melt=densityMelt)
    t_eq_laminar = z_eq_laminar / velocity_laminar
    t_eq_turbulent = z_eq_turbulent / velocity_turbulent
    z_eq_laminar_list.append(z_eq_laminar)
    z_eq_turbulent_list.append(z_eq_turbulent)
    t_eq_laminar_list.append(t_eq_laminar)
    t_eq_turbulent_list.append(t_eq_turbulent)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(etas, [i * 100 for i in radius_laminar_list], linewidth=2.0, label="Stable Radius (Laminar) (cm)")
ax1.plot(etas, [i * 100 for i in radius_turbulent_list], linewidth=2.0, label="Stable Radius (Turbulent) (cm)")
ax1.axvspan(10 ** (-3.5), 10 ** (-1), color='red', alpha=0.2, label="Magma Ocean Viscosity Range")
ax1.set_title("Maximum Stable Droplet Radius vs. Dynamic Viscosity")
ax1.set_xlabel("Dynamic Viscosity (Pa s)")
ax1.set_ylabel("Max. Stable Droplet Radius (cm)")
ax1.grid()
ax1.legend(loc='upper left')

fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
ax2_twin = ax2.twiny()
ax2_2 = fig2.add_subplot(212)
ax2_2_twin = ax2_2.twiny()
ax2.plot(etas, velocity_laminar_list, linewidth=2.0, label="Settling Velocity (Laminar)")
ax2_2.plot(etas, velocity_turbulent_list, linewidth=2.0, label="Setting Velocity (Turbulent)")
ax2.set_title("Stable Settling Velocity vs. Dynamic Viscosity")
ax2_2.set_xlabel("Dynamic Viscosity (Pa s)")
ax2.set_ylabel("Settling Velocity (m/s)")
ax2_2.set_ylabel("Settling Velocity (m/s)")
ax2_twin.set_xticks(radius_laminar_list)
ax2_2.set_xlim(ax2.get_xlim())
ax2_2_twin.set_xticks(radius_turbulent_list)
ax1.axvspan(10 ** (-3.5), 10 ** (-1), color='red', alpha=0.2, label="Magma Ocean Viscosity Range")
ax2.grid()
ax2_2.grid()
ax2.legend(loc='upper left')
ax2_2.legend(loc='upper left')


plt.show()