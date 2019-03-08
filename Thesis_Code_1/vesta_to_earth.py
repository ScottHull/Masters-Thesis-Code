from math import pi, log, sqrt
import pandas as pd
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

def getMeshSizeTimestep(z_eq, t_eq, magma_ocean_depth):
    total_mesh_cells = magma_ocean_depth / z_eq
    time_to_bottom = t_eq * total_mesh_cells
    return total_mesh_cells, time_to_bottom


def linearizeAdiabat(total_mesh_cells, magma_ocean_depth, adiabat):
    min_temp = adiabat[0]
    max_temp = adiabat[-1]
    diff = max_temp - min_temp
    step = diff / total_mesh_cells
    return step

def linearizeHydrostatic(total_mesh_cells, magma_ocean_depth, hydrostatic):
    min_pressure = hydrostatic[0]
    max_pressure = hydrostatic[-1]
    diff = max_pressure - min_pressure
    step = diff / total_mesh_cells
    return step


def calcInitialMolesInDropletEarth(moles_core_of_vesta, droplet_radius, vesta_core_radius):
    droplet_volume = (4/3) * pi * (droplet_radius**3)
    vesta_core_volume = (4/3) * pi * (vesta_core_radius**3)
    droplets_dispered = vesta_core_volume / droplet_volume
    moles_in_droplet = moles_core_of_vesta / droplets_dispered
    return moles_in_droplet


def iterReverseD(obj_concs, cell_concs, index=0, iterReverseDList=[]):

    if index < len(list(obj_concs)):
        obj = list(obj_concs)[index]
        cell_concs_range = list(cell_concs)[0:index + 1]
        avg_cell_concs_range = sum(cell_concs_range) / (len(cell_concs_range))
        # print(cell_concs[index], avg_cell_concs_range)
        avg_D = obj / avg_cell_concs_range
        iterReverseDList.append(avg_D)
        return iterReverseD(obj_concs=obj_concs, cell_concs=cell_concs, index=(index + 1), iterReverseDList=iterReverseDList)
    else:
        return iterReverseDList

def forIterReverseD(obj_concs, cell_concs):
    iterReverseDList = []
    for index in range(len(list(obj_concs))):
        if index + 1 < len(list(obj_concs)):
            obj = list(obj_concs)[index]
            cell_concs_range = list(cell_concs)[0:index + 1]
            avg_cell_concs_range = sum(cell_concs_range) / (len(cell_concs_range))
            # print(cell_concs[index], avg_cell_concs_range)
            avg_D = obj / avg_cell_concs_range
            iterReverseDList.append(avg_D)
        else:
            return iterReverseDList

def calcDiffusionLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def meltLengthWidth(diff_length, droplet_radius):
    length_width = (2 * droplet_radius) + (2 * diff_length)
    return length_width

def recalcConcentration(predicted_d, original_moles_silicate, original_moles_metal, volume_mesh, radius_object):
    volume_obj = (4 / 3) * pi * (radius_object ** 3)

    original_conc_silicate = original_moles_silicate / volume_mesh
    original_conc_metal = original_moles_metal / volume_obj
    concs_mesh = [original_conc_silicate]
    concs_objs = [original_conc_metal]
    moles_mesh = [original_moles_silicate]
    moles_objs = [original_moles_metal]
    verify_D = []
    for index, d in enumerate(list(predicted_d)):
        old_moles_obj = moles_objs[index - 1]
        old_moles_cell = moles_mesh[index - 1]

        adj_matrix = (old_moles_cell /
                      (1 + (3 * volume_mesh * (
                                  (4 * pi * (radius_object ** 3) * d) ** (-1)))))
        adj_object = (old_moles_obj /
                      (1 + (4 * pi * (radius_object ** 3) * d) * (
                                  (3 ** (-1)) * (volume_mesh**(-1)))))
        adj_moles = adj_matrix - adj_object

        # adjust the moles of the element in the object and matrix, respectively
        new_moles_obj = old_moles_obj + adj_moles
        new_moles_cell = old_moles_cell - adj_moles
        new_obj_conc = new_moles_obj / volume_obj
        new_mesh_conc = new_moles_cell / volume_mesh
        check_D = new_obj_conc / new_mesh_conc
        moles_objs.append(new_moles_obj)
        moles_mesh.append(new_moles_cell)
        concs_mesh.append(new_mesh_conc)
        concs_objs.append(new_obj_conc)
        verify_D.append(check_D)

    return concs_mesh[1:], concs_objs[1:], moles_mesh[1:], moles_objs[1:], verify_D[1:]



radius_earth = 520 * 1000  # 1000 km
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
etas = [10 ** (-4.0), 10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0), 10 ** (-0.5), 10 ** (0)]
droplet_radius_samples = [0.001, 0.005, 0.01, 0.015, 0.02]



adiabatic_depths_earth, adiabatic_earth = adiabat(
    current_depth=0 / 1000,  # begin at surface
    current_temperature=2000,  # degrees K
    depth_increment=depth_increment,  # 10 km interval
    depth=radius_earth,  # to depth of 250 km
    gravity=surface_gravity,  # m/s^2, Earth
    thermal_expansion=6 * 10**(-5),
    heat_capacity=10**3,
    depth_gradient=[0],  # returns this list at index=0
    adiabat_gradient=[2000],  # returns this list at index=1
)

hydrostatic_depths_earth, hydrostat_earth = hydrostatic(
    current_depth=0 / 1000,
    depth_increment=depth_increment,
    depth=radius_earth,
    gravity=surface_gravity,
    current_pressure=0,
    density_melt=densityMelt,
    depth_gradient=[0],
    pressure_gradient=[0],
)

adiabatic_depths_vesta, adiabatic_vesta = adiabat(
    current_depth=0 / 1000,  # begin at surface
    current_temperature=2000,  # degrees K
    depth_increment=depth_increment,  # 10 km interval
    depth=265 * 1000,  # to depth of 250 km
    gravity=0.25,  # m/s^2, Vesta
    thermal_expansion=6 * 10**(-5),
    heat_capacity=10**3,
    depth_gradient=[0],  # returns this list at index=0
    adiabat_gradient=[2000],  # returns this list at index=1
)

hydrostatic_depths_vesta, hydrostat_vesta = hydrostatic(
    current_depth=0 / 1000,
    depth_increment=depth_increment,
    depth=265 * 1000,
    gravity=0.25,
    current_pressure=0,
    density_melt=densityMelt,
    depth_gradient=[0],
    pressure_gradient=[0],
)


sample_f_vals = []
for j in droplet_radius_samples:
    temp_sample_list = []
    for i in etas:
        temp_sample_list.append(frictionCoeff(density_melt=densityMelt, density_droplet=densityDroplet, dynamic_viscosity=i, gravity=surface_gravity, droplet_radius=j))
    sample_f_vals.append(temp_sample_list)

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
ax1.set_title("Maximum Stable Droplet Radius vs. Dynamic Viscosity (Earth)")
ax1.set_xlabel("Dynamic Viscosity (Pa s)")
ax1.set_ylabel("Max. Stable Droplet Radius (cm)")
ax1.grid()
ax1.legend(loc='upper left')
ax1.set_xscale('log')

fig2 = plt.figure()
ax2 = fig2.add_subplot(211)
ax2_2 = fig2.add_subplot(212)
ax2.plot(etas, velocity_laminar_list, linewidth=2.0, label="Settling Velocity (Laminar)")
ax2_2.plot(etas, velocity_turbulent_list, linewidth=2.0, label="Setting Velocity (Turbulent)")
ax2.set_title("Stable Settling Velocity vs. Dynamic Viscosity (Earth)")
ax2_2.set_xlabel("Dynamic Viscosity (Pa s)")
ax2.set_ylabel("Settling Velocity (m/s)")
ax2_2.set_ylabel("Settling Velocity (m/s)")
ax1.axvspan(10 ** (-3.5), 10 ** (-1), color='red', alpha=0.2, label="Magma Ocean Viscosity Range")
ax2.grid()
ax2_2.grid()
ax2.legend(loc='upper left')
ax2_2.legend(loc='upper left')
ax2.set_xscale('log')
ax2_2.set_xscale('log')

# recall that viscosity is independent of stable radius and settling velocity in the turbulent regime.
# therefore, does not need to be plotted here.
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot([i * 100 for i in radius_laminar_list], velocity_laminar_list, linewidth=2.0, label="Settling Velocity (Laminar)")
ax3.set_title("Stable Settling Velocity vs. Stable Droplet Radius (Earth)")
ax3.set_xlabel("Stable Droplet Radius (cm)")
ax3.set_ylabel("Settling Velocity (m/s)")
ax3.axvspan(radius_laminar_list[etas.index(10 ** (-3.5))] * 100, radius_laminar_list[etas.index(10 ** (-1.0))] * 100,
            color='red', alpha=0.2, label="Magma Ocean Viscosity Range")
ax3.grid()
ax3.legend(loc='upper left')

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for index, i in enumerate(sample_f_vals):
    ax4.plot(etas, i, linewidth=2.0, label="Droplet Radius: {} cm".format(droplet_radius_samples[index] * 100))
ax4.axvspan(10 ** (-3.5), 10 ** (-1), color='red', alpha=0.2, label="Magma Ocean Viscosity Range")
ax4.set_title("Friction Coefficient (f) vs. Dynamic Viscosity for Possible Droplet Radius (Earth)")
ax4.set_xlabel("Dynamic Viscosity (Pa s)")
ax4.set_ylabel("Friction Coefficient (f)")
ax4.legend(loc='upper right')
ax4.grid()
ax4.set_xscale('log')

fig5 = plt.figure()
ax5 = fig5.add_subplot(211)
ax5_2 = fig5.add_subplot(212)
ax5.plot(adiabatic_depths_earth, adiabatic_earth, linewidth=2.0, label="Earth Magma Ocean Adiabat")
ax5.plot(adiabatic_depths_vesta, adiabatic_vesta, linewidth=2.0, label="Vesta Magma Ocean Adiabat")
ax5_2.plot(hydrostatic_depths_earth, hydrostat_earth, linewidth=2.0, label="Earth Magma Ocean Hydrostatic Pressure Gradient")
ax5_2.plot(hydrostatic_depths_vesta, hydrostat_vesta, linewidth=2.0, label="Vesta Magma Ocean Hydrostatic Pressure Gradient")
ax5.grid()
ax5_2.grid()
ax5.set_title("Adiabatic and Hyrostatic Gradients Over Depth for Earth and Vesta Magma Oceans")
ax5_2.set_xlabel("Depth (km)")
ax5.set_ylabel("Adiabatic Temperature (degK)")
ax5_2.set_ylabel("Hydrostatic Pressure Gradient (GPa)")
ax5.legend(loc='upper right')
ax5_2.legend(loc='upper right')


def cottrellModel(pressure, temperature, fO2, nbo_t=2.6):
    # eg. Cottrell et al. (2009) metal-silicate W partitioning model:
    # log(D) = alpha + beta * (delta IW) + chi * (nbo/t) + delta * (1/T) + epsilon(P/T)


    coeffs = {
        'alpha': 0,
        'beta': 0,
        'chi': 0,
        'delta': 0,
        'epsilon': 0
    }

    if 0.0 <= pressure <= 2:
        coeffs['alpha'] = 1.11
        coeffs['beta'] = -1.18
        coeffs['chi'] = -0.85
        coeffs['delta'] = 1680
        coeffs['epsilon'] = 487
    else:
        coeffs['alpha'] = 1.05
        coeffs['beta'] = -1.10
        coeffs['chi'] = -0.84
        coeffs['delta'] = 3588
        coeffs['epsilon'] = -102

    alpha = coeffs['alpha']
    beta = coeffs['beta']
    chi = coeffs['chi']
    delta = coeffs['delta']
    epsilon = coeffs['epsilon']

    logD = alpha + (beta * fO2) + (chi * nbo_t) + (delta * (1 / temperature)) + (epsilon * (pressure / temperature))
    D = 10 ** logD

    return D

cottrell = []
for index, i in enumerate(adiabatic_depths_earth):
    temp = adiabatic_earth[index]
    pressure = hydrostat_earth[index]
    fO2 = -2.25
    D = cottrellModel(pressure=pressure, temperature=temp, fO2=fO2)
    cottrell.append(D)

earth_droplet_radius = rFromWeberTurbulent(density_melt=densityMelt, density_droplet=densityDroplet,
                                           gravity=surface_gravity)
earth_velocity_turbulent = turbulentVelocity(gravity=surface_gravity, droplet_radius=earth_droplet_radius,
                                           density_droplet=densityDroplet, density_melt=densityMelt)
z_eq_10_minus2 = z_eq(radius_droplet=earth_droplet_radius, dynamic_viscosity=10**(-2), settling_velocity=earth_velocity_turbulent,
                     diffusion_coeff=diffusivity, density_melt=densityMelt)
t_eq_10_minus2 = z_eq_10_minus2 / earth_velocity_turbulent
rounded_z_eq = 40
rounded_t_eq = rounded_z_eq / earth_velocity_turbulent
total_mesh_cells, total_time = getMeshSizeTimestep(z_eq=40, t_eq=rounded_t_eq, magma_ocean_depth=radius_earth)
adiabat_step = linearizeAdiabat(total_mesh_cells=total_mesh_cells, magma_ocean_depth=radius_earth, adiabat=adiabatic_earth)
hydrostatic_step = linearizeHydrostatic(total_mesh_cells=total_mesh_cells, magma_ocean_depth=radius_earth, hydrostatic=hydrostat_earth)
moles_in_droplets = calcInitialMolesInDropletEarth(moles_core_of_vesta=1.0704 * 10**17,
                                                   droplet_radius=earth_droplet_radius, vesta_core_radius=113 * 1000)

print(earth_droplet_radius, earth_velocity_turbulent, z_eq_10_minus2, t_eq_10_minus2)
print(rounded_z_eq, rounded_t_eq)
print(total_mesh_cells, total_time, adiabat_step, hydrostatic_step)
print(moles_in_droplets)

print(hydrostat_earth[-1])

diffusion_length = calcDiffusionLength(chem_diffusivity=diffusivity, droplet_radius=earth_droplet_radius,
                                       settling_velocity=earth_velocity_turbulent)
length_width = meltLengthWidth(diff_length=diffusion_length, droplet_radius=earth_droplet_radius)

df_fO2_225 = pd.read_csv("earth_models/earth_fO2-225.csv")
predicted_D_fO2_225 = df_fO2_225['D']
depths_fO2_225 = df_fO2_225['z-depth']
cell_temps = df_fO2_225['cell_temperature']
cell_pressures = df_fO2_225['cell_pressure']
recalc_concs_mesh_fO2_08, recalc_concs_objs_fO2_08, recalc_moles_mesh_fO2_08, recalc_moles_objs_fO2_08, recalc_verify_D_fO2_08 = \
    recalcConcentration(predicted_d=predicted_D_fO2_225, original_moles_silicate=0.034937625,
                                          original_moles_metal=1.9195542474206624 * 10**(-6), volume_mesh=(rounded_z_eq * (length_width**2)),
                                          radius_object=earth_droplet_radius)
reverse_recalc_concs_fO2_225 = forIterReverseD(obj_concs=recalc_concs_objs_fO2_08, cell_concs=recalc_concs_mesh_fO2_08)

cottrell = []
for index, i in enumerate(list(depths_fO2_225)):
    temp = list(cell_temps)[index]
    pressure = list(cell_pressures)[index]
    fO2 = -2.25
    D = cottrellModel(pressure=pressure, temperature=temp, fO2=fO2)
    cottrell.append(D)

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.plot([i / 1000 for i in list(depths_fO2_225)[:-1]], cottrell[:-1], linestyle="--")
ax6.plot([i / 1000 for i in list(depths_fO2_225)[:-1]], reverse_recalc_concs_fO2_225)


plt.show()




# need to fix the linearized hydrostatic and adiabatic gradients

