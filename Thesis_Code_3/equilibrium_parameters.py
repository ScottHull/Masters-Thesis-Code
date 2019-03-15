from math import pi, log, sqrt
import pandas as pd
import numpy as np
from math import pi, log, exp
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt



def z_eq(droplet_radius, dynamic_viscosity, settling_velocity, chem_diffusivity, density_silicate):
    kinematic_viscosity = dynamic_viscosity / density_silicate
    z_eq_num = ((2 * droplet_radius)**(3/2)) * (kinematic_viscosity**(1/6)) * (settling_velocity**(1/2))
    z_eq_den = chem_diffusivity**(2/3)
    z = 4.624 * (z_eq_num / z_eq_den)
    return z

def calcDiffusionLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def rFromWeberTurbulent(density_melt, density_droplet, gravity, Cd=0.2, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = sqrt(weber_number * (((16 * gravity * (density_droplet - density_melt)) / (3 * Cd * surface_tension)) *
                        ((density_droplet - density_melt) / density_melt))**(-1))
    return r

def turbulentVelocity(gravity, droplet_radius, density_droplet, density_melt, Cd=0.2):
    diameter = 2 * droplet_radius
    v = sqrt(((4 * gravity * diameter) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v



etas = [10 ** (-4.0), 10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0), 10 ** (-0.5), 10 ** (0)]

earth_magma_ocean_depth = 520 * 1000
earth_magma_ocean_depth_increment = 260 * 10
earth_surface_gravity = 9.8

vesta_magma_ocean_depth = 262.8 * 1000
vesta_magma_ocean_depth_increment = 400
vesta_surface_gravity = 0.25

silicate_density = 3750
metal_density = 7800
w_diffusivity = 10**-8
thermal_expansivity = 6 * 10**(-5)
heat_capacity = 10**3

r_stable_vesta = rFromWeberTurbulent(density_melt=silicate_density, density_droplet=metal_density,
                                     gravity=vesta_surface_gravity)
r_stable_earth = rFromWeberTurbulent(density_melt=silicate_density, density_droplet=metal_density,
                                     gravity=earth_surface_gravity)

velocity_vesta = turbulentVelocity(gravity=vesta_surface_gravity, droplet_radius=r_stable_vesta,
                                   density_droplet=metal_density, density_melt=silicate_density)
velocity_earth = turbulentVelocity(gravity=earth_surface_gravity, droplet_radius=r_stable_earth,
                                   density_droplet=metal_density, density_melt=silicate_density)

diff_length_vesta = calcDiffusionLength(chem_diffusivity=w_diffusivity, droplet_radius=r_stable_vesta,
                                        settling_velocity=velocity_vesta)
diff_length_earth = calcDiffusionLength(chem_diffusivity=w_diffusivity, droplet_radius=r_stable_earth,
                                        settling_velocity=velocity_earth)

z_eq_vesta = []
z_eq_earth = []
volume_eq_vesta = []
volume_eq_earth = []
for i in etas:
    z_vesta = z_eq(droplet_radius=r_stable_vesta, dynamic_viscosity=i, settling_velocity=velocity_vesta,
                   chem_diffusivity=w_diffusivity, density_silicate=silicate_density)
    z_earth = z_eq(droplet_radius=r_stable_earth, dynamic_viscosity=i, settling_velocity=velocity_earth,
                   chem_diffusivity=w_diffusivity, density_silicate=silicate_density)
    v_eq_vesta = ((2 * (r_stable_vesta + diff_length_vesta))**2) * z_vesta
    v_eq_earth = ((2 * (r_stable_earth + diff_length_earth)) ** 2) * z_earth
    z_eq_vesta.append(z_vesta)
    z_eq_earth.append(z_earth)
    volume_eq_vesta.append(v_eq_vesta)
    volume_eq_earth.append(v_eq_earth)

fig1 = plt.figure()
ax1_0 = fig1.add_subplot(111)
ax1_1 = fig1.add_subplot(211)
ax1_2 = fig1.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax1_0.spines['top'].set_color('none')
ax1_0.spines['bottom'].set_color('none')
ax1_0.spines['left'].set_color('none')
ax1_0.spines['right'].set_color('none')
ax1_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax1_0.xaxis.labelpad = 20
ax1_0.yaxis.labelpad = 20
ax1_0.set_xlabel("$\eta$ (Pa s))")
ax1_0.set_ylabel("z$_{eq}$ (m)")
ax1_0.set_title("Chemical Equilibrium Sinking Distance for Vestian and Earth Magma Ocean")
ax1_1.plot(etas, z_eq_vesta, linewidth=2.0, color='black', label="Vesta")
ax1_2.plot(etas, z_eq_earth, linewidth=2.0, color='black', label="Earth")
ax1_1.axvspan((10**-3.5), (10**-1), color='red', alpha=0.2, label='Preferred Viscosity Range')
ax1_2.axvspan((10**-3.5), (10**-1), color='red', alpha=0.2, label='Preferred Viscosity Range')
ax1_1.grid()
ax1_2.grid()
ax1_1.legend(loc='lower right')
ax1_2.legend(loc='lower right')
ax1_1.set_xscale('log')
ax1_2.set_xscale('log')


fig2 = plt.figure()
ax2_0 = fig2.add_subplot(111)
ax2_1 = fig2.add_subplot(211)
ax2_2 = fig2.add_subplot(212)
# Turn off axis lines and ticks of the big subplot
ax2_0.spines['top'].set_color('none')
ax2_0.spines['bottom'].set_color('none')
ax2_0.spines['left'].set_color('none')
ax2_0.spines['right'].set_color('none')
ax2_0.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax2_0.xaxis.labelpad = 20
ax2_0.yaxis.labelpad = 20
ax2_0.set_xlabel("$\eta$ (Pa s))")
ax2_0.set_ylabel("Volume (m$^3$)")
ax2_0.set_title("Volume of Modeled Reacted Silicate for Vestian and Earth Magma Ocean")
ax2_1.plot(etas, volume_eq_vesta, linewidth=2.0, color='black', label="Vesta")
ax2_2.plot(etas, volume_eq_earth, linewidth=2.0, color='black', label="Earth")
ax2_1.axvspan((10**-3.5), (10**-1), color='red', alpha=0.2, label='Preferred Viscosity Range')
ax2_2.axvspan((10**-3.5), (10**-1), color='red', alpha=0.2, label='Preferred Viscosity Range')
ax2_1.grid()
ax2_2.grid()
ax2_1.legend(loc='lower right')
ax2_2.legend(loc='lower right')
ax2_1.set_xscale('log')
ax2_2.set_xscale('log')

plt.show()
