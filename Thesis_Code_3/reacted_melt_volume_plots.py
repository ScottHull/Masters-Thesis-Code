from math import pi, log, sqrt
import pandas as pd
import numpy as np
import matplotlib as mpl; mpl.use("Qt5Agg")
import matplotlib.pyplot as plt
from decimal import Decimal

plt.rcParams.update({'font.size': 16})



def rFromWeberTurbulent(density_melt, density_droplet, gravity, Cd=0.2, weber_number=10, surface_tension=1):
    # Rubie et al 2003 suggests We=10 for stable droplet sizes
    # Rubie et al 2003 suggests a surface tension/metal-silicate interface energy of sigma=1 N/m^2
    r = sqrt(weber_number * (((16 * gravity * (density_droplet - density_melt)) / (3 * Cd * surface_tension)) *
                        ((density_droplet - density_melt) / density_melt))**(-1))
    return r

def z_eq(radius_droplet, dynamic_viscosity, settling_velocity, diffusion_coeff, density_melt):
    kinematic_viscosity = dynamic_viscosity / density_melt
    num = ((2 * radius_droplet)**(3/2)) * (kinematic_viscosity**(1/6)) * (settling_velocity**(1/2))
    den = diffusion_coeff**(2/3)
    z = 4.624 * (num/den)
    return z

def turbulentVelocity(gravity, droplet_radius, density_droplet, density_melt, Cd=0.2):
    diameter = 2 * droplet_radius
    v = sqrt(((4 * gravity * diameter) / (3 * Cd)) * ((density_droplet - density_melt) / density_melt))
    return v

def calcDiffusionLength(chem_diffusivity, droplet_radius, settling_velocity):
    l = sqrt((2 * chem_diffusivity * droplet_radius) / settling_velocity)
    return l

def meltLengthWidth(diff_length, droplet_radius):
    length_width = (2 * droplet_radius) + (2 * diff_length)
    return length_width


etas = [10 ** (-3.5), 10 ** (-3.0), 10 ** (-2.5), 10 ** (-2.0), 10 ** (-1.5), 10 ** (-1.0)]# radii = list(np.arange(0.001, 0.02 + 0.001, 0.001))
radii = [0.001, 0.01, 0.1]
gravity_vesta = 0.25
gravity_earth = 9.8
density_droplet = 7800
density_melt = 3750
diffusion_coeff = 10**-8

velocities_vesta = []
z_eqs_vesta = []
velocities_earth = []
z_eqs_earth = []
volume_reacted_melt = []
v_eq_vesta = []
v_eq_earth = []

for eta in etas:
    velocities_vesta_temp = []
    z_eqs_vesta_temp = []
    velocities_earth_temp = []
    z_eqs_earth_temp = []
    v_eq_vesta_temp = []
    v_eq_earth_temp = []
    for i in radii:
        v_earth = turbulentVelocity(gravity=gravity_earth, density_melt=density_melt, density_droplet=density_droplet,
                                    droplet_radius=i)
        z_earth = z_eq(radius_droplet=i, dynamic_viscosity=eta, settling_velocity=v_earth, diffusion_coeff=diffusion_coeff,
                       density_melt=density_melt)
        v_vesta = turbulentVelocity(gravity=gravity_vesta, density_melt=density_melt, density_droplet=density_droplet,
                                    droplet_radius=i)
        z_vesta = z_eq(radius_droplet=i, dynamic_viscosity=eta, settling_velocity=v_vesta, diffusion_coeff=diffusion_coeff,
                       density_melt=density_melt)
        diff_length_vesta = calcDiffusionLength(chem_diffusivity=diffusion_coeff, droplet_radius=i,
                                                settling_velocity=v_vesta)
        diff_length_earth = calcDiffusionLength(chem_diffusivity=diffusion_coeff, droplet_radius=i,
                                                settling_velocity=v_earth)
        print(diff_length_vesta)
        volume_vesta = ((2 * (i + diff_length_vesta)) ** 2) * z_vesta
        volume_earth = ((2 * (i + diff_length_earth)) ** 2) * z_earth
        velocities_vesta_temp.append(v_vesta)
        z_eqs_vesta_temp.append(z_vesta)
        velocities_earth_temp.append(v_earth)
        z_eqs_earth_temp.append(z_earth)
        v_eq_vesta_temp.append(volume_vesta)
        v_eq_earth_temp.append(volume_earth)
    velocities_vesta.append(velocities_vesta_temp)
    z_eqs_vesta.append(z_eqs_vesta_temp)
    velocities_earth.append(velocities_earth_temp)
    z_eqs_earth.append(z_eqs_earth_temp)
    v_eq_vesta.append(v_eq_vesta_temp)
    v_eq_earth.append(v_eq_earth_temp)

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
for index, i in enumerate(etas):
    ax1.plot(radii, v_eq_vesta[index], linewidth=2.0, label="$\eta$={:.2E}".format(Decimal(str(i))))
ax1.grid()
ax1.set_title("Volume of Equilibrated Melt on Vesta for Droplet Radius Per Modeling Interval")
ax1.set_xlabel("Radius (m)")
ax1.set_ylabel("Volume (m$^3$)")
ax1.legend(loc='upper left')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for index, i in enumerate(etas):
    ax2.plot(radii, v_eq_earth[index], linewidth=2.0, label="$\eta$={:.2E}".format(Decimal(str(i))))
ax2.grid()
ax2.set_title("Volume of Equilibrated Melt on Earth for Droplet Radius Per Modeling Interval")
ax2.set_xlabel("Radius (m)")
ax2.set_ylabel("Volume (m$^3$)")
ax2.legend(loc='upper left')

plt.show()